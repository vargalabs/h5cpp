/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */

#pragma once
#include "H5capi.hpp"
#include "H5Tmeta.hpp"
#include "H5cout.hpp"
#include "H5Zpipeline_threaded.hpp"
#include <memory>
#include <string>
#include <variant>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <ostream>

namespace h5 {
	struct pt_t;
}
std::ostream& operator<<(std::ostream& os, const h5::pt_t& pt);

namespace h5::impl {
    // pt_t::pipeline selects between the synchronous basic pipeline (default,
    // bytewise-identical to pre-241 behavior) and the parallel threaded pipeline.
    // Both alternatives are indirect-owned through unique_ptr so that the
    // variant remains move-assignable regardless of the underlying pipeline's
    // move semantics (threaded_pipeline_t deletes moves because it owns
    // std::jthread workers and atomics; basic_pipeline_t inherits a manually-
    // written move-assign from pipeline_t<Derived> that suppresses the
    // implicit move-ctor needed by variant assignment).
    using pt_pipeline_t = std::variant<
        std::unique_ptr<impl::basic_pipeline_t>,
        std::unique_ptr<impl::threaded_pipeline_t>
    >;
}

// packet table template specialization with inheritance
namespace h5 {
	struct pt_t {
		pt_t();
		pt_t( const h5::ds_t& handle ); // conversion ctor — synchronous pipeline
		pt_t( const h5::ds_t& handle, h5::filter::threads workers ); // threaded pipeline
		// deep copy with own cache memory — always uses the synchronous pipeline,
		// since the threaded pipeline owns workers that cannot be duplicated.
		pt_t( const h5::pt_t& pt ) : h5::pt_t(pt.ds) {
		};
		~pt_t();

		pt_t& operator=( h5::pt_t&& pt ){
            // prevent self assign
            if (this == &pt) return *this;
            if(H5Iis_valid(this->ds)){ // flush and close dataset
                this->flush();
                free(this->fill_value);
            }

            this->ds = std::move(pt.ds);
            this->dxpl = std::move(pt.dxpl);
            this->pipeline = std::move(pt.pipeline);

            this->block_size = pt.block_size;
            this->element_size = pt.element_size;
            this->N = pt.N; this->n = pt.n; this->rank = pt.rank;
            this->ptr = pt.ptr;  this->fill_value = pt.fill_value;

            pt.ptr = nullptr; pt.fill_value = nullptr;
            pt.N=0; pt.n=0; pt.rank=0;
		    for(hsize_t i=0; i<rank; i++){
			    this->offset[i] = pt.offset[i];
                this->current_dims[i] = pt.current_dims[i];
                this->chunk_dims[i] = pt.chunk_dims[i];
                this->count[i] = pt.count[i];
            }
			return *this;
		}
		friend std::ostream& ::operator<<(std::ostream &os, const h5::pt_t& pt);
		template<class T>
		friend void append( h5::pt_t& ds, const T& ref);
		friend void flush(h5::pt_t&);
		// resets the packet-table dimension tracker so the same pt_t can be reused
		// for a fresh logical session (e.g. start-of-day re-init in streaming sinks).
		void reset();
		private:
		void init(const h5::ds_t& ds_);
		void flush();

		template<class T> inline std::enable_if_t<h5::meta::is_scalar<T>::value,
		void> append( const T* ptr );
		template<class T> inline std::enable_if_t< h5::meta::is_scalar<T>::value && !std::is_pointer_v<T>,
		void> append( const T& ref );
		template<class T> inline std::enable_if_t< !h5::meta::is_scalar<T>::value,
		void> append( const T& ref );
		void append( const std::string& ref );

		// Resolve the variant to a reference to the live pipeline_t<Derived> CRTP
		// instance, regardless of which alternative is active. Both alternatives
		// expose the same set of fields (chunk0, block_size, ds, dxpl, ...) inherited
		// from pipeline_t<Derived>, so the visiting lambda can use the result
		// duck-typed across alternatives.
		template <typename F>
		decltype(auto) visit_pipeline(F&& f) {
			return std::visit([&](auto& p) -> decltype(auto) {
				return std::forward<F>(f)(*p);   // both alternatives are unique_ptr
			}, pipeline);
		}

		impl::pt_pipeline_t pipeline;
		h5::dxpl_t dxpl;
		h5::ds_t ds;
		h5::dt_t<void> dt;

		hsize_t offset[H5CPP_MAX_RANK], current_dims[H5CPP_MAX_RANK],
			chunk_dims[H5CPP_MAX_RANK], count[H5CPP_MAX_RANK];
		size_t block_size,element_size,N,n,rank;
		void *ptr, *fill_value;
	};
}


/* initialized to invalid state
 * */
inline h5::pt_t::pt_t() :
	pipeline{std::make_unique<impl::basic_pipeline_t>()},
	dxpl{H5Pcreate(H5P_DATASET_XFER)},ds{H5I_UNINIT},n{0},fill_value{nullptr}{
		for(hsize_t i=0; i<H5CPP_MAX_RANK; i++ )
			count[i] = 1, offset[i] = 0;
	}

// conversion ctor — synchronous pipeline (default)
inline
h5::pt_t::pt_t( const h5::ds_t& handle ) : pt_t() {
	/*default ctor has an invalid state -- skip initialization */
	if( !is_valid(handle) ) return;
	init(handle);
}

// conversion ctor — threaded pipeline with N compression workers
inline
h5::pt_t::pt_t( const h5::ds_t& handle, h5::filter::threads workers ) : pt_t() {
	if( !is_valid(handle) ) return;
	auto threaded = std::make_unique<impl::threaded_pipeline_t>();
	threaded->set_worker_count(workers.n);
	pipeline.emplace<std::unique_ptr<impl::threaded_pipeline_t>>(std::move(threaded));
	init(handle);
}

inline
h5::pt_t::~pt_t(){
	/*default ctor has an invalid state -- skip flushing cache */
	if( !h5::is_valid( ds ) )
		return;
	this->flush();
	free(this->fill_value);
}

inline
void h5::pt_t::init( const h5::ds_t& handle ){
	try {
		// Re-open with zero HDF5 chunk cache: the default DAPL allocates ~1MB per H5Dopen2
		// call which accumulates in malloc arenas when pt_t is used in loops, even though
		// the memory is logically freed on close.
		hid_t raw = static_cast<hid_t>(handle);
		hid_t fid = H5Iget_file_id(raw);
		ssize_t nlen = H5Iget_name(raw, nullptr, 0);
		std::vector<char> dname(static_cast<size_t>(nlen) + 1);
		H5Iget_name(raw, dname.data(), dname.size());
		hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
		H5Pset_chunk_cache(dapl, 0, 0, H5D_CHUNK_CACHE_W0_DEFAULT);
		ds = h5::ds_t{H5Dopen2(fid, dname.data(), dapl)};
		H5Pclose(dapl);
		H5Fclose(fid);
		dt = h5::dt_t<void>{H5Dget_type(static_cast<hid_t>(ds))};
		h5::sp_t file_space = h5::get_space( handle );
		rank = h5::get_simple_extent_dims( file_space, current_dims, nullptr );

		h5::dcpl_t dcpl = h5::get_dcpl( ds );
		h5::dt_t<void*> type = h5::get_type<void*>( ds );
		hsize_t size = h5::get_size( type );
		this->fill_value = h5::get_fill_value(dcpl, type, size);
		visit_pipeline([&](auto& p) {
			p.set_cache(dcpl, size);
			this->ptr = p.chunk0;
			this->block_size = p.block_size;
			this->element_size = p.element_size;
			this->N = p.n;
			p.ds = ds; p.dxpl = dxpl;
		});
		h5::get_chunk_dims( dcpl, chunk_dims );
		for(hsize_t i=1; i<rank; i++)
			current_dims[i] = chunk_dims[i];
	} catch ( ... ){
		throw h5::error::io::packet_table::misc( H5CPP_ERROR_MSG("CTOR: unable to create handle from dataset..."));
	}
}
template<class T> inline std::enable_if_t< h5::meta::is_scalar<T>::value,
void> h5::pt_t::append( const T* ptr ) try {
	//PTR: write directly chunk size from provided buffer/ptr
	*offset = *current_dims;
	*current_dims += *chunk_dims;
	h5::set_extent(ds, current_dims);
	visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, ptr); });
} catch( const std::runtime_error& err ){
	throw h5::error::io::dataset::append( err.what() );
}
inline void h5::pt_t::append( const std::string& ref ) {
	static_cast<const char**>( ptr )[n++] = ref.data();
	if( n != N ) return;

	*offset = *current_dims;
	*current_dims += *chunk_dims;
	h5::set_extent(ds, current_dims);

	hsize_t block = 1, count = n;
	h5::sp_t mem_space{H5Screate_simple(static_cast<int>(rank), &count, nullptr )};
	h5::sp_t file_space{H5Dget_space( static_cast<::hid_t>(ds) )};
	h5::select_all( mem_space );
	H5Sselect_hyperslab( static_cast<hid_t>(file_space), H5S_SELECT_SET, offset, nullptr, &block, &count);
	
	H5Dwrite( static_cast<hid_t>( ds ), 
		dt, mem_space, file_space, static_cast<hid_t>(dxpl), ptr);
	n = 0;
}
template <>
inline void h5::pt_t::append( const char* ref ) {
	static_cast<const char**>( ptr )[n++] = ref;
	if( n != N ) return;

	*offset = *current_dims;
	*current_dims += *chunk_dims;
	h5::set_extent(ds, current_dims);

	hsize_t block = 1, count = n;
	h5::sp_t mem_space{H5Screate_simple(static_cast<int>(rank), &count, nullptr )};
	h5::sp_t file_space{H5Dget_space( static_cast<::hid_t>(ds) )};
	h5::select_all( mem_space );
	H5Sselect_hyperslab( static_cast<hid_t>(file_space), H5S_SELECT_SET, offset, NULL, &block, &count);

	H5Dwrite( static_cast<hid_t>( ds ),
		dt, mem_space, file_space, static_cast<hid_t>(dxpl), ptr);
	n = 0;
}


template<class T> inline std::enable_if_t< h5::meta::is_scalar<T>::value && !std::is_pointer_v<T>,
void> h5::pt_t::append( const T& ref ) try {
//SCALAR: store inbound data directly in pipeline cache
	static_cast<T*>( ptr )[n++] = ref;
	if( n != N ) return;

	n = 0;
	*offset = *current_dims;
	*current_dims += *chunk_dims;
	h5::set_extent(ds, current_dims);
	visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, ptr); });
} catch( const std::runtime_error& err ){
	throw h5::error::io::dataset::append( err.what() );
}

template<class T> inline std::enable_if_t< !h5::meta::is_scalar<T>::value,
void> h5::pt_t::append( const T& ref ) try {
	// Iterator-only containers (list, forward_list, deque, set, …) have no .data().
	// Delegate element-by-element to the scalar append overload.
	using access = h5::meta::access_traits_t<T>;
	if constexpr (access::kind == h5::meta::access_t::iterators) {
		for (const auto& elem : ref)
			this->append(elem);
	} else {

	auto dims = meta::size( ref );

	*offset = *current_dims;
	*current_dims += 1;
	h5::set_extent(ds, current_dims);
	auto ptr_ = meta::data( ref );
	auto dims_ = meta::size( ref );

	switch( dims_.size() ){
		case 1: // vector
			if( dims[0] * element_size == block_size )
				visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, (void*) ptr_ ); });
			else throw h5::error::io::packet_table::write(
					H5CPP_ERROR_MSG("dimension mismatch: "
						+ std::to_string( dims[0] * element_size) + " != " + std::to_string(block_size) ));
			break;
		case 2: //matrix
			if( dims[0] * dims[1] * element_size == block_size )
				visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, (void*) ptr_ ); });
			else throw h5::error::io::packet_table::write(
					H5CPP_ERROR_MSG("dimension mismatch: "
						+ std::to_string( dims[0] * dims[1] * element_size) + " != " + std::to_string(block_size) ));
			break;
		case 3: // cube
			if( dims[0] * dims[1] * dims[2] * element_size == block_size )
				visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, (void*) ptr_ ); });
			else throw h5::error::io::packet_table::write(
					H5CPP_ERROR_MSG("dimension mismatch: "
						+ std::to_string( dims[0] * dims[1] * dims[2] * element_size) + " != " + std::to_string(block_size) ));
			;break;
		default:
			throw h5::error::io::packet_table::misc( H5CPP_ERROR_MSG("objects with rank > 2 are not supported... "));
	}
	} // end else (non-iterator path)
} catch( const std::runtime_error& err ){
	throw h5::error::io::dataset::append( err.what() );
}

inline
void h5::pt_t::flush(){
	if( n == 0 ) return;
	*offset = *current_dims;
	*current_dims += *chunk_dims;
	h5::set_extent(ds, current_dims);

	if( H5Tis_variable_str(this->dt)) {
		hsize_t block = 1, count = n;
	 	h5::sp_t mem_space{H5Screate_simple(static_cast<int>(rank), &count, nullptr )};
		h5::sp_t file_space{H5Dget_space( static_cast<::hid_t>(ds) )};
		h5::select_all( mem_space );
		H5Sselect_hyperslab( static_cast<hid_t>(file_space), H5S_SELECT_SET, offset, nullptr, &block, &count);

		H5Dwrite( static_cast<hid_t>( ds ),
			dt, mem_space, file_space, static_cast<hid_t>(dxpl), ptr);
	} else {
		// the remainder of last chunk must be set to fill_value; arbitrary type size supported
		for(hsize_t i=0; i<(N-n); i++)
			for(size_t j=0; j < element_size; j++)
				static_cast<char*>( ptr )[(n + i) * element_size + j] = static_cast<char*>( fill_value )[ j ];
    	visit_pipeline([&](auto& p){ p.write_chunk(offset, block_size, ptr); });
	}
	n = 0;
}

inline void h5::pt_t::reset() {
	std::memset(current_dims, 0, H5CPP_MAX_RANK * sizeof(hsize_t));
}

namespace h5 {
	/** @ingroup io-append
	 * @brief extends HDF5 dataset along the first/slowest growing dimension, then writes passed object to the newly created space
	 * @param pt packet_table descriptor
	 * @param ref T type const reference to object appended
	 * @tparam T dimensions must match the dimension of HDF5 space upto rank-1
	 */

	template<class T> inline
	void append( h5::pt_t& pt, const T& ref){
		pt.append( ref );
	}

	inline void flush(h5::pt_t& pt) try {
		pt.flush();
        //TODO: find better mechanism for deprecating code: #pragma message("not implemented: do not call pt_t::flush() ...")
		// for now
	} catch ( const std::runtime_error& e){
		throw h5::error::io::dataset::close( e.what() );
	}
	/** @ingroup io-append
	 * @brief zeros the packet-table dimension tracker so the same pt_t can be
	 * reused for a fresh logical session. Does not shrink the underlying
	 * dataset on disk; the caller is responsible for any HDF5-level cleanup.
	 * @param pt packet_table descriptor
	 */
	inline void reset(h5::pt_t& pt) try {
		pt.reset();
	} catch ( const std::runtime_error& e){
		throw h5::error::io::dataset::write( e.what() );
	}
}

inline std::ostream& operator<<(std::ostream &os, const h5::pt_t& pt) {
    os << std::dec;
	os <<"packet table:\n"
		 "------------------------------------------\n";
    if( !H5Iis_valid(pt.ds)) {
        os << "ds: H5I_UNINIT" <<std::endl;
        return os;
    }
	std::vector<size_t> current_dims(pt.current_dims, pt.current_dims + pt.rank);
	std::vector<size_t> chunk_dims(pt.chunk_dims, pt.chunk_dims + pt.rank);
	std::vector<size_t> count(pt.count, pt.count + pt.rank);
	std::vector<size_t> offset(pt.offset, pt.offset + pt.rank);

	os << "rank: " << pt.rank << " N:" << pt.N <<" n:" << pt.n << "\n";
	os << "element size: " << pt.element_size << " block size: " << pt.block_size << "\n";
	os << "current dims: " << current_dims << std::endl;
	os << "chunk dims: " << chunk_dims << std::endl;
	os << "offset : " <<  offset <<  " count : " << count << std::endl;
	os << "ds: "<< static_cast<hid_t>( pt.ds ) <<" dxpl: "<< static_cast<hid_t>( pt.dxpl ) << std::endl;
	os << "fill value: " << std::hex << pt.fill_value << " buffer: " << pt.ptr;
	os << "\n\n";
    os << std::dec;
	return os;
}
