/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

#define H5CPP_DAPL_HIGH_THROUGHPUT "h5cpp_dapl_highthroughput"

namespace h5::impl {
	inline herr_t dapl_pipeline_close( const char *name, size_t size, void *ptr ){
		delete *static_cast< impl::pipeline_t<impl::basic_pipeline_t>**>( ptr );
		return 0;
	}
	inline ::herr_t dapl_pipeline_set(::hid_t dapl ) {
		//TODO: see why H5CPP pipeline crashes with 1.12.x  
#if H5_VERSION_LE(1,10,6) 
		// ignore if already set 
		if( H5Pexist(dapl, H5CPP_DAPL_HIGH_THROUGHPUT) ) return 0;
		using pipeline = impl::pipeline_t<impl::basic_pipeline_t>;
	 	pipeline* ptr = new pipeline();
		return H5Pinsert2(dapl, H5CPP_DAPL_HIGH_THROUGHPUT, sizeof( pipeline* ), &ptr,
				nullptr, nullptr, nullptr, nullptr, nullptr, dapl_pipeline_close);
#else
			return 0;
#endif
	}
}

namespace h5 {
// DATA ACCESS PROPERTY LISTS
#if H5_VERSION_GE(1,10,0)
	using efile_prefix 		   = impl::dapl_call< impl::dapl_args<hid_t,const char*>,H5Pset_efile_prefix>;
	using virtual_view         = impl::dapl_call< impl::dapl_args<hid_t,H5D_vds_view_t>,H5Pset_virtual_view>;
	using virtual_printf_gap   = impl::dapl_call< impl::dapl_args<hid_t,hsize_t>,H5Pset_virtual_printf_gap>;
#endif
	using chunk_cache          = impl::dapl_call< impl::dapl_args<hid_t,size_t, size_t, double>,H5Pset_chunk_cache>;
	//using num_threads  	   = impl::dapl_call< impl::dapl_args<hid_t, unsigned char>,impl::dapl_threads>;
	namespace flag {
		using high_throughput  = impl::dapl_call< impl::dapl_args<hid_t>,impl::dapl_pipeline_set>;
	}
	const static flag::high_throughput high_throughput;

	namespace impl {
		// Heap-allocated, intentionally never deleted.  HDF5's atexit handler
		// closes the underlying property list.  Avoids the static-destruction-
		// order race with HDF5's library shutdown that manifests as the
		// "HDF5: infinite loop closing library" cascade on Windows MSVC.
		inline const h5::dapl_t& _dapl_singleton() {
			static h5::dapl_t* p = new h5::dapl_t(static_cast<h5::dapl_t>(H5Pcreate(H5P_DATASET_ACCESS)));
			return *p;
		}
		inline const h5::dapl_t& _default_dapl_singleton() {
			static h5::dapl_t* p = new h5::dapl_t(static_cast<h5::dapl_t>(H5Pcreate(H5P_DATASET_ACCESS)));
			return *p;
		}
	}
	inline const h5::dapl_t& dapl = impl::_dapl_singleton();
	//const static h5::dapl_t default_dapl = high_throughput;
	inline const h5::dapl_t& default_dapl = impl::_default_dapl_singleton();
}

