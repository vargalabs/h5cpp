/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

#include <string>
#include <stdexcept>
#include <cstring>
#include <zlib.h>

#if !defined(H5CPP_DISABLE_LIBDEFLATE)
#if defined(H5CPP_HAS_LIBDEFLATE)
#include <libdeflate.h>
#elif defined(__has_include)
#if __has_include(<libdeflate.h>)
#include <libdeflate.h>
#define H5CPP_HAS_LIBDEFLATE 1
#endif
#endif
#endif
namespace h5::impl::filter {
	// TODO: figure something out to map c++ filters to C calls? 
	template<class Derived>
	struct filter_t {
		htri_t can_apply(::hid_t dcpl, ::hid_t type, ::hid_t space){
			return static_cast<Derived*>(this)->can_apply_impl(dcpl,type,space);
		}
		herr_t set_local(::hid_t dcpl, ::hid_t type, ::hid_t space){
			return static_cast<Derived*>(this)->set_local_impl(dcpl,type,space);
		}
		size_t callback(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf){
			return 0;
		}
		size_t apply( void* dst, const void* src, size_t size){
			return static_cast<Derived*>(this)->apply(dst,src,size);
		}
		int version;
		unsigned id;
		unsigned encoder_present;
		unsigned decoder_present;
		std::string name;
	};

	using call_t = size_t (*)(void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] );
	inline size_t mock( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		memcpy(dst,src,size);
		return size;
	}
	inline unsigned compression_level(size_t n, const unsigned params[]) {
		return n > 0 ? params[0] : H5CPP_DEFAULT_COMPRESSION;
	}
	inline size_t decompressed_size_hint(size_t input_size, size_t n, const unsigned params[]) {
		return n > 1 ? params[1] : input_size;
	}
	inline size_t deflate_bound(size_t size) {
		return static_cast<size_t>(compressBound(static_cast<uLong>(size)));
	}
	inline size_t zlib_deflate_encode(void* dst, const void* src, size_t size, unsigned level) {
#if defined(H5CPP_HAS_LIBDEFLATE)
		libdeflate_compressor* compressor = libdeflate_alloc_compressor(static_cast<int>(level));
		if (!compressor)
			return 0;
		const size_t nbytes = libdeflate_zlib_compress(
			compressor, src, size, dst, deflate_bound(size));
		libdeflate_free_compressor(compressor);
		return nbytes;
#else
		uLongf nbytes = static_cast<uLongf>(deflate_bound(size));
		const int status = compress2(
			static_cast<Bytef*>(dst), &nbytes,
			static_cast<const Bytef*>(src), static_cast<uLong>(size),
			static_cast<int>(level));
		return status == Z_OK ? static_cast<size_t>(nbytes) : 0;
#endif
	}
	inline size_t zlib_deflate_decode(void* dst, const void* src, size_t compressed_size, size_t output_size) {
#if defined(H5CPP_HAS_LIBDEFLATE)
		libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
		if (!decompressor)
			return 0;
		size_t actual_size = 0;
		const libdeflate_result status = libdeflate_zlib_decompress(
			decompressor, src, compressed_size, dst, output_size, &actual_size);
		libdeflate_free_decompressor(decompressor);
		return status == LIBDEFLATE_SUCCESS ? actual_size : 0;
#else
		uLongf actual_size = static_cast<uLongf>(output_size);
		const int status = uncompress(
			static_cast<Bytef*>(dst), &actual_size,
			static_cast<const Bytef*>(src), static_cast<uLong>(compressed_size));
		return status == Z_OK ? static_cast<size_t>(actual_size) : 0;
#endif
	}
	inline size_t deflate( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		if (flags & H5Z_FLAG_REVERSE)
			return zlib_deflate_decode(dst, src, size, decompressed_size_hint(size, n, params));
		return zlib_deflate_encode(dst, src, size, compression_level(n, params));
	}
	inline size_t scaleoffset( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t gzip( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		return deflate(dst, src, size, flags, n, params);
	}
	inline size_t szip( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t nbit( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t fletcher32( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[]){
		memcpy(dst,src,size);
		return size;
	}

	inline size_t add( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t shuffle( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t jpeg( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t disperse( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		memcpy(dst,src,size);
		return size;
	}
	inline size_t error( void* dst, const void* src, size_t size, unsigned flags, size_t n, const unsigned params[] ){
		throw std::runtime_error("invalid filter");
		return size;
	}
	inline call_t get_callback( H5Z_filter_t filter_id ){

		switch( filter_id ){
			case H5Z_FILTER_DEFLATE: return filter::gzip;
			case H5Z_FILTER_SHUFFLE: return filter::shuffle;
			case H5Z_FILTER_FLETCHER32: return filter::fletcher32;
			case H5Z_FILTER_SZIP: return filter::szip;
			case H5Z_FILTER_NBIT: return filter::nbit;
			case H5Z_FILTER_SCALEOFFSET: return filter::scaleoffset;
			default:
					return filter::error;
		}
	}
}
