/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <hdf5.h>

#if defined(DLIB_MATRIx_HEADER) || defined(H5CPP_USE_DLIB)
namespace h5::dlib {
// dlib template:
// const dlib::matrix<short int, 0, 0, dlib::memory_manager_stateless_kernel_1<char>, dlib::row_major_layout>&
		template<class T> using rowmat = ::dlib::matrix<T, 0, 0,
			::dlib::memory_manager_stateless_kernel_1<char>,
			::dlib::row_major_layout>;
		template <class Object, class T = typename impl::decay<Object>::type>
			using is_supported = std::bool_constant<std::is_same_v<Object,h5::dlib::rowmat<T>>>;
}
namespace h5::meta {
		template <class T> struct is_contiguous<h5::dlib::rowmat<T>> : std::true_type {};
}

namespace h5::impl {
	// 1.) object -> H5T_xxx
	template <class T> struct decay<h5::dlib::rowmat<T>>{ using type = T; };

	// get read access to datastaore
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::dlib::is_supported<Object>::value,
	const T*> data(const Object& ref ){
			return &ref(0,0);
	}
	// read write access
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::dlib::is_supported<Object>::value,
	T*> data( Object& ref ){
			return &ref(0,0);
	}
	template<class T> struct rank<h5::dlib::rowmat<T>> : public std::integral_constant<size_t,2>{};



	template<class T>
	inline std::array<size_t,2> size( const dlib::rowmat<T>& ref ){
				return { (hsize_t)ref.nc(),(hsize_t)ref.nr()};
	}



	template <class T> struct get<h5::dlib::rowmat<T>> {
		static inline h5::dlib::rowmat<T> ctor( std::array<size_t,2> dims ){
			return h5::dlib::rowmat<T>( dims[1], dims[0] );
	}};
}
#endif
