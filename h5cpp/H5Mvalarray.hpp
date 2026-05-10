/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

#include <valarray>

namespace h5::valarray {
		template<class T> using array = ::std::valarray<T>;
		template <class Object, class T = typename impl::decay<Object>::type>
			using is_supported = std::bool_constant<std::is_same_v<Object,h5::valarray::array<T>>>;
}

namespace h5::meta {
    template <class T> struct is_contiguous<h5::valarray::array<T>> : std::true_type {};
}

namespace h5::impl {
	// 1.) object -> H5T_xxx
	template <class T> struct decay<h5::valarray::array<T>>{ using type = T; };
	// get read access to datastaore
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::valarray::is_supported<Object>::value,
	const T*> data(const Object& ref ){
			return std::begin(ref);
	}
	// read write access
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::valarray::is_supported<Object>::value,
	T*> data( Object& ref ){
			return std::begin(ref);
	}
	template<class T> struct rank<h5::valarray::array<T>> : public std::integral_constant<size_t,1>{};
	template <class T> inline std::array<size_t,1> size( const h5::valarray::array<T>& ref ){ return { (hsize_t)ref.size() };}
	template <class T> struct get<h5::valarray::array<T>> {
		static inline h5::valarray::array<T> ctor( std::array<size_t,1> dims ){
			return h5::valarray::array<T>( dims[0] );
	}};
}
