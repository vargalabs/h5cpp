/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <vector>

#if defined(XTENSOR_HPP) || defined(H5CPP_USE_XTENSOR)

namespace h5::impl {
	// 1.) object -> H5T_xxx
	template <class T> struct decay<xt::xarray<T>>{ using type = T; };
	template <class T, size_t N> struct decay<xt::xtensor<T, N>>{ using type = T; };

	// get read access to datastore
	template <class T> inline
	const T* data(const xt::xarray<T>& ref ){
			return ref.data();
	}
	template <class T> inline
	T* data(xt::xarray<T>& ref ){
			return ref.data();
	}
	template <class T, size_t N> inline
	const T* data(const xt::xtensor<T, N>& ref ){
			return ref.data();
	}
	template <class T, size_t N> inline
	T* data(xt::xtensor<T, N>& ref ){
			return ref.data();
	}

	// rank
	template<class T, size_t N> struct rank<xt::xtensor<T, N>> : public std::integral_constant<size_t, N>{};

	// determine rank and dimensions
	template <class T> inline h5::count_t size( const xt::xarray<T>& ref ){
		h5::count_t count;
		count.rank = static_cast<int>(ref.dimension());
		for(size_t i = 0; i < ref.dimension(); ++i)
			count[i] = ref.shape()[i];
		return count;
	}
	template <class T, size_t N> inline std::array<size_t,N> size( const xt::xtensor<T, N>& ref ){
		std::array<size_t,N> dims;
		for(size_t i = 0; i < N; ++i)
			dims[i] = ref.shape()[i];
		return dims;
	}

	// CTOR-s
	template <class T> struct get<xt::xarray<T>> {
		static inline xt::xarray<T> ctor( h5::count_t count ){
			std::vector<size_t> shape(count.rank);
			for(int i = 0; i < count.rank; ++i)
				shape[i] = count[i];
			return xt::xarray<T>::from_shape(shape);
	}};
	template <class T, size_t N> struct get<xt::xtensor<T, N>> {
		static inline xt::xtensor<T, N> ctor( std::array<size_t,N> dims ){
			return xt::xtensor<T, N>::from_shape(dims);
	}};
}

namespace h5::meta {
    template <class T> struct is_contiguous<xt::xarray<T>> : std::true_type {};
    template <class T, size_t N> struct is_contiguous<xt::xtensor<T, N>> : std::true_type {};
}

#endif
