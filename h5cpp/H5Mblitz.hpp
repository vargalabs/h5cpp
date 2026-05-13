/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#ifndef  H5CPP_BLITZ_HPP 
#define  H5CPP_BLITZ_HPP

#if defined(BZ_ARRAY_H) || defined(H5CPP_USE_BLITZ)
namespace h5::blitz {
		template<class T, int N> using array = ::blitz::Array<T,N>;

		template <class Object, class T = typename impl::decay<Object>::type, int N = Object::rank_>
		using is_supported = std::bool_constant<std::is_same_v<Object,h5::blitz::array<T,N>>>;
}

namespace h5::meta {
		template <class T, int N> struct is_contiguous<h5::blitz::array<T,N>> : std::true_type {};
}

namespace h5::impl {
	// 1.) object -> H5T_xxx
	template <class T, int N> struct decay<h5::blitz::array<T,N>>{ using type = T; };

	// get read access to datastore
	template <class Object, class T = typename impl::decay<Object>::type, int N = Object::rank_> inline
	std::enable_if_t< h5::blitz::is_supported<Object>::value,
	const T*> data(const Object& ref ){
			return ref.data();
	}
	// read write access
	template <class Object, class T = typename impl::decay<Object>::type, int N = Object::rank_> inline
	std::enable_if_t< h5::blitz::is_supported<Object>::value,
	T*> data( Object& ref ){
			return ref.data();
	}
	// rank
	template<class T, int N> struct rank<h5::blitz::array<T,N>> : public std::integral_constant<size_t,N>{};
	// determine rank and dimensions
	template <class T, int N>
	inline std::array<size_t,N> size( const h5::blitz::array<T,N>& ref ){
		std::array<size_t,N> dims;
		for(int i=0; i<N; ++i) dims[i] = ref.extent(i);
		return dims;
	}

	template <class T, int N> struct get<h5::blitz::array<T,N>> {
		static inline h5::blitz::array<T,N> ctor( std::array<size_t,N> dims ){
			::blitz::TinyVector<int,N> shape;
			for(int i=0; i<N; ++i) shape[i] = static_cast<int>(dims[i]);
			return h5::blitz::array<T,N>(shape);
	}};
}
#endif
#endif
