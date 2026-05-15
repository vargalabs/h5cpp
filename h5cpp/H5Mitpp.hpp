/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <hdf5.h>

#if defined(MAT_H) || defined(H5CPP_USE_ITPP_MATRIX)
namespace h5::itpp {
		template<class T> using rowmat = ::itpp::Mat<T>;
		template <class Object, class T = typename impl::decay<Object>::type>
			using is_supported = std::bool_constant<std::is_same_v<Object,h5::itpp::rowmat<T>>>;
}
namespace h5::meta {
		template <class T> struct is_contiguous<h5::itpp::rowmat<T>> : std::true_type {};

		// Register types so generic access_traits_t fallbacks don't create ambiguous partial specializations
		template <class T> struct detail::has_explicit_access_traits<h5::itpp::rowmat<T>> : std::true_type {};
}
namespace h5::impl {
	// 1.) object -> H5T_xxx
	template <class T> struct detail::has_explicit_decay<h5::itpp::rowmat<T>> : std::true_type {};
	template <class T> struct decay<h5::itpp::rowmat<T>>{ using type = T; };

	// get read access to datastaore
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::itpp::is_supported<Object>::value,
	const T*> data(const Object& ref ){
			return ref._data();
	}
	// read write access
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::itpp::is_supported<Object>::value,
	T*> data( Object& ref ){
			return ref._data();
	}
	template<class T> struct rank<h5::itpp::rowmat<T>> : public std::integral_constant<size_t,2>{};
	template <class T> inline std::array<size_t,2> size( const h5::itpp::rowmat<T>& ref ){ return {(hsize_t)ref.cols(),(hsize_t)ref.rows()};}
	template <class T> struct get<h5::itpp::rowmat<T>> {
		static inline h5::itpp::rowmat<T> ctor( std::array<size_t,2> dims ){
			return h5::itpp::rowmat<T>( dims[1], dims[0] );
	}};
}
namespace h5::meta {
		template <class T> struct access_traits_t<h5::itpp::rowmat<T>> {
			using element_t = T;
			static constexpr access_t kind = access_t::contiguous;
			static constexpr bool is_trivially_packable = true;
			static auto data(const h5::itpp::rowmat<T>& c) noexcept { return h5::impl::data(c); }
			static auto size(const h5::itpp::rowmat<T>& c) noexcept { return h5::impl::size(c); }
			static std::size_t bytes(const h5::itpp::rowmat<T>& c) noexcept {
				auto s = size(c); std::size_t n = 1;
				for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
				return n * sizeof(element_t);
			}
		};
}
#endif

#if defined(VEC_H) || defined(H5CPP_USE_ITPP_VECTOR)
namespace h5::itpp {
		template<class T> using rowvec = ::itpp::Vec<T>;
		template <class Object, class T = typename impl::decay<Object>::type>
			using is_supported_v = std::bool_constant<std::is_same_v<Object,h5::itpp::rowvec<T>>>;
}
namespace h5::meta {
		template <class T> struct is_contiguous<h5::itpp::rowvec<T>> : std::true_type {};

		// Register types so generic access_traits_t fallbacks don't create ambiguous partial specializations
		template <class T> struct detail::has_explicit_access_traits<h5::itpp::rowvec<T>> : std::true_type {};
}
namespace h5::impl {
	template <class T> struct detail::has_explicit_decay<h5::itpp::rowvec<T>> : std::true_type {};
	template <class T> struct decay<h5::itpp::rowvec<T>>{ using type = T; };
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::itpp::is_supported_v<Object>::value,
	const T*> data(const Object& ref ){
			return ref._data();
	}
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::itpp::is_supported_v<Object>::value,
	T*> data( Object& ref ){
			return ref._data();
	}
	template<class T> struct rank<h5::itpp::rowvec<T>> : public std::integral_constant<size_t,1>{};
	template <class T> inline std::array<size_t,1> size( const h5::itpp::rowvec<T>& ref ){ return { (hsize_t)ref.size() };}
	template <class T> struct get<h5::itpp::rowvec<T>> {
		static inline h5::itpp::rowvec<T> ctor( std::array<size_t,1> dims ){
			return h5::itpp::rowvec<T>( dims[0] );
	}};
}
namespace h5::meta {
		template <class T> struct access_traits_t<h5::itpp::rowvec<T>> {
			using element_t = T;
			static constexpr access_t kind = access_t::contiguous;
			static constexpr bool is_trivially_packable = true;
			static auto data(const h5::itpp::rowvec<T>& c) noexcept { return h5::impl::data(c); }
			static auto size(const h5::itpp::rowvec<T>& c) noexcept { return h5::impl::size(c); }
			static std::size_t bytes(const h5::itpp::rowvec<T>& c) noexcept {
				auto s = size(c); std::size_t n = 1;
				for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
				return n * sizeof(element_t);
			}
		};
}
#endif
