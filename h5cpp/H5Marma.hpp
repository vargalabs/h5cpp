/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

//#include "H5Tmeta.hpp"

#if defined(ARMA_INCLUDES) || defined(H5CPP_USE_ARMADILLO)
namespace h5::arma {
		template<class T> using rowvec = ::arma::Row<T>;
		template<class T> using colvec = ::arma::Col<T>;
		template<class T> using colmat = ::arma::Mat<T>;
		template<class T> using cube   = ::arma::Cube<T>;

		// is_linalg_type := filter
		template <class Object, class T = typename impl::decay<Object>::type> using is_supported =
		std::bool_constant<std::is_same_v<Object,h5::arma::cube<T>> || std::is_same_v<Object,h5::arma::colmat<T>>
			|| std::is_same_v<Object,h5::arma::rowvec<T>> ||  std::is_same_v<Object,h5::arma::colvec<T>>>;
}

namespace h5::meta {
    template <class T> struct is_contiguous<h5::arma::rowvec<T>> : std::true_type {};
    template <class T> struct is_contiguous<h5::arma::colvec<T>> : std::true_type {};
    template <class T> struct is_contiguous<h5::arma::colmat<T>> : std::true_type {};
    template <class T> struct is_contiguous<h5::arma::cube<T>> : std::true_type {};

    // Register types so generic access_traits_t fallbacks don't create ambiguous partial specializations
    template <class T> struct detail::has_explicit_access_traits<h5::arma::rowvec<T>> : std::true_type {};
    template <class T> struct detail::has_explicit_access_traits<h5::arma::colvec<T>> : std::true_type {};
    template <class T> struct detail::has_explicit_access_traits<h5::arma::colmat<T>> : std::true_type {};
    template <class T> struct detail::has_explicit_access_traits<h5::arma::cube<T>> : std::true_type {};
}

namespace h5::impl {
	// 1.) object -> H5T_xxx

	// Register Armadillo types so the structural decay fallback in H5Mstl.hpp
	// does not create an ambiguous partial-specialisation with these explicit specs.
	template <class T> struct detail::has_explicit_decay<h5::arma::rowvec<T>> : std::true_type {};
	template <class T> struct detail::has_explicit_decay<h5::arma::colvec<T>> : std::true_type {};
	template <class T> struct detail::has_explicit_decay<h5::arma::colmat<T>> : std::true_type {};
	template <class T> struct detail::has_explicit_decay<h5::arma::cube<T>>   : std::true_type {};

	template <class T> struct decay<h5::arma::rowvec<T>>{ using type = T; };
	template <class T> struct decay<h5::arma::colvec<T>>{ using type = T; };
	template <class T> struct decay<h5::arma::colmat<T>>{ using type = T; };
	template <class T> struct decay<h5::arma::cube<T>>{ using type = T; };

	// get read access to datastaore
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::arma::is_supported<Object>::value,
	const T*> data( const Object& ref ){
			return ref.memptr();
	}

	// read write access
	template <class Object, class T = typename impl::decay<Object>::type> inline
	std::enable_if_t< h5::arma::is_supported<Object>::value,
	T*> data( Object& ref ){
			return ref.memptr();
	}

	// rank
	template<class T> struct rank<h5::arma::rowvec<T>> : public std::integral_constant<size_t,1>{};
	template<class T> struct rank<h5::arma::colvec<T>> : public std::integral_constant<size_t,1>{};
	template<class T> struct rank<h5::arma::colmat<T>> : public std::integral_constant<size_t,2>{};
	template<class T> struct rank<h5::arma::cube<T>> : public std::integral_constant<size_t,3>{};

	// determine rank and dimensions
	template <class T> inline std::array<size_t,1> size( const h5::arma::rowvec<T>& ref ){ return {ref.n_elem};}
	template <class T> inline std::array<size_t,1> size( const h5::arma::colvec<T>& ref ){ return {ref.n_elem};}
	template <class T> inline std::array<size_t,2> size( const h5::arma::colmat<T>& ref ){ return {ref.n_rows,ref.n_cols};}
	template <class T> inline std::array<size_t,3> size( const h5::arma::cube<T>& ref ){ return {ref.n_slices,ref.n_cols,ref.n_rows};}

	// CTOR-s 
	template <class T> struct get<h5::arma::rowvec<T>> {
		static inline  h5::arma::rowvec<T> ctor( std::array<size_t,1> dims ){
			return h5::arma::rowvec<T>( dims[0] );
	}};
	template <class T> struct get<h5::arma::colvec<T>> {
		static inline h5::arma::colvec<T> ctor( std::array<size_t,1> dims ){
			return h5::arma::colvec<T>( dims[0] );
	}};
	template <class T> struct get<h5::arma::colmat<T>> {
		static inline h5::arma::colmat<T> ctor( std::array<size_t,2> dims ){
			return h5::arma::colmat<T>( dims[0], dims[1] );
	}};
	template <class T> struct get<h5::arma::cube<T>> {
		static inline h5::arma::colmat<T> ctor( std::array<size_t,3> dims ){
			return h5::arma::colmat<T>( dims[2], dims[0], dims[1] );
	}};
}

namespace h5::meta {
    template <class T> struct access_traits_t<h5::arma::rowvec<T>> {
        using element_t = T;
        static constexpr access_t kind = access_t::contiguous;
        static constexpr bool is_trivially_packable = true;
        static auto data(const h5::arma::rowvec<T>& c) noexcept { return h5::impl::data(c); }
        static auto size(const h5::arma::rowvec<T>& c) noexcept { return h5::impl::size(c); }
        static std::size_t bytes(const h5::arma::rowvec<T>& c) noexcept {
            auto s = size(c); std::size_t n = 1;
            for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
            return n * sizeof(element_t);
        }
    };
    template <class T> struct access_traits_t<h5::arma::colvec<T>> {
        using element_t = T;
        static constexpr access_t kind = access_t::contiguous;
        static constexpr bool is_trivially_packable = true;
        static auto data(const h5::arma::colvec<T>& c) noexcept { return h5::impl::data(c); }
        static auto size(const h5::arma::colvec<T>& c) noexcept { return h5::impl::size(c); }
        static std::size_t bytes(const h5::arma::colvec<T>& c) noexcept {
            auto s = size(c); std::size_t n = 1;
            for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
            return n * sizeof(element_t);
        }
    };
    template <class T> struct access_traits_t<h5::arma::colmat<T>> {
        using element_t = T;
        static constexpr access_t kind = access_t::contiguous;
        static constexpr bool is_trivially_packable = true;
        static auto data(const h5::arma::colmat<T>& c) noexcept { return h5::impl::data(c); }
        static auto size(const h5::arma::colmat<T>& c) noexcept { return h5::impl::size(c); }
        static std::size_t bytes(const h5::arma::colmat<T>& c) noexcept {
            auto s = size(c); std::size_t n = 1;
            for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
            return n * sizeof(element_t);
        }
    };
    template <class T> struct access_traits_t<h5::arma::cube<T>> {
        using element_t = T;
        static constexpr access_t kind = access_t::contiguous;
        static constexpr bool is_trivially_packable = true;
        static auto data(const h5::arma::cube<T>& c) noexcept { return h5::impl::data(c); }
        static auto size(const h5::arma::cube<T>& c) noexcept { return h5::impl::size(c); }
        static std::size_t bytes(const h5::arma::cube<T>& c) noexcept {
            auto s = size(c); std::size_t n = 1;
            for (std::size_t i = 0; i < s.size(); ++i) n *= s[i];
            return n * sizeof(element_t);
        }
    };
}
#endif
