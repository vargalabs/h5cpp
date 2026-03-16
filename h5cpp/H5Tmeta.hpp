
/* SPDX-License-Identifier: MIT
 * This file is part of H5CPP.
 * Copyright (c) 2018-2024 Varga Consulting, Toronto, ON, Canada.
 * Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. */

#pragma once

#include <type_traits>
#include <string>
#include <string_view>
#include <array>
#include <vector>
#include <cstddef>
#include <tuple>
#include <complex>
#include "H5meta.hpp"
#include "H5Iall.hpp"

//FIXME: move it elsewhere
#define H5CPP_supported_elementary_types "supported elementary types ::= pod_struct | float | double |  [signed](int8 | int16 | int32 | int64)"

// stl detection with templates, this probably should stay until concepts become mainstream
namespace h5::meta {
    template<class T>  struct is_array : std::is_array<T> {};
    template<class T, std::size_t N>  struct is_array<std::array<T, N>> : std::true_type {};
    template<class T> using is_stl = /* will use concepts once become available */
        std::disjunction<meta::has_iterator<T>, meta::has_data<T>,meta::has_size<T>,meta::has_value_type<T>>;
    template<class... Ts> using enable_or = std::enable_if_t<std::disjunction_v<Ts...>, int>;
    template<class... Ts> using enable_and = std::enable_if_t<std::conjunction_v<Ts...>, int>;
    template<class T> struct decay { using type = typename meta::value<T>::type; };
    template<class T, std::size_t N> struct decay<T[N]> { using type = typename std::remove_all_extents<T[N]>::type; };

    /* std::array<...,R> size<T>(){} template has to compute rank R at compile 
     * time, these templates, and their respective specializations aid to accomplish that*/
    template<class T> struct rank : public std::integral_constant<std::size_t, meta::has_size<T>::value> {}; // definition
    template<class U> struct rank<U[]>  : public std::integral_constant<std::size_t, rank<U>::value + 1U>{};
    template<class U, std::size_t N> struct rank<U[N]> : public std::integral_constant<std::size_t, rank<U>::value + 1U>{};
    template<class U, std::size_t N> struct rank<U*[N]> : public std::integral_constant<std::size_t, rank<U*>::value>{};
    template<std::size_t N> struct rank<char[N]> : std::integral_constant<std::size_t, 0U> {}; // character literals
    template<class T, std::size_t N> using is_rank = std::integral_constant<bool, rank<T>::value == N>;
    /** helpers for is_rank<T>; define behavior via rank<T> specializations */
    template<class T> using is_scalar = is_rank<T, 0U>;
    template<class T> using is_vector = is_rank<T, 1U>;
    template<class T> using is_matrix = is_rank<T, 2U>;
    template<class T> using is_cube = is_rank<T, 3U>;

    namespace impl {
        template<class T> struct is_string : std::false_type {};
        template<class char_t, class traits_t, class alloc_t> struct is_string<std::basic_string<char_t, traits_t, alloc_t>> : std::true_type {};
        template<class char_t, class traits_t> struct is_string<std::basic_string_view<char_t, traits_t>> : std::true_type {};
    }
    template<class T, class D = typename meta::decay<T>::type> using is_string = std::disjunction<impl::is_string<T>, impl::is_string<D>>;

    /* Objects may reside in continuous memory region such as vectors, matrices, POD structures can be saved/loaded in a single transfer,
     * the rest needs to be handled on a member variable bases*/
    template<class T> struct is_contiguous : std::integral_constant<bool, compat::is_pod<T>::value> {};
    template<class T, class... Ts> struct is_contiguous<std::vector<T, Ts...>> : std::integral_constant<bool, compat::is_pod<T>::value> {};
    template<class T, class traits_t, class alloc_t> struct is_contiguous<std::basic_string<T, traits_t, alloc_t>> : std::true_type {};
    template<class T, class traits_t> struct is_contiguous<std::basic_string_view<T, traits_t>> : std::true_type {};
    template <std::size_t N> struct is_contiguous <const char*[N]> : std::false_type {};

    template <class T> struct is_contiguous <std::complex<T>> : std::true_type{};
    template <class T, std::size_t N> struct is_contiguous <std::array<T,N>> : std::integral_constant<bool, compat::is_pod<T>::value>{};

    template <class T, class... Ts> struct is_linalg : std::false_type {};
    template <class C, class T, class... Cs> struct is_valid : std::false_type {};

    // DEFAULT CASE
    template <class T> struct rank<T*>: public std::integral_constant<std::size_t,1U>{}; // pointers are treated as rank-1 views; pointer indirection does not accumulate tensor rank
    template<class T> std::enable_if_t<!std::is_array<T>::value, const T*> data(const T& ref) { return &ref; }
    template<class T> std::enable_if_t<meta::has_size<T>::value, std::array<std::size_t, 1U>> size(const T& ref) { return {ref.size()}; } 
    template <class T, class... Ts> struct get {
        static inline T ctor( std::array<std::size_t,0> dims ){
            return T(); }};
    // ARRAYS
    namespace detail { // recursive replacement for the former raw-array staircase overloads
        template<class T> constexpr T* array_data(T& ref) noexcept { return &ref; }
        template<class T> constexpr const T* array_data(const T& ref) noexcept { return &ref; }
        template<class T, std::size_t N> constexpr std::remove_all_extents_t<T>* array_data(T (&ref)[N]) noexcept { return array_data(ref[0]); }
        template<class T, std::size_t N> constexpr const std::remove_all_extents_t<T>* array_data(const T (&ref)[N]) noexcept { return array_data(ref[0]); }
    } // namespace detail

    template<class T, std::size_t N>  constexpr std::remove_all_extents_t<T>* data(T (&ref)[N]) noexcept { return detail::array_data(ref); }
    template<class T, std::size_t N> constexpr const std::remove_all_extents_t<T>* data(const T (&ref)[N]) noexcept { return detail::array_data(ref); }
    template <class T, std::size_t N> std::array<std::size_t, 1U> size(const T(&)[N]) {
        constexpr auto extents = h5::meta::get_extent<T[N]>();
        std::size_t product = 1;
        for (std::size_t i = 0; i < extents.size(); ++i) product *= extents[i];
        return {product};
    }
    //STD::STRING
    template<class T, class traits_t, class alloc_t> struct rank<std::basic_string<T, traits_t, alloc_t>> : std::integral_constant<std::size_t, 0U> {};
    template<class T, class traits_t> struct rank<std::basic_string_view<T, traits_t>> : std::integral_constant<std::size_t, 0U> {};    
    template <class T, class... Ts> std::array<std::size_t,1U> size( const std::basic_string<T, Ts...>& ref ){ return{ref.size()}; }
    template<class T, class... Ts> struct get<std::basic_string<T, Ts...>> {
        static inline std::basic_string<T, Ts...> ctor(std::array<std::size_t, 1U> dims) {
            return std::basic_string<T, Ts...>(dims[0], T{});
        }
    };

    // STD::INITIALIZER_LIST<T>
    template<std::size_t N0> struct rank<std::initializer_list<char[N0]>> : std::integral_constant<std::size_t, 1U> {};
    template<std::size_t N0, std::size_t N1> struct rank<std::initializer_list<char[N0][N1]>> : std::integral_constant<std::size_t, 2U> {};    
    template<class T> struct rank<std::initializer_list<T>>: public std::integral_constant<std::size_t,1U>{};

    // STD::VECTOR<T>
    template<class T, class... Ts> struct rank<std::vector<T, Ts...>> : std::integral_constant<std::size_t, 1U> {};
    template <class T, class... Ts> std::array<std::size_t,1U> size(const std::vector<T, Ts...>& ref) {
        if constexpr (h5::meta::rank<T>::value != 0) {
            std::size_t total = 0;
            for (const auto& x : ref)
                total += h5::meta::size(x)[0];
            return {total};
        } else return {ref.size()};
    }
    // STD::ARRAY<T>
    // 3.) read access
    template <class T, std::size_t N> inline const T* data( const std::array<T,N>& ref ){ return ref.data(); }
    template <class T, std::size_t N> inline T* data( std::array<T,N>& ref ){ return ref.data(); }
    template <class T, std::size_t N> inline std::array<std::size_t,1U> size(const std::array<T,N>& ref) {
        if constexpr (h5::meta::rank<T>::value != 0) {
            std::size_t total = 0;
            for (const auto& x : ref)
                total += h5::meta::size(x)[0];
            return {total};
        } else return {N};
    }


    template <class T, std::size_t N> struct rank<std::array<T,N>> : public std::integral_constant<std::size_t,1U> {};
    // END STD::ARRAY

    template <class T> void get_fields( T& sp ){}
    template <class T> void get_field_names( T& sp ){}
    template <class T> void get_field_attributes( T& sp ){}

    // NON_CONTIGUOUS 
    template <class T> struct member {
        using type = std::tuple<void>;
        static constexpr std::size_t size = 0U;
    };
    template <class T> using csc_t = std::tuple< // compressed sparse column: indices, indptr, values
            std::vector<unsigned long>, std::vector<unsigned long>, std::vector<T>>;
    constexpr std::tuple<const char*, const char*, const char*> csc_names = {"indices", "indptr","data"};
// HID_T
//    template <class T> struct decay<h5::dt_t<T>>{ using type = T; };
}
namespace h5::meta::linalg {
    /*types accepted by BLAS/LAPACK*/
    using blas = std::tuple<float,double,std::complex<float>,std::complex<double>>;
}
namespace h5::meta {
    // what handles may have attributes:  h5::acreate, h5::aread, h5::awrite
    template <class T, class... Ts> struct has_attribute : std::false_type {};
    template <> struct has_attribute<h5::gr_t> : std::true_type {};
    template <> struct has_attribute<h5::ds_t> : std::true_type {};
    template <> struct has_attribute<h5::ob_t> : std::true_type {};
//    template <class T> struct has_attribute<h5::dt_t<T>> : std::true_type {};
    template <class T, class... Ts> struct is_location : std::false_type {};
    template <> struct is_location<h5::gr_t> : std::true_type {};
    template <> struct is_location<h5::fd_t> : std::true_type {};
    //template <class T> struct 
}
