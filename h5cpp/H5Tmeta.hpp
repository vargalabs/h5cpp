
/*
 * Copyright (c) 2018 - 2021 vargaconsulting, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#ifndef  H5CPP_META456_HPP 
#define  H5CPP_META456_HPP

#include "H5Iall.hpp"
#include "H5meta.hpp"
#include <type_traits>
#include <string>
#include <string_view>
#include <array>
#include <vector>
#include <deque>
#include <forward_list>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>
#include <any>
#include <tuple>
#include <complex>

//FIXME: move it elsewhere
#define H5CPP_supported_elementary_types "supported elementary types ::= pod_struct | float | double |  [signed](int8 | int16 | int32 | int64)"

// stl detection with templates, this probably should stay until concepts become mainstream
namespace h5::meta {
 
    template <class T, class... Ts> struct is_array : public std::is_array<T>{};
    template <class T, size_t N> struct is_array <std::array<T,N>> : std::true_type{};

    template<class T> using is_stl = /* will use concepts once become available */
        std::disjunction<meta::has_iterator<T>, meta::has_data<T>,meta::has_size<T>,meta::has_value_type<T>>;

    template <class... Args> using enable_or = std::enable_if<std::disjunction_v<Args...>>;
    template <class... Args> using enable_and = std::enable_if<std::conjunction_v<Args...>>;
    template <class T, class... Ts> struct decay {
        using type = typename meta::value<T>::type; };
    template <class T, size_t N> struct decay<T[N]>{ // support for array types
        using type = typename std::remove_all_extents<T>::type; };

    /* std::array<...,R> size<T>(){} template has to compute rank R at compile 
     * time, these templates, and their respective specializations aid to accomplish that*/
    template<class T, class... Ts> struct rank : public std::integral_constant<int,
        meta::has_size<T>::value> {}; // definition
    template<class U> struct rank<U[]>  : public std::integral_constant<std::size_t, rank<U>::value + 1>{};
    template<class U, std::size_t N> struct rank<U[N]> : public std::integral_constant<std::size_t, rank<U>::value + 1>{};
    template<class U, std::size_t N> struct rank<U*[N]> : public std::integral_constant<std::size_t, rank<U*>::value>{};
    template <size_t N> struct rank<char[N]> : public std::integral_constant<int,0>{}; // character literals

    template<class T, int N, class... Ts> using is_rank = std::integral_constant<bool, rank<T, Ts...>::value == N >;
    // helpers for is_rank<T>, don't need specialization, instead define 'rank'
    template<class T, class... Ts> using is_scalar = is_rank<T,0,Ts...>; // numerical | pod 
    template<class T, class... Ts> using is_vector = is_rank<T,1,Ts...>;
    template<class T, class... Ts> using is_matrix = is_rank<T,2,Ts...>;
    template<class T, class... Ts> using is_cube   = is_rank<T,3,Ts...>;

    template <class T, class D=typename meta::decay<T>::type>
    using is_string = typename std::integral_constant<bool,
        std::is_same<T,std::basic_string<char>>::value || std::is_same<D, std::basic_string<char>>::value || 
        std::is_same<T,std::basic_string<wchar_t>>::value || std::is_same<D, std::basic_string<wchar_t>>::value || 
        std::is_same<T,std::basic_string<char16_t>>::value || std::is_same<D, std::basic_string<char16_t>>::value || 
        std::is_same<T,std::basic_string<char32_t>>::value || std::is_same<D, std::basic_string<char32_t>>::value ||
        std::is_same<T,std::basic_string_view<char>>::value || std::is_same<D, std::basic_string<char>>::value || 
        std::is_same<T,std::basic_string_view<wchar_t>>::value || std::is_same<D, std::basic_string_view<wchar_t>>::value || 
        std::is_same<T,std::basic_string_view<char16_t>>::value || std::is_same<D, std::basic_string_view<char16_t>>::value || 
        std::is_same<T,std::basic_string_view<char32_t>>::value || std::is_same<D, std::basic_string_view<char32_t>>::value>;


    /* Objects may reside in continuous memory region such as vectors, matrices, POD structures can be saved/loaded in a single transfer,
     * the rest needs to be handled on a member variable bases*/
    template <class T, class... Ts> struct is_contiguous : std::integral_constant<bool, std::is_pod<T>::value> {};
    template <class T, class... Ts> struct is_contiguous <std::basic_string<T,Ts...>> : std::true_type {};
    template <class T, class... Ts> struct is_contiguous <std::basic_string_view<T,Ts...>> : std::true_type {};
    template <size_t N> struct is_contiguous <const char*[N]> : std::false_type {};

    template <class T> struct is_contiguous <std::complex<T>> : std::true_type{};
    template <class T, class... Ts> struct is_contiguous <std::vector<T,Ts...>> :
        std::integral_constant<bool, std::is_pod<T>::value>{};
    template <class T, size_t N> struct is_contiguous <std::array<T,N>> :
        std::integral_constant<bool, std::is_pod<T>::value>{};

    template <class T, class... Ts> struct is_linalg : std::false_type {};
    template <class C, class T, class... Cs> struct is_valid : std::false_type {};

    // ------------------------------------------------------------
    // Capability-based type families (#86)
    // Additive predicate vocabulary used by the upcoming #87/#88/#89/#90
    // refactors. Predicates only; no storage/shape/access traits, no recursive
    // contiguity, no datatype synthesis, no dataset I/O changes. C++17 floor.
    // Reuses existing detection primitives from H5meta.hpp; does not redeclare
    // has_data, has_size, has_direct_access, has_value_type, has_iterator,
    // or has_const_iterator.
    // ------------------------------------------------------------

    template <class T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

    // Internal expression-detector aliases (consistent with H5meta.hpp `_f` style)
    template <class T> using key_type_f    = typename T::key_type;
    template <class T> using mapped_type_f = typename T::mapped_type;
    template <class T> using key_compare_f = typename T::key_compare;
    template <class T> using hasher_f      = typename T::hasher;
    template <class T> using resize_f      = decltype(std::declval<T&>().resize(std::declval<std::size_t>()));

    // text-like
    template <class T> struct is_fixed_text_like : std::false_type {};
    template <std::size_t N> struct is_fixed_text_like<char[N]>       : std::true_type {};
    template <std::size_t N> struct is_fixed_text_like<const char[N]> : std::true_type {};

    template <class T> struct is_vl_text_like : std::false_type {};
    template <> struct is_vl_text_like<char*>       : std::true_type {};
    template <> struct is_vl_text_like<const char*> : std::true_type {};
    template <class Tr> struct is_vl_text_like<std::basic_string_view<char, Tr>> : std::true_type {};
    template <class Tr, class A> struct is_vl_text_like<std::basic_string<char, Tr, A>> : std::true_type {};

    template <class T> struct is_text_like
        : std::bool_constant<is_fixed_text_like<remove_cvref_t<T>>::value
                          || is_vl_text_like<remove_cvref_t<T>>::value> {};

    // array-like (built-in arrays + std::array, cv/ref-stable)
    namespace detail_capabilities {
        template <class T> struct is_array_like_impl : std::false_type {};
        template <class T, std::size_t N> struct is_array_like_impl<T[N]>            : std::true_type {};
        template <class T, std::size_t N> struct is_array_like_impl<std::array<T,N>> : std::true_type {};
    }
    template <class T> struct is_array_like
        : detail_capabilities::is_array_like_impl<remove_cvref_t<T>> {};

    // iterable — derived from existing has_iterator (begin && end)
    template <class T> struct is_iterable : has_iterator<remove_cvref_t<T>> {};

    // resizable
    template <class T> struct is_resizable : compat::is_detected<resize_f, remove_cvref_t<T>> {};

    // sequential / associative / unordered / set / map / stl_like
    // Text-like types (std::string, string_view) are excluded so they remain
    // classified under is_text_like rather than leaking into is_stl_like.
    template <class T> struct is_sequential_like
        : std::bool_constant<is_iterable<T>::value
                          && compat::is_detected<value_type_f,  remove_cvref_t<T>>::value
                          && !compat::is_detected<key_type_f,    remove_cvref_t<T>>::value
                          && !compat::is_detected<mapped_type_f, remove_cvref_t<T>>::value
                          && !is_text_like<T>::value> {};

    template <class T> struct is_associative_like
        : std::bool_constant<is_iterable<T>::value
                          && compat::is_detected<key_type_f,    remove_cvref_t<T>>::value
                          && compat::is_detected<key_compare_f, remove_cvref_t<T>>::value> {};

    template <class T> struct is_unordered_like
        : std::bool_constant<is_iterable<T>::value
                          && compat::is_detected<key_type_f, remove_cvref_t<T>>::value
                          && compat::is_detected<hasher_f,   remove_cvref_t<T>>::value> {};

    template <class T> struct is_set_like
        : std::bool_constant<compat::is_detected<key_type_f,    remove_cvref_t<T>>::value
                          && compat::is_detected<value_type_f,  remove_cvref_t<T>>::value
                          && !compat::is_detected<mapped_type_f, remove_cvref_t<T>>::value> {};

    template <class T> struct is_map_like
        : std::bool_constant<compat::is_detected<key_type_f,    remove_cvref_t<T>>::value
                          && compat::is_detected<mapped_type_f, remove_cvref_t<T>>::value
                          && compat::is_detected<value_type_f,  remove_cvref_t<T>>::value> {};

    template <class T> struct is_stl_like
        : std::bool_constant<is_sequential_like<T>::value
                          || is_associative_like<T>::value
                          || is_unordered_like<T>::value> {};

    // enumerated / bitfield / opaque
    template <class T> struct is_enumerated_like : std::is_enum<remove_cvref_t<T>> {};

    template <class T> struct is_bitfield_like : std::false_type {};
    template <class A> struct is_bitfield_like<std::vector<bool, A>> : std::true_type {};

    template <class T> struct is_opaque_like : std::false_type {};
    template <> struct is_opaque_like<void*>        : std::true_type {};
    template <> struct is_opaque_like<const void*>  : std::true_type {};
    template <> struct is_opaque_like<void**>       : std::true_type {};
    template <> struct is_opaque_like<const void**> : std::true_type {};

    // strict pointer-data capability:
    //   true iff T::data() exists AND its return type is a pointer.
    // Distinct from has_direct_access (which only detects existence of .data()).
    template <class T> struct has_data_pointer
        : std::bool_constant<std::is_pointer_v<
              compat::detected_or_t<void, data_f, remove_cvref_t<T>>>> {};
    // ------------------------------------------------------------
    // End capability-based type families (#86)
    // ------------------------------------------------------------

    // DEFAULT CASE
    template <class T> struct rank<T*>: public std::integral_constant<size_t,1>{};
    template <class T, class... Ts>
        typename std::enable_if<!std::is_array<T>::value, const T*>::type data(const T& ref ){ return &ref; };
    template <class T, class... Ts> 
        typename std::enable_if<meta::has_size<T>::value, std::array<size_t,1>
        >::type size(const T& ref){
        return {ref.size()};
    };
    template <class T, size_t N>
        std::array<size_t,1> size(const T(&ref)[N]){ return {N};};
    template <class T, class... Ts> struct get {
        static inline T ctor( std::array<size_t,0> dims ){
            return T(); }};
    // ARRAYS

    //already defined line 58: template <class T, int N> struct rank<T[N]> : public std::rank<T[N]>{};

    template <class T,int N0> const T* data( const T(&ref)[N0]){ return &ref[0];};
    template <class T,int N1,int N0> const T* data( const T(&ref)[N1][N0]){ return &ref[0][0];};
    template <class T,int N2,int N1,int N0> const T* data( const T(&ref)[N2][N1][N0]){ return &ref[0][0][0];};
    template <class T,int N3,int N2,int N1,int N0> const T* data( const T(&ref)[N3][N2][N1][N0]){ return &ref[0][0][0][0];};
    template <class T,int N4,int N3,int N2,int N1,int N0> const T* data( const T(&ref)[N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0];};
    template <class T,int N5,int N4,int N3,int N2,int N1,int N0> const T* data( const T(&ref)[N5][N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0][0];};
    template <class T,int N6,int N5,int N4,int N3,int N2,int N1,int N0> const T* data( const T(&ref)[N6][N5][N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0][0][0];};

    template <class T,int N0>  T* data(T(&ref)[N0]){ return &ref[0];};
    template <class T,int N1,int N0>T* data(T(&ref)[N1][N0]){ return &ref[0][0];};
    template <class T,int N2,int N1,int N0>T* data(T(&ref)[N2][N1][N0]){ return &ref[0][0][0];};
    template <class T,int N3,int N2,int N1,int N0>T* data(T(&ref)[N3][N2][N1][N0]){ return &ref[0][0][0][0];};
    template <class T,int N4,int N3,int N2,int N1,int N0>T* data(T(&ref)[N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0];};
    template <class T,int N5,int N4,int N3,int N2,int N1,int N0>T* data(T(&ref)[N5][N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0][0];};
    template <class T,int N6,int N5,int N4,int N3,int N2,int N1,int N0>T* data(T(&ref)[N6][N5][N4][N3][N2][N1][N0]){ return &ref[0][0][0][0][0][0][0];};

    template <class T, int N> std::array<size_t, std::rank<T[N]>::value>
        size(const T* ref ){ return  h5::meta::get_extent<T[N]>(); };
   // template <class T, int N> struct get {
   //     static inline T[N] ctor( std::array<size_t,0> dims ){
   //         return T[N](); }};


    //STD::STRING
    template<> struct rank<std::basic_string<char>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string<wchar_t>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string<char16_t>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string<char32_t>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string_view<char>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string_view<wchar_t>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string_view<char16_t>>: public std::integral_constant<size_t,0>{};
    template<> struct rank<std::basic_string_view<char32_t>>: public std::integral_constant<size_t,0>{};
  
   // template <class T, class... Ts> T const* data(const std::basic_string<T, Ts...>& ref ){
   //     return ref.data();
   // }
    template <class T, class... Ts> std::array<size_t,1> size( const std::basic_string<T, Ts...>& ref ){ return{ref.size()}; }
    template <class T, class... Ts> struct get<std::basic_string<T,Ts...>> {
        static inline std::basic_string<T,Ts...> ctor( std::array<size_t,1> dims ){
            return std::basic_string<T,Ts...>(); }};

// STD::INITIALIZER_LIST<T>
    template<int N0> struct rank<std::initializer_list<char[N0]>>: public std::integral_constant<size_t,1>{};
    template<int N0, int N1> struct rank<std::initializer_list<char[N0][N1]>>: public std::integral_constant<size_t,2>{};
    template<class T> struct rank<std::initializer_list<T>>: public std::integral_constant<size_t,1>{};

// STD::VECTOR<T>
    template<class T> struct rank<std::vector<T>>: public std::integral_constant<size_t,1>{};
    
    template <class T, class... Ts>
    std::array<size_t,1> size( const std::vector<T, Ts...>& ref ){ return{ref.size()}; }

// STD::ARRAY<T>
    // 3.) read access
    template <class T, size_t N> inline const T* data( const std::array<T,N>& ref ){ return ref.data(); }
    template <class T, size_t N> inline T* data( std::array<T,N>& ref ){ return ref.data(); }

    template <class T, size_t N> inline typename std::array<size_t,1> size( const std::array<T,N>& ref ){ return {N}; }
    template <class T, size_t N> struct rank<std::array<T,N>> : public std::integral_constant<int,1> {};
// END STD::ARRAY

    template <class T> void get_fields( T& sp ){}
    template <class T> void get_field_names( T& sp ){}
    template <class T> void get_field_attributes( T& sp ){}

// NON_CONTIGUOUS 
    template <class T> struct member {
        using type = std::tuple<void>;
        static constexpr size_t size = 0;
    };
    template <class T> using csc_t = std::tuple< //compresses sparse row: index, colptr, values
            std::vector<unsigned long>, std::vector<unsigned long>, std::vector<T>>;
    const constexpr std::tuple<const char*, const char*, const char*> 
        csc_names = {"indices", "indptr","data"};
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

#endif
