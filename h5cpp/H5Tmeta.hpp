#pragma once

#include <hdf5.h>
#include <array>
#include <cstring>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace h5 {
    template<class T> struct is_registered : std::false_type {}; // is_supported<> is reserved for c++26 static reflection
    template<class T> hid_t inline register_struct();            // so is make_type<>()
}
// TODO: create detector for STL container types: direct access, iterable and container: direct access + iterable
namespace h5::meta {
    /** canonical type helpers */
    template<class T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
    /** detection idiom */
    template<class...> using void_t = void;
    template<class default_t, class always_void, template<class...> class op, class... args_t>
    struct detector {
        using value_t = std::false_type;
        using type = default_t;
    };
    template<class default_t, template<class...> class op, class... args_t>
    struct detector<default_t, void_t<op<args_t...>>, op, args_t...> {
        using value_t = std::true_type;
        using type = op<args_t...>;
    };
    template<template<class...> class op, class... args_t> using is_detected = typename detector<void, void, op, args_t...>::value_t;
    template<template<class...> class op, class... args_t> inline constexpr bool is_detected_v = is_detected<op, args_t...>::value;
    template<class default_t, template<class...> class op, class... args_t> using detected_or_t = typename detector<default_t, void_t<op<args_t...>>, op, args_t...>::type;

    /** type registry */
    template<class T> using register_struct_expr_t = decltype(h5::register_struct<T>());
    template<class T> struct has_register_struct : is_detected<register_struct_expr_t, remove_cvref_t<T>> {};
    template<class T> inline constexpr bool has_register_struct_v = has_register_struct<T>::value;    

    /** primitive expression detectors */
    template<class T> using value_type_t = typename T::value_type;
    template<class T> using key_type_t = typename T::key_type;
    template<class T> using mapped_type_t = typename T::mapped_type;
    template<class T> using key_compare_t = typename T::key_compare;
    template<class T> using hasher_t = typename T::hasher;

    template<class T> using size_expr_t = decltype(std::declval<const T&>().size());
    template<class T> using data_expr_t = decltype(std::declval<const T&>().data());
    template<class T> using begin_expr_t = decltype(std::declval<const T&>().begin());
    template<class T> using end_expr_t = decltype(std::declval<const T&>().end());
    /** semantic grouping: enumerated */
    template<class T> struct is_enumerated_like : std::is_enum<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_enumerated_like_v = is_enumerated_like<T>::value;
    /** semantic grouping: bitfield */
    template<class T> struct is_bitfield : std::false_type {};
    template<class alloc_t> struct is_bitfield<std::vector<bool, alloc_t>> : std::true_type {}; /*!< packed boolean proxy container treated as bitfield storage */
    template<class T> struct is_bitfield_like : is_bitfield<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_bitfield_like_v = is_bitfield_like<T>::value;
    /** semantic grouping: opaque */
    template<class T> struct is_opaque : std::false_type {};
    template<> struct is_opaque<void*> : std::true_type {};
    template<> struct is_opaque<const void*> : std::true_type {};
    template<> struct is_opaque<void**> : std::true_type {};
    template<> struct is_opaque<const void**> : std::true_type {};
    template<class T> struct is_opaque_like : is_opaque<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_opaque_like_v = is_opaque_like<T>::value;
    /** semantic grouping: compound */
    template<class T> struct is_compound_like : h5::is_registered<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_compound_like_v = is_compound_like<T>::value;
    /** semantic grouping: fixed-size text */
    template<class T> struct fixed_text_like : std::false_type {};
    template<std::size_t N> struct fixed_text_like<char[N]> : std::true_type {};
    template<std::size_t N> struct fixed_text_like<const char[N]> : std::true_type {};
    /** semantic grouping: variable-length text */
    template<class T> struct vl_text_like : std::false_type {};
    template<> struct vl_text_like<char*> : std::true_type {};                   // raw c-string pointers
    template<> struct vl_text_like<const char*> : std::true_type {};
    template<> struct vl_text_like<char* const> : std::true_type {};
    template<> struct vl_text_like<const char* const> : std::true_type {};
    template<> struct vl_text_like<char* volatile> : std::true_type {};
    template<> struct vl_text_like<const char* volatile> : std::true_type {};
    template<> struct vl_text_like<char* const volatile> : std::true_type {};
    template<> struct vl_text_like<const char* const volatile> : std::true_type {};
    template<class traits_t> struct vl_text_like<std::basic_string_view<char, traits_t>> : std::true_type {};               // std::string_view family
    template<class traits_t, class alloc_t> struct vl_text_like<std::basic_string<char, traits_t, alloc_t>> : std::true_type {}; // std::basic_string family
    /** semantic grouping: text */
    template<class T> struct text_like : std::bool_constant<fixed_text_like<T>::value || vl_text_like<T>::value> {};
    template<class T> struct is_text_like : text_like<std::remove_reference_t<T>> {};
    template<class T> inline constexpr bool is_text_like_v = is_text_like<T>::value;
    /** semantic grouping: fixed extent array-like */
    template<class T> struct is_array : std::false_type {};
    template<class T, std::size_t N> struct is_array<T[N]> : std::true_type {};            // built-in arrays
    template<class T, std::size_t N> struct is_array<std::array<T, N>> : std::true_type {};// std::array family
    template<class T> struct array_like : is_array<std::remove_reference_t<T>> {};         //  references are semantically irrelevant
    template<class T> inline constexpr bool is_array_v = array_like<T>::value;
    /** semantic grouping: is_iterable, has_direct_access or is_blas_like*/
    template<class T> struct is_iterable : std::bool_constant<is_detected_v<begin_expr_t, remove_cvref_t<T>> && is_detected_v<end_expr_t, remove_cvref_t<T>>>{};
    template<class T> inline constexpr bool is_iterable_v = is_iterable<T>::value;
    template<class T> struct has_direct_access : std::bool_constant< std::is_pointer_v<detected_or_t<void, data_expr_t, remove_cvref_t<T>> >> {};          
    template<class T> inline constexpr bool has_direct_access_v = has_direct_access<T>::value;
    template<class T> struct is_blas_like : std::false_type{};
    template<class T> inline constexpr bool is_blas_like_v = is_blas_like<T>::value;

    /** semantic grouping: map-like / set-like */
    template<class T> struct is_map_like : std::bool_constant<is_detected_v<key_type_t, remove_cvref_t<T>> && is_detected_v<mapped_type_t, remove_cvref_t<T>> && is_detected_v<value_type_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_map_like_v = is_map_like<T>::value;
    template<class T> struct is_set_like : std::bool_constant<is_detected_v<key_type_t, remove_cvref_t<T>> && is_detected_v<value_type_t, remove_cvref_t<T>> && !is_detected_v<mapped_type_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_set_like_v = is_set_like<T>::value;
    /** semantic grouping: sequential containers */
    template<class T> struct is_sequential_like : std::bool_constant<is_iterable_v<T> && is_detected_v<value_type_t, remove_cvref_t<T>> && !is_detected_v<key_type_t, remove_cvref_t<T>> && !is_detected_v<mapped_type_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_sequential_like_v = is_sequential_like<T>::value;
    /** semantic grouping: ordered associative containers */
    template<class T> struct is_associative_like : std::bool_constant<is_iterable_v<T> && is_detected_v<key_type_t, remove_cvref_t<T>> && is_detected_v<key_compare_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_associative_like_v = is_associative_like<T>::value;
    /** semantic grouping: unordered associative containers */
    template<class T> struct is_unordered_like : std::bool_constant<is_iterable_v<T> && is_detected_v<key_type_t, remove_cvref_t<T>> && is_detected_v<hasher_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_unordered_like_v = is_unordered_like<T>::value;
    /** semantic grouping: stl-like containers */
    template<class T> struct is_stl_like : std::bool_constant<is_sequential_like_v<T> || is_associative_like_v<T> || is_unordered_like_v<T>> {};
    template<class T> inline constexpr bool is_stl_like_v = is_stl_like<T>::value;

    /** recursively compute semantic leaf element type */
    template<class T, class = void> struct decay { using type = remove_cvref_t<T>; };
    template<class T> struct decay<T, std::enable_if_t<fixed_text_like<remove_cvref_t<T>>::value>> { using type = char; };
    template<class T> struct decay<T, std::enable_if_t<vl_text_like<remove_cvref_t<T>>::value>> {  using type = char; };
    template<class T, std::size_t N> struct decay<T[N], void> { using type = typename decay<T>::type; };
    template<class T, std::size_t N> struct decay<const T[N], void> { using type = typename decay<T>::type; };
    template<class T, std::size_t N> struct decay<std::array<T,N>, void> { using type = typename decay<T>::type; };
    
    template<class T, class alloc_t> struct decay<std::vector<T,alloc_t>, void> { using type = typename decay<T>::type; };
    template<class T, class traits_t> struct decay<std::basic_string_view<T,traits_t>, std::enable_if_t<std::is_same_v<T,char>>> { using type = char; };
    template<class T, class traits_t, class alloc_t> struct decay<std::basic_string<T,traits_t,alloc_t>, std::enable_if_t<std::is_same_v<T,char>>> { using type = char; };
    template<class alloc_t> struct decay<std::vector<bool,alloc_t>, void> { using type = bool; };
    template<class T> using decay_t = typename decay<T>::type;