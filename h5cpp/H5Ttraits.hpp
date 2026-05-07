#pragma once

#include <hdf5.h>
#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace h5::meta {
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

    /** primitive expression detectors */
    template<class T> using value_type_t = typename T::value_type;
    template<class T> using key_type_t = typename T::key_type;
    template<class T> using mapped_type_t = typename T::mapped_type;

    template<class T> using size_expr_t = decltype(std::declval<const T&>().size());
    template<class T> using data_expr_t = decltype(std::declval<const T&>().data());
    template<class T> using begin_expr_t = decltype(std::declval<const T&>().begin());
    template<class T> using end_expr_t = decltype(std::declval<const T&>().end());

    template<class T> using rows_expr_t = decltype(std::declval<const T&>().rows());
    template<class T> using cols_expr_t = decltype(std::declval<const T&>().cols());

    template<class T> using memptr_expr_t = decltype(std::declval<const T&>().memptr());

    template<class T> using n_rows_expr_t = decltype(std::declval<const T&>().n_rows);
    template<class T> using n_cols_expr_t = decltype(std::declval<const T&>().n_cols);
    template<class T> using n_slices_expr_t = decltype(std::declval<const T&>().n_slices);

    /** canonical type helpers */
    template<class T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

    /** semantic grouping */
    template<class T> struct is_c_string_like : std::false_type {};
    template<> struct is_c_string_like<char*> : std::true_type {};
    template<> struct is_c_string_like<const char*> : std::true_type {};
    template<class T> struct is_std_string_like : std::false_type {};
    template<> struct is_std_string_like<std::string> : std::true_type {};
    template<class T> struct is_string_like : std::disjunction<is_c_string_like<remove_cvref_t<T>>,is_std_string_like<remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_string_like_v = is_string_like<T>::value;
    template<class T> struct is_string_view_like : std::false_type {};
    template<class char_t, class traits_t> struct is_string_view_like<std::basic_string_view<char_t, traits_t>> : std::true_type {};
    template<class T> inline constexpr bool is_string_view_like_v = is_string_view_like<remove_cvref_t<T>>::value;
    template<class T> struct is_std_vector_bool : std::false_type {};
    template<class alloc_t> struct is_std_vector_bool<std::vector<bool, alloc_t>> : std::true_type {};
    template<class T> inline constexpr bool is_std_vector_bool_v = is_std_vector_bool<remove_cvref_t<T>>::value;
    template<class T> struct is_char_array : std::false_type {};
    template<std::size_t N> struct is_char_array<char[N]> : std::true_type {};
    template<std::size_t N> struct is_char_array<const char[N]> : std::true_type {};
    template<class T> inline constexpr bool is_char_array_v = is_char_array<remove_cvref_t<T>>::value;
    template<class T> struct is_fixed_extent_array_like : std::false_type {};
    template<class T, std::size_t N> struct is_fixed_extent_array_like<T[N]> : std::true_type {};
    template<class T, std::size_t N> struct is_fixed_extent_array_like<std::array<T, N>> : std::true_type {};
    template<class T> inline constexpr bool is_fixed_extent_array_like_v = is_fixed_extent_array_like<remove_cvref_t<T>>::value;
    template<class T> struct is_sequence_like : std::bool_constant<
        is_detected_v<value_type_t, remove_cvref_t<T>> && is_detected_v<size_expr_t, remove_cvref_t<T>> && // 
        is_detected_v<begin_expr_t, remove_cvref_t<T>> && is_detected_v<end_expr_t, remove_cvref_t<T>> &&  // sequence must have iterators
        !is_string_like_v<T> && !is_string_view_like_v<T> // but strings are special class
    > {};
    template<class T> inline constexpr bool is_sequence_like_v = is_sequence_like<T>::value;
    template<class T>
    struct is_contiguous_sequence_like : std::bool_constant<
        is_sequence_like_v<T> && ( is_detected_v<data_expr_t, remove_cvref_t<T>> || is_detected_v<memptr_expr_t, remove_cvref_t<T>> )> {};
    template<class T> inline constexpr bool is_contiguous_sequence_like_v = is_contiguous_sequence_like<T>::value;
    template<class T> struct is_matrix_like : std::bool_constant<is_detected_v<value_type_t, remove_cvref_t<T>> &&
        ((is_detected_v<rows_expr_t, remove_cvref_t<T>> && is_detected_v<cols_expr_t, remove_cvref_t<T>>) ||
            (is_detected_v<n_rows_expr_t, remove_cvref_t<T>> && is_detected_v<n_cols_expr_t, remove_cvref_t<T>>)) &&
        ( is_detected_v<data_expr_t, remove_cvref_t<T>> || is_detected_v<memptr_expr_t, remove_cvref_t<T>>) && !is_string_like_v<T>> {};
    template<class T> inline constexpr bool is_matrix_like_v = is_matrix_like<T>::value;
    template<class T> struct is_cube_like : std::bool_constant<
        is_detected_v<value_type_t, remove_cvref_t<T>> && is_detected_v<n_rows_expr_t, remove_cvref_t<T>> && is_detected_v<n_cols_expr_t, remove_cvref_t<T>> &&
        is_detected_v<n_slices_expr_t, remove_cvref_t<T>> && ( is_detected_v<data_expr_t, remove_cvref_t<T>> || is_detected_v<memptr_expr_t, remove_cvref_t<T>> ) && !is_string_like_v<T>> {};
    template<class T> inline constexpr bool is_cube_like_v = is_cube_like<T>::value;
    template<class T> struct is_map_like : std::bool_constant<
        is_detected_v<key_type_t, remove_cvref_t<T>> && is_detected_v<mapped_type_t, remove_cvref_t<T>> && is_detected_v<value_type_t, remove_cvref_t<T>>> {};
    template<class T> inline constexpr bool is_map_like_v = is_map_like<T>::value;
    /** recursive contiguity true iff T presents a single recursively flat direct memory image.
     *  any pointer, string, string-view, or proxy-like leaf breaks contiguity. */
    template<class T, class = void> struct is_contiguous : std::false_type {};
    template<class T> inline constexpr bool is_contiguous_v = is_contiguous<T>::value;
    template<class T> struct is_contiguous<T, std::enable_if_t<
        (std::is_arithmetic<remove_cvref_t<T>>::value || std::is_enum<remove_cvref_t<T>>::value) &&
        !std::is_pointer<remove_cvref_t<T>>::value && !is_string_like_v<T> && !is_string_view_like_v<T>>> : std::true_type {};
    template<class T> struct is_contiguous<T, std::enable_if_t<std::is_pointer<remove_cvref_t<T>>::value>> : std::false_type {};
    template<class T> struct is_contiguous<T, std::enable_if_t<is_string_like_v<T> || is_string_view_like_v<T>>> : std::false_type {};
    template<class T, std::size_t N> struct is_contiguous<T[N], void> : is_contiguous<T> {};
    template<class T, std::size_t N> struct is_contiguous<std::array<T, N>, void> : is_contiguous<T> {};
    template<class T, class alloc_t> struct is_contiguous<std::vector<T, alloc_t>, void> : std::bool_constant<!is_std_vector_bool_v<std::vector<T, alloc_t>> && is_contiguous<T>::value> {};
    template<class T> struct is_contiguous<T, std::enable_if_t<
        is_contiguous_sequence_like_v<T> && !is_matrix_like_v<T> && !is_cube_like_v<T> &&
        !is_map_like_v<T> && !is_std_vector_bool_v<T> >> : std::bool_constant<is_contiguous<typename remove_cvref_t<T>::value_type>::value> {};
    template<class T> struct is_contiguous<T, std::enable_if_t<is_matrix_like_v<T>>> 
        : std::bool_constant<is_contiguous<typename remove_cvref_t<T>::value_type>::value> {};
    template<class T> struct is_contiguous<T, std::enable_if_t<is_cube_like_v<T>>> 
        : std::bool_constant<is_contiguous<typename remove_cvref_t<T>::value_type>::value> {};
    template<class T> struct is_contiguous<T, std::enable_if_t<is_map_like_v<T>>> : std::false_type {};


    /** trait enums */
    enum class storage_kind_t { 
        scalar, compound, fixed_array, fixed_string, vl_string, sequence_of_vl_string, opaque, unsupported };
    enum class shape_kind_t { scalar, fixed_extents, runtime_extents, explicit_only };
    enum class access_kind_t { object_address, contiguous_payload, pointer_table, gather, custom, unsupported };
    enum class layout_kind_t { none, contiguous, strided, packed, triangular, diagonal, symmetric, hermitian, banded, transposed };
    enum class matrix_order_t { none, row_major, column_major };
    enum class triangle_kind_t { none, upper, lower };
    enum class diagonal_kind_t { implicit, explicit_nonunit, unit };
    

    /** storage traits */
    template<class T, class = void> struct storage_traits_t {
        using user_t = T;
        using value_t = remove_cvref_t<T>;
        using element_t = value_t;
        static constexpr storage_kind_t kind =
            std::is_scalar<value_t>::value ? storage_kind_t::scalar : storage_kind_t::unsupported;
    };
    template<> struct storage_traits_t<char*> {
        using user_t = char*;
        using value_t = char*;
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::vl_string;
    };
    template<> struct storage_traits_t<const char*> {
        using user_t = const char*;
        using value_t = const char*;
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::vl_string;
    };
    template<> struct storage_traits_t<std::string> {
        using user_t = std::string;
        using value_t = std::string;
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::vl_string;
    };
    template<class char_t, class traits_t> struct storage_traits_t<std::basic_string_view<char_t, traits_t>> {
        using user_t = std::basic_string_view<char_t, traits_t>; // FIXME: have to think through all cases
        using value_t = user_t;
        using element_t = char_t;
        static constexpr storage_kind_t kind = storage_kind_t::unsupported;
    };
    template<std::size_t N> struct storage_traits_t<char[N]> {
        using user_t = char[N];
        using value_t = char[N];
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::fixed_string;
    };
    template<std::size_t N> struct storage_traits_t<const char[N]> {
        using user_t = const char[N];
        using value_t = const char[N];
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::fixed_string;
    };
    template<class T, std::size_t N> struct storage_traits_t<T[N], std::enable_if_t<!std::is_same<char, std::remove_cv_t<T>>::value>> {
        using user_t = T[N];
        using value_t = T[N];
        using element_t = typename storage_traits_t<T>::element_t;
        static constexpr storage_kind_t kind = storage_kind_t::fixed_array;
    };
    template<class T, std::size_t N> struct storage_traits_t<std::array<T, N>> {
        using user_t  = std::array<T, N>;
        using value_t = std::array<T, N>;
        using element_t = typename storage_traits_t<T>::element_t;
        static constexpr storage_kind_t kind = storage_kind_t::fixed_array;
    };
    template<class T, class alloc_t> struct storage_traits_t<std::vector<T, alloc_t>> {
        using user_t = std::vector<T, alloc_t>;
        using value_t = std::vector<T, alloc_t>;
        using element_t = typename storage_traits_t<T>::element_t;
        static constexpr storage_kind_t kind = storage_traits_t<T>::kind; // HDF5 storage class, not the outer C++ container category.
    };
    template<class alloc_t> struct storage_traits_t<std::vector<std::string, alloc_t>> {
        using user_t = std::vector<std::string, alloc_t>;
        using value_t = std::vector<std::string, alloc_t>;
        using element_t = char;
        static constexpr storage_kind_t kind = storage_kind_t::sequence_of_vl_string;
    };

    /** shape traits */
    template<class T, class = void> struct shape_traits_t {
        using value_t = remove_cvref_t<T>;
        static constexpr shape_kind_t kind = shape_kind_t::scalar;
        static constexpr std::size_t rank = 0;
    };
    template<class T, std::size_t N> struct shape_traits_t<T[N]> {
        using value_t = T[N];
        static constexpr shape_kind_t kind = shape_kind_t::fixed_extents;
        static constexpr std::size_t rank = 1 + shape_traits_t<T>::rank;
    };
    template<class T, std::size_t N> struct shape_traits_t<std::array<T, N>> {
        using value_t = std::array<T, N>;
        static constexpr shape_kind_t kind = shape_kind_t::fixed_extents;
        static constexpr std::size_t rank = 1 + shape_traits_t<T>::rank;
    };
    template<class T, class alloc_t> struct shape_traits_t<std::vector<T, alloc_t>> {
        using value_t = std::vector<T, alloc_t>;
        static constexpr shape_kind_t kind = shape_kind_t::runtime_extents;
        static constexpr std::size_t rank = 1 + shape_traits_t<T>::rank; /* recursively defined*/
    };

    /** access traits */
    template<class T, class = void> struct access_traits_t {
        using value_t = remove_cvref_t<T>;
        using pointer_t = const value_t*;
        static constexpr access_kind_t kind = access_kind_t::object_address;
        static pointer_t data(const value_t& ref) { return &ref; }
    };
    template<class T, std::size_t N> struct access_traits_t<T[N]> {
        using value_t = T[N];
        using pointer_t = const T*;
        static constexpr access_kind_t kind = access_kind_t::contiguous_payload;
        static pointer_t data(const value_t& ref) { return &ref[0]; }
    };
    template<class T, std::size_t N> struct access_traits_t<std::array<T, N>> {
        using value_t = std::array<T, N>;
        using pointer_t = const T*;
        static constexpr access_kind_t kind = access_kind_t::contiguous_payload;
        static pointer_t data(const value_t& ref) { return ref.data(); }
    };
  
    template<class alloc_t> struct access_traits_t<std::vector<std::string, alloc_t>> {
        using value_t   = std::vector<std::string, alloc_t>;
        using pointer_t = const char* const*;
        static constexpr access_kind_t kind = access_kind_t::pointer_table;
        /** Holds transient c_str() pointers into the source strings; this object must outlive any HDF5 call using data(...). */
        struct marshalled_t { std::vector<const char*> table; };
        static marshalled_t marshal(const value_t& ref) {
            marshalled_t out;
            out.table.reserve(ref.size());
            for (const auto& s : ref)
                out.table.push_back(s.c_str());
            return out;
        }
        static pointer_t data(const marshalled_t& ref) {
            return ref.table.data();
        }
    };
    template<class T> struct access_traits_t<T, std::enable_if_t<is_matrix_like_v<T>>> {
        using value_t = remove_cvref_t<T>;
        using element_t = typename value_t::value_type;
        using pointer_t = const element_t*;
        static constexpr access_kind_t kind = is_contiguous_v<element_t> ? access_kind_t::contiguous_payload : access_kind_t::gather;

        static pointer_t data(const value_t& ref) {
            if constexpr (is_detected_v<data_expr_t, value_t>)
                return ref.data();
            else return ref.memptr();
        }
    };

    template<class T>
    struct access_traits_t<T, std::enable_if_t<is_cube_like_v<T>>> {
        using value_t = remove_cvref_t<T>;
        using element_t = typename value_t::value_type;
        using pointer_t = const element_t*;
        static constexpr access_kind_t kind = is_contiguous_v<element_t> ? access_kind_t::contiguous_payload : access_kind_t::gather;

        static pointer_t data(const value_t& ref) {
            if constexpr (is_detected_v<data_expr_t, value_t>)
                return ref.data();
            else return ref.memptr();
        }
    };
    template<class T> struct access_traits_t<T, std::enable_if_t< is_sequence_like_v<T> && is_detected_v<data_expr_t, remove_cvref_t<T>> &&
        !is_matrix_like_v<T> && !is_cube_like_v<T> && !is_map_like_v<T> && !is_std_vector_bool_v<T> >> {
        using value_t   = remove_cvref_t<T>;
        using element_t = typename value_t::value_type;
        using pointer_t = const element_t*;

        static constexpr access_kind_t kind =
            is_contiguous_v<element_t> ? access_kind_t::contiguous_payload : access_kind_t::gather;
        static pointer_t data(const value_t& ref) { return ref.data(); }
    };
    template<class T> struct access_traits_t<T, std::enable_if_t<is_pod_like_v<T> && has_type_descriptor_t<T>::value>> {
        using value_t   = remove_cvref_t<T>;
        using pointer_t = const value_t*;
        static constexpr access_kind_t kind = access_kind_t::object_address;
        static pointer_t data(const value_t& ref) { return &ref; }
    };    
    /** layout traits (for BLAS like objects) */ 
    template<class T, class = void> struct layout_traits_t {
        static constexpr matrix_order_t order = matrix_order_t::none;
        static constexpr bool is_transposed = false;
        static constexpr triangle_kind_t triangle = triangle_kind_t::none;
        static constexpr bool is_strided = false;
        static constexpr bool is_packed = false;
    };

    /** dataset defaults */
    template<class T> struct dataset_defaults_t {
        static constexpr bool recommend_chunking = false;
        static constexpr bool supports_unlimited = false;
        static constexpr bool requires_vl_heap = false;
    };

    template<class T, class alloc_t> struct dataset_defaults_t<std::vector<T, alloc_t>> {
        static constexpr bool recommend_chunking = false;
        static constexpr bool supports_unlimited = true;
        static constexpr bool requires_vl_heap = false;
    };
    template<class alloc_t> struct dataset_defaults_t<std::vector<std::string, alloc_t>> {
        static constexpr bool recommend_chunking = true;
        static constexpr bool supports_unlimited = true;
        static constexpr bool requires_vl_heap = true;
    };

    /** public convenience aliases */
    template<class T> using storage_element_t = typename storage_traits_t<T>::element_t;
    template<class T> inline constexpr storage_kind_t storage_kind_v = storage_traits_t<T>::kind;
    template<class T> inline constexpr shape_kind_t shape_kind_v = shape_traits_t<T>::kind;
    template<class T> inline constexpr access_kind_t access_kind_v = access_traits_t<T>::kind;
} // namespace h5::meta