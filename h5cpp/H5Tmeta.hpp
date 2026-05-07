#pragma once

#include <hdf5.h>
#include <array>
#include <cstring>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>
#include <stdexcept>

namespace h5 {
    template<class T> struct is_registered : std::false_type {}; // is_supported<> is reserved for c++26 static reflection
    template<class T> hid_t inline register_struct();            // so is make_type<>()
}
// TODO: create detector for STL container types: direct access, iterable and container: direct access + iterable
namespace h5::meta {

    /** canonical type helpers */
    template<class T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
    template<class> inline constexpr bool unreachable_v = false;
    template<class> inline constexpr bool dependent_failure_v = false;
    template<class> inline constexpr bool substitution_failure_v = false;
    template<class> inline constexpr bool instantiation_error_v = false;
    template<class> inline constexpr bool not_supported_v = false;
    namespace blas {
        template<class T> size_t n_elements(const remove_cvref_t<T>& value);
        template<class T> size_t extent(const remove_cvref_t<T>& value);
        template<class T> size_t data(remove_cvref_t<T>& value);
        template<class T> size_t data(const remove_cvref_t<T>& value);
    }    
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
    template<class T> using resize_expr_t = decltype(std::declval<T&>().resize(std::declval<std::size_t>()));
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
    template<class T> struct is_fixed_text_like : fixed_text_like<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_fixed_text_like_v = is_fixed_text_like<T>::value;    
    /** semantic grouping: sentinel text are C-like text */
    template<class T> struct sentinel_text_like : std::false_type {};
    template<class T> struct sentinel_text_like<T*> : std::bool_constant<std::is_same_v<std::remove_cv_t<T>, char>> {};
    template<class T> struct is_sentinel_text_like : sentinel_text_like<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_sentinel_text_like_v = is_sentinel_text_like<T>::value;
    /** semantic grouping: dynamic text is std::string variant*/
    template<class T> struct dynamic_text_like : std::false_type {};
    template<class traits_t> struct dynamic_text_like<std::basic_string_view<char, traits_t>> : std::true_type {};
    template<class traits_t, class alloc_t> struct dynamic_text_like<std::basic_string<char, traits_t, alloc_t>> : std::true_type {};
    template<class T> struct is_dynamic_text_like : dynamic_text_like<remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_dynamic_text_like_v = is_dynamic_text_like<T>::value;
    /** semantic grouping: variable length string */
    template<class T> struct vl_text_like : std::bool_constant<is_sentinel_text_like_v<T> || is_dynamic_text_like_v<T>> {};
    template<class T> struct is_vl_text_like : vl_text_like<T> {};
    template<class T> inline constexpr bool is_vl_text_like_v = is_vl_text_like<T>::value;    
    /** semantic grouping: text */
    template<class T> struct text_like : std::bool_constant<is_fixed_text_like_v<T> || is_sentinel_text_like_v<T> || is_dynamic_text_like_v<T>>{};
    template<class T> struct is_text_like : text_like<T> {};
    template<class T> inline constexpr bool is_text_like_v = is_text_like<T>::value;

    /** semantic grouping: fixed extent array-like */
    template<class T> struct is_builtin_array : std::false_type {};
    template<class T, std::size_t N> struct is_builtin_array<T[N]> : std::true_type {};
    template<class T> inline constexpr bool is_builtin_array_v = is_builtin_array<remove_cvref_t<T>>::value;
    template<class T> struct is_array : std::false_type {};
    template<class T, std::size_t N> struct is_array<T[N]> : std::true_type {};            // built-in arrays
    template<class T, std::size_t N> struct is_array<std::array<T, N>> : std::true_type {};// std::array family
    template<class T> struct array_like : is_array<remove_cvref_t<T>> {};         //  references are semantically irrelevant
    template<class T> inline constexpr bool is_array_v = array_like<T>::value;
    /** semantic grouping: is_iterable, has_direct_access, has_resize or is_blas_like*/
    template<class T> struct is_iterable : std::bool_constant<is_detected_v<begin_expr_t, remove_cvref_t<T>> && is_detected_v<end_expr_t, remove_cvref_t<T>>>{};
    template<class T> inline constexpr bool is_iterable_v = is_iterable<T>::value;
    template<class T> struct has_direct_access : std::bool_constant< std::is_pointer_v<detected_or_t<void, data_expr_t, remove_cvref_t<T>> >> {};          
    template<class T> inline constexpr bool has_direct_access_v = has_direct_access<T>::value;
    template<class T> struct is_resizable : is_detected<resize_expr_t, remove_cvref_t<T>> {};
    template<class T> inline constexpr bool is_resizable_v = is_resizable<T>::value;
    template<class T> struct has_data : is_detected<data_expr_t, remove_cvref_t<T>> {};
    template<class T> inline constexpr bool has_data_v = has_data<T>::value;
    template<class T> struct has_size : is_detected<size_expr_t, remove_cvref_t<T>> {};
    template<class T> inline constexpr bool has_size_v = has_size<T>::value;
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
    template<class T> inline constexpr bool is_stl_like_v = is_stl_like<remove_cvref_t<T>>::value;

    /** recursively compute semantic leaf element type; FIXME: std::complex<T> must be a leaf type */
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

    /** bulk-transfer leaf: scalar, enum, or registered compound */
    template<class T> struct is_contiguous_leaf : std::bool_constant<std::is_arithmetic_v<remove_cvref_t<T>> || std::is_enum_v<remove_cvref_t<T>> || h5::is_registered<remove_cvref_t<T>>::value> {};
    template<class T> inline constexpr bool is_contiguous_leaf_v = is_contiguous_leaf<T>::value;
    /** elements that break inline recursive embedding */
    template<class T> struct breaks_inline_embedding : std::false_type {};
    template<class T> struct breaks_inline_embedding<T*> : std::true_type {};
    template<class char_t, class traits_t> struct breaks_inline_embedding<std::basic_string_view<char_t, traits_t>> : std::true_type {};
    template<class char_t, class traits_t, class alloc_t> struct breaks_inline_embedding<std::basic_string<char_t, traits_t, alloc_t>> : std::true_type {};
    template<class alloc_t> struct breaks_inline_embedding<std::vector<bool, alloc_t>> : std::true_type {};
    template<class T> inline constexpr bool breaks_inline_embedding_v = breaks_inline_embedding<remove_cvref_t<T>>::value;
    /** recursive bulk-transfer legality */
    template<class T, class = void> struct has_legal_recursive_path : std::false_type {};
    template<class T> struct has_legal_recursive_path<T, std::enable_if_t<is_contiguous_leaf_v<T>>> : std::true_type {}; // leaf
    /** built-in and std::array preserve inline embedding */
    template<class T, std::size_t N> struct has_legal_recursive_path<T[N], void> : has_legal_recursive_path<T> {};
    template<class T, std::size_t N> struct has_legal_recursive_path<std::array<T, N>, void> : has_legal_recursive_path<T> {};
    /** direct-access containers are legal only at level 1 and only iff their element recursively do so */
    template<class T> struct has_legal_recursive_path<T, std::enable_if_t<has_direct_access_v<T> && is_detected_v<value_type_t, remove_cvref_t<T>> && !breaks_inline_embedding_v<T>>> : has_legal_recursive_path<typename remove_cvref_t<T>::value_type> {};
    template<class T> inline constexpr bool has_legal_recursive_path_v = has_legal_recursive_path<T>::value;
    /** top-level special case: text and BLAS objects are transferable as a single blob */
    template<class T> struct is_top_level_contiguous : std::bool_constant<is_text_like_v<T> || is_blas_like_v<T>> {};
    template<class T> inline constexpr bool is_top_level_contiguous_v = is_top_level_contiguous<T>::value;
    template<class T> struct is_recursively_contiguous : std::bool_constant<is_top_level_contiguous_v<T> || has_legal_recursive_path_v<T>> {};
    template<class T> inline constexpr bool is_recursively_contiguous_v = is_recursively_contiguous<T>::value;

    /** rank - recursively nested containers */
    template<class T, class = void> struct rank : std::integral_constant<std::size_t, 0> {};
    template<class T> struct rank<T, std::enable_if_t<is_text_like_v<T>>> : std::integral_constant<std::size_t, 0> {};
    template<class T> struct rank<T, std::enable_if_t<is_map_like_v<T> || is_set_like_v<T> || is_associative_like_v<T> || is_unordered_like_v<T>>> : std::integral_constant<std::size_t, 0> {};
    template<class T> struct rank<T, std::enable_if_t<std::is_array_v<T>>> : std::integral_constant<std::size_t, std::rank_v<remove_cvref_t<T>>> {};
    template<class T> struct rank<T, std::enable_if_t<is_sequential_like_v<T> >> : std::integral_constant<std::size_t,1 + rank<typename remove_cvref_t<T>::value_type>::value> {};
    template<class T>  inline constexpr std::size_t rank_v = rank<T>::value;

    /** extent - recursively nested */
    template<std::size_t... dims> struct extents_t {};
    template<class T, std::size_t n, std::size_t... dims> struct extents_i : extents_i<T, n - 1, std::extent_v<T, n - 1>, dims...> {};
    template<class T, std::size_t... dims> struct extents_i<T, 0, dims...> { using type = extents_t<dims...>; };
    template<class T> using extents_of_t = typename extents_i<remove_cvref_t<T>, std::rank_v<remove_cvref_t<T>>>::type;
    template<std::size_t... dims> constexpr auto to_array(extents_t<dims...>) { return std::array<hsize_t, sizeof...(dims)>{ static_cast<hsize_t>(dims)... }; }
    template<class T> constexpr auto extents() { return to_array(extents_of_t<T>{});}

    /** trait enums */
    enum class storage_t { scalar, enumerated, bitfield, compound, array, fixed_text, vl_text, opaque, unsupported }; /*!< what am I storing? */
    enum class shape_t { scalar, fixed_extents, runtime_extents, explicit_only };                                     /*!< how do I obtain extents? */
    enum class access_t { direct, indirect, custom, unsupported };                                                    /*!< how do I transfer data? */
    
    enum class linalg_layout_t {none, contiguous, strided, packed };                                                  /*!< how is payload laid out? */
    enum class linalg_order_t { none, row_major, column_major };
    enum class matrix_structure_t { none, general, triangular, diagonal, symmetric, hermitian, banded };
    enum class matrix_triangle_t { none, upper, lower };
    enum class matrix_diagonal_t { none, explicit_nonunit, unit };
    enum class matrix_transpose_t { none, transposed, conjugate_transposed };

    /** storage traits */
    template<class T, class = void> struct storage_traits_t {
        using user_t = T; using value_t = remove_cvref_t<T>; using element_t = decay_t<T>;
        static constexpr storage_t kind = is_opaque_like_v<value_t> ? storage_t::opaque : is_bitfield_like_v<value_t> ? storage_t::bitfield :
            is_enumerated_like_v<value_t> ? storage_t::enumerated : is_compound_like_v<value_t> ? storage_t::compound :
            fixed_text_like<value_t>::value? storage_t::fixed_text :  vl_text_like<value_t>::value   ? storage_t::vl_text :
            is_array_v<value_t> ? storage_t::array : std::is_scalar<value_t>::value ? storage_t::scalar : storage_t::unsupported;
    };
    template<class T> inline constexpr storage_t storage_v = storage_traits_t<T>::kind;
    /** shape traits */
    template<class T, class = void> struct shape_traits_t {
        using user_t = T; using value_t = remove_cvref_t<T>; using element_t = decay_t<T>;
        static constexpr shape_t kind = std::is_scalar<value_t>::value || is_enumerated_like_v<value_t> ? shape_t::scalar :
            fixed_text_like<value_t>::value ? shape_t::fixed_extents : is_array_v<value_t> ? shape_t::fixed_extents :
            vl_text_like<value_t>::value ? shape_t::runtime_extents : is_blas_like_v<value_t> ? shape_t::runtime_extents :
            (is_iterable_v<value_t> && is_detected_v<size_expr_t, value_t>) ? shape_t::runtime_extents : shape_t::explicit_only;
    };
    template<class T> inline constexpr shape_t shape_v = shape_traits_t<T>::kind;

    /** access traits */
    template<class T, class = void>
    struct access_traits_t {
        using user_t = T; using value_t = remove_cvref_t<T>; using element_t = decay_t<value_t>; using size_type = std::size_t;
        static constexpr access_t kind = access_t::unsupported;
        static constexpr std::size_t rank = 0; /*!< rank of the access-visible payload */
        static size_type n_elements(const value_t&) = delete; /*!< logical element count */
        static size_type size(const value_t&) = delete; /*!< n_elements * sizeof(element_t) */
        static std::array<std::size_t, rank> extent(const value_t&) = delete; /*!< returns rank sized extents */
        static const element_t* data(const value_t&) = delete; /*!< borrowed payload pointer */
        static element_t* data(value_t&) = delete; /*!< borrowed mutable payload pointer */
        static void scatter(const value_t&, element_t*) = delete; /*!< packs exactly size(x) leaf elements into canonical flat order */
        static void gather(value_t&, const element_t*) = delete;  /*!< reconstructs object from exactly size(x) leaf elements in canonical flat order */    
    };
    
    template<class T> struct access_traits_t<T, std::enable_if_t<is_recursively_contiguous_v<T> >> {
        using user_t = T; using value_t = remove_cvref_t<T>; using element_t = decay_t<value_t>; using size_type = std::size_t;
        static constexpr access_t kind = access_t::direct;
        static constexpr std::size_t rank = rank_v<value_t>;
        static size_type n_elements(const value_t& value) {
            if constexpr (is_builtin_array_v<value_t>) {
                constexpr auto extents = to_array(extents_of_t<value_t>{});
                std::size_t product = 1; for (std::size_t i = 0; i < extents.size(); ++i) product *= extents[i];
                return product;
            } else if constexpr (is_blas_like_v<value_t> ) 
                return blas::n_elements<value_t>(value);
            else if constexpr (rank == 0) // note: string is 1 element/line
                return 1; 
            else if constexpr (has_size_v<value_t>)
                return value.size();
            else static_assert(not_supported_v<T>, "Could not infer size for type");
        }
        static size_type size(const value_t& value) {
            if constexpr (is_dynamic_text_like_v<value_t> && has_size_v<value_t>)
                return value.size() * sizeof(element_t);
            if constexpr (is_sentinel_text_like_v<value_t>)
                return std::strlen(value);
            if constexpr (is_fixed_text_like_v<value_t> || is_builtin_array_v<value_t> )
                return sizeof(value_t);
            else return n_elements(value) * sizeof(element_t);
        }
        static std::array<std::size_t, rank> extent(const value_t& value){
            if constexpr (is_blas_like_v<value_t> )
                return blas::extent<value_t>(value);
            else if constexpr (is_builtin_array_v<value_t>)
                return to_array(extents_of_t<value_t>{});
            else if constexpr (rank == 0) 
                return {};
            else if constexpr (rank == 1 && has_size_v<value_t>)
                return { value.size() };
            else static_assert(not_supported_v<T>, "access_traits_t<T>::extent(): could not infer extents for type");
        }
        template<class U> static auto data_impl(U& value) ->decltype(auto) {
            using data_ptr_t = std::conditional_t<std::is_const_v<std::remove_reference_t<U>>, const element_t*, element_t*>;
            if constexpr (is_blas_like_v<value_t>) {
                return blas::data<value_t>(value);
            } else if constexpr (rank == 0) {
                if constexpr (is_dynamic_text_like_v<value_t>) return value.empty() ? nullptr : value.data();
                if constexpr (is_fixed_text_like_v<value_t> || is_sentinel_text_like_v<value_t>) return value;
                return reinterpret_cast<data_ptr_t>(&value);
            } else if constexpr (is_builtin_array_v<value_t>) // T[]..
                return reinterpret_cast<data_ptr_t>(&value);
            else if constexpr (has_data_v<value_t>) // vector like contoners
                return value.data();
            else  static_assert(not_supported_v<T>, "access_traits_t<T>::data(): could not obtain mutable data pointer");
        }

        static element_t* data(value_t& value) { return data_impl(value); }
        static const element_t* data(const value_t& value) { return data_impl(value); }
        static void scatter(const value_t&, element_t*) = delete;
        static void gather(value_t&, const element_t*) = delete;        
    };
    
    template<class T> struct access_traits_t<T, std::enable_if_t<!is_recursively_contiguous_v<T> >> {
        using user_t = T; using value_t = remove_cvref_t<T>; using element_t = decay_t<value_t>; using size_type = std::size_t;
        static constexpr access_t kind = access_t::indirect;
        static constexpr std::size_t rank = rank_v<value_t>;
    };



} // namespace h5::meta

/*
direct, indirect, custom, unsupported
reference, address, accessor, indirect, custom, reference
contiguous, indirect, custom, unsupported

template struct<class T> type {
    using storage = storage_trait_t<T>;
    using shape = shape_traits<T>;
    using access = access_traits<T>;

};

void read([...]){
    if constexpr (access::kind == access_traits_t::direct){
        using meta::type
    } else if constexpr (access::kind == access_traits_t::indirect){
    } [..]
}


*/