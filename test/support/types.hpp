/*
 * Copyright (c) 2019-2026 vargaconsulting, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 *
 * Cartesian product I/O type machinery for h5cpp test suite.
 * Issue #114.
 */
#pragma once

#include <complex>
#include <array>
#include <vector>
#include <deque>
#include <forward_list>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <tuple>
#include <type_traits>
#include <h5cpp/H5meta.hpp>
#include <h5cpp/H5Tmeta.hpp>

namespace h5::test {

    // ----------------------------------------------------------------
    // pod_t — simple POD struct usable as HDF5 compound element
    // ----------------------------------------------------------------
    struct pod_t {
        int    a;
        float  b;
        double c;

        bool operator==(const pod_t& o) const { return a == o.a && b == o.b && c == o.c; }
        bool operator!=(const pod_t& o) const { return !(*this == o); }
        bool operator< (const pod_t& o) const { return a < o.a; }
    };

    // ----------------------------------------------------------------
    // Numeric element types
    // ----------------------------------------------------------------
    using numerical_t = std::tuple<
        unsigned char, unsigned short, unsigned int, unsigned long long int,
        char, short, int, long long int, float, double>;

    // element_t = numerical types + pod_t
    using element_t = typename h5::impl::cat<numerical_t, std::tuple<pod_t>>::type;

    // ----------------------------------------------------------------
    // C-array family: ranks 1–7
    // Rank-7 is the documented supported upper bound (issue #115).
    // N=2 keeps array sizes at 2/4/8/16/32/64/128 elements — manageable.
    // ----------------------------------------------------------------
    template <class T, size_t N = 2>
    using array_t = std::tuple<
        T[N],
        T[N][N],
        T[N][N][N],
        T[N][N][N][N],
        T[N][N][N][N][N],
        T[N][N][N][N][N][N],
        T[N][N][N][N][N][N][N]>;

    // ----------------------------------------------------------------
    // Sequential containers with fixed inner extent
    // ----------------------------------------------------------------
    template <class T, size_t N = 4>
    using sequential_t = std::tuple<
        std::array<T, N>,
        std::vector<T>,
        std::vector<std::array<T, N>>,
        std::vector<std::vector<T>>>;

    // ----------------------------------------------------------------
    // Linear / set-like containers
    // ----------------------------------------------------------------
    template <class T>
    using linear_t = std::tuple<
        std::list<T>,
        std::forward_list<T>,
        std::deque<T>,
        std::set<T>,
        std::multiset<T>,
        std::unordered_set<T>,
        std::unordered_multiset<T>>;

    // ----------------------------------------------------------------
    // Map-like containers
    // ----------------------------------------------------------------
    template <class K, class V = K>
    using map_t = std::tuple<
        std::map<K, V>,
        std::multimap<K, V>,
        std::unordered_map<K, V>,
        std::unordered_multimap<K, V>>;

    // ----------------------------------------------------------------
    // Single-param container aliases
    // ----------------------------------------------------------------
    template <class T> using vec_t     = std::vector<T>;
    template <class T> using list_t    = std::list<T>;
    template <class T> using deque_t   = std::deque<T>;
    template <class T> using fwdlist_t = std::forward_list<T>;

    // ----------------------------------------------------------------
    // TMP utilities
    // ----------------------------------------------------------------

    // filter_t<Tuple, Pred> — keeps only types for which Pred<T>::value is true
    namespace detail {
        template <class Tuple, template <class> class Pred, class Out>
        struct filter_impl;

        template <template <class> class Pred, class Out>
        struct filter_impl<std::tuple<>, Pred, Out> {
            using type = Out;
        };

        template <class Head, class... Tail, template <class> class Pred, class... Out>
        struct filter_impl<std::tuple<Head, Tail...>, Pred, std::tuple<Out...>> {
            using type = typename filter_impl<
                std::tuple<Tail...>,
                Pred,
                typename std::conditional<
                    Pred<Head>::value,
                    std::tuple<Out..., Head>,
                    std::tuple<Out...>>::type>::type;
        };
    }

    template <class Tuple, template <class> class Pred>
    using filter_t = typename detail::filter_impl<Tuple, Pred, std::tuple<>>::type;

    // bind_all_t<C, Tuple> — given a single-arg template C and a tuple of types,
    // produce a tuple<C<T1>, C<T2>, ...>
    namespace detail {
        template <template <class> class C, class Tuple>
        struct bind_all_impl;

        template <template <class> class C, class... Ts>
        struct bind_all_impl<C, std::tuple<Ts...>> {
            using type = std::tuple<C<Ts>...>;
        };
    }

    template <template <class> class C, class Tuple>
    using bind_all_t = typename detail::bind_all_impl<C, Tuple>::type;

    // ----------------------------------------------------------------
    // is_supported<T> predicate
    // ----------------------------------------------------------------
    template <class T>
    struct is_supported : std::integral_constant<bool,
        h5::meta::storage_representation_v<T> != h5::meta::storage_representation_t::unsupported> {};

    // ----------------------------------------------------------------
    // io_types — the final filtered type list used in TEST_CASE_TEMPLATE_APPLY
    //
    // Concatenates: numerical scalars + c-arrays of int + sequential
    // containers of int + linear containers of int + maps of int
    // ----------------------------------------------------------------
    using _scalars_t  = numerical_t;
    using _carrays_t  = array_t<int, 2>;
    using _seq_t      = sequential_t<int, 4>;
    using _linear_t   = linear_t<int>;
    using _map_t      = map_t<int, double>;

    using _all_io_candidates_t = typename h5::impl::cat<
        _scalars_t,
        _carrays_t,
        _seq_t,
        _linear_t,
        _map_t>::type;

    using io_types = filter_t<_all_io_candidates_t, is_supported>;

    // ----------------------------------------------------------------
    // name<T> specialisations  (in h5::test namespace)
    // ----------------------------------------------------------------

    template <class T, class...>
    struct name {
        static constexpr char const* value = "T";
    };

    // scalar arithmetic types
    template <> struct name<char>                     { static constexpr char const* value = "char"; };
    template <> struct name<unsigned char>            { static constexpr char const* value = "unsigned char"; };
    template <> struct name<short>                    { static constexpr char const* value = "short"; };
    template <> struct name<unsigned short>           { static constexpr char const* value = "unsigned short"; };
    template <> struct name<int>                      { static constexpr char const* value = "int"; };
    template <> struct name<unsigned int>             { static constexpr char const* value = "unsigned int"; };
    template <> struct name<long long int>            { static constexpr char const* value = "long long int"; };
    template <> struct name<unsigned long long int>   { static constexpr char const* value = "unsigned long long int"; };
    template <> struct name<float>                    { static constexpr char const* value = "float"; };
    template <> struct name<double>                   { static constexpr char const* value = "double"; };
    template <> struct name<pod_t>                    { static constexpr char const* value = "pod_t"; };

    // C-array ranks 1–7
    template <class T, size_t N>             struct name<T[N]>         { static constexpr char const* value = "T[i]"; };
    template <class T, size_t I, size_t J>   struct name<T[I][J]>      { static constexpr char const* value = "T[i,j]"; };
    template <class T, size_t I, size_t J, size_t K> struct name<T[I][J][K]> { static constexpr char const* value = "T[i,j,k]"; };
    template <class T, size_t A, size_t B, size_t C, size_t D>
    struct name<T[A][B][C][D]>                  { static constexpr char const* value = "T[i,j,k,l]"; };
    template <class T, size_t A, size_t B, size_t C, size_t D, size_t E>
    struct name<T[A][B][C][D][E]>               { static constexpr char const* value = "T[i,j,k,l,m]"; };
    template <class T, size_t A, size_t B, size_t C, size_t D, size_t E, size_t F>
    struct name<T[A][B][C][D][E][F]>            { static constexpr char const* value = "T[i,j,k,l,m,n]"; };
    template <class T, size_t A, size_t B, size_t C, size_t D, size_t E, size_t F, size_t G>
    struct name<T[A][B][C][D][E][F][G]>         { static constexpr char const* value = "T[i,j,k,l,m,n,o]"; };

    // std::array
    template <class T, size_t N> struct name<std::array<T,N>>                   { static constexpr char const* value = "std::array<T,N>"; };
    template <class T, size_t N> struct name<std::array<std::vector<T>,N>>       { static constexpr char const* value = "std::array<std::vector<T>,N>"; };

    // std::vector
    template <class T> struct name<std::vector<T>>                               { static constexpr char const* value = "std::vector<T>"; };
    template <class T> struct name<std::vector<std::vector<T>>>                  { static constexpr char const* value = "std::vector<std::vector<T>>"; };
    template <class T, size_t N> struct name<std::vector<std::array<T,N>>>       { static constexpr char const* value = "std::vector<std::array<T,N>>"; };

    // sequential non-contiguous
    template <class T> struct name<std::deque<T>>                                { static constexpr char const* value = "std::deque<T>"; };
    template <class T> struct name<std::forward_list<T>>                         { static constexpr char const* value = "std::forward_list<T>"; };
    template <class T> struct name<std::list<T>>                                 { static constexpr char const* value = "std::list<T>"; };

    // ordered sets / maps
    template <class T> struct name<std::set<T>>                                  { static constexpr char const* value = "std::set<T>"; };
    template <class T> struct name<std::multiset<T>>                             { static constexpr char const* value = "std::multiset<T>"; };
    template <class K, class V> struct name<std::map<K,V>>                       { static constexpr char const* value = "std::map<K,V>"; };
    template <class K, class V> struct name<std::multimap<K,V>>                  { static constexpr char const* value = "std::multimap<K,V>"; };

    // unordered sets / maps
    template <class T> struct name<std::unordered_set<T>>                        { static constexpr char const* value = "std::unordered_set<T>"; };
    template <class T> struct name<std::unordered_multiset<T>>                   { static constexpr char const* value = "std::unordered_multiset<T>"; };
    template <class K, class V> struct name<std::unordered_map<K,V>>             { static constexpr char const* value = "std::unordered_map<K,V>"; };
    template <class K, class V> struct name<std::unordered_multimap<K,V>>        { static constexpr char const* value = "std::unordered_multimap<K,V>"; };

} // namespace h5::test

// ----------------------------------------------------------------
// std::hash specialisation for pod_t so it can be used in unordered
// containers if ever needed
// ----------------------------------------------------------------
namespace std {
    template <>
    struct hash<h5::test::pod_t> {
        size_t operator()(const h5::test::pod_t& p) const noexcept {
            size_t h1 = std::hash<int>{}(p.a);
            size_t h2 = std::hash<float>{}(p.b);
            size_t h3 = std::hash<double>{}(p.c);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}
