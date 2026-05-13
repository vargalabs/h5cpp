/*
 * Copyright (c) 2018-2026 vargaconsulting, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 *
 * Typed random data generator for the h5cpp test suite.
 * Ported and adapted from HDFGroup-mailinglist/tuple-2022-aprl-06/H5Uall.hpp.
 * Issue #117.
 *
 * Interface: h5::test::fill(T& obj, lower, upper, min, max)
 *   lower/upper — value range for scalar elements
 *   min/max     — size range for dynamically-sized containers
 *
 * fill() rather than get() handles C arrays uniformly: functions cannot
 * return arrays in C++, but can fill them via reference.
 */
#pragma once

#include <array>
#include <vector>
#include <deque>
#include <list>
#include <forward_list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <string>
#include <type_traits>
#include <algorithm>
#include <random>
#include "types.hpp"

namespace h5::test {
namespace detail {

    // Detection idioms not present in the current h5cpp headers
    template <class T, class = void> struct has_reserve       : std::false_type {};
    template <class T, class = void> struct has_emplace_back  : std::false_type {};
    template <class T, class = void> struct has_emplace_front : std::false_type {};
    template <class T, class = void> struct has_emplace       : std::false_type {};

    template <class T> struct has_reserve<T,
        std::void_t<decltype(std::declval<T&>().reserve(size_t{}))>>
        : std::true_type {};
    template <class T> struct has_emplace_back<T,
        std::void_t<decltype(std::declval<T&>().emplace_back(
            std::declval<typename T::value_type>()))>>
        : std::true_type {};
    template <class T> struct has_emplace_front<T,
        std::void_t<decltype(std::declval<T&>().emplace_front(
            std::declval<typename T::value_type>()))>>
        : std::true_type {};
    template <class T> struct has_emplace<T,
        std::void_t<decltype(std::declval<T&>().emplace(
            std::declval<typename T::value_type>()))>>
        : std::true_type {};

    template <class T>
    inline constexpr bool is_stl_container_v =
        has_emplace_back<T>::value ||
        has_emplace_front<T>::value ||
        has_emplace<T>::value;

    // Uniform random in [lo, hi]
    template <class T>
    T rng_uniform(T lo, T hi) {
        using dist_t = std::conditional_t<std::is_integral<T>::value,
            std::uniform_int_distribution<std::conditional_t<sizeof(T) < sizeof(int), int, T>>,
            std::uniform_real_distribution<T>>;
        thread_local std::mt19937_64 eng{std::random_device{}()};
        dist_t dist(static_cast<typename dist_t::result_type>(lo), static_cast<typename dist_t::result_type>(hi));
        return static_cast<T>(dist(eng));
    }

    inline size_t rng_size(size_t mn, size_t mx) {
        thread_local std::mt19937_64 eng{std::random_device{}()};
        return std::uniform_int_distribution<size_t>{mn, mx}(eng);
    }

    inline std::string rng_string(size_t mn, size_t mx) {
        static constexpr char alpha[] =
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        thread_local std::mt19937_64 eng{std::random_device{}()};
        size_t n = std::uniform_int_distribution<size_t>{mn, mx}(eng);
        std::uniform_int_distribution<size_t> cd{0, sizeof(alpha) - 2};
        std::string s;
        s.reserve(n);
        for (size_t i = 0; i < n; ++i) s += alpha[cd(eng)];
        return s;
    }

} // namespace detail


// ----------------------------------------------------------------
// C-array overloads — must be separate signatures (can't return arrays)
// ----------------------------------------------------------------
template <class T, size_t N>
void fill(T (&arr)[N], size_t lower, size_t upper, size_t min = 5, size_t max = 10);

template <class T, size_t N, size_t M>
void fill(T (&arr)[N][M], size_t lower, size_t upper, size_t min = 5, size_t max = 10);

template <class T, size_t N, size_t M, size_t P>
void fill(T (&arr)[N][M][P], size_t lower, size_t upper, size_t min = 5, size_t max = 10);


// ----------------------------------------------------------------
// Primary dispatcher — handles everything that isn't a C array
// ----------------------------------------------------------------
template <class T>
void fill(T& obj, size_t lower, size_t upper, size_t min = 5, size_t max = 10) {
    if constexpr (std::is_same_v<T, bool>) {
        thread_local std::mt19937_64 eng{std::random_device{}()};
        obj = std::bernoulli_distribution{
            static_cast<double>(lower) / static_cast<double>(upper)}(eng);

    } else if constexpr (std::is_arithmetic_v<T>) {
        obj = detail::rng_uniform<T>(static_cast<T>(lower), static_cast<T>(upper));

    } else if constexpr (std::is_same_v<T, std::string>) {
        obj = detail::rng_string(lower, upper);

    } else if constexpr (std::is_same_v<T, pod_t>) {
        obj.a = detail::rng_uniform<int>(static_cast<int>(lower), static_cast<int>(upper));
        obj.b = detail::rng_uniform<float>(static_cast<float>(lower), static_cast<float>(upper));
        obj.c = detail::rng_uniform<double>(static_cast<double>(lower), static_cast<double>(upper));

    } else if constexpr (std::is_array_v<T>) {
        // multidimensional C array passed as non-reference (shouldn't happen, but safe)
        for (auto& elem : obj)
            fill(elem, lower, upper, min, max);

    } else if constexpr (detail::is_stl_container_v<T>) {
        using elem_t = typename T::value_type;
        size_t n = detail::rng_size(min, max);
        obj.clear();
        if constexpr (detail::has_reserve<T>::value)
            obj.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            elem_t e{};
            fill(e, lower, upper, min, max);
            if constexpr (detail::has_emplace_back<T>::value)
                obj.emplace_back(std::move(e));
            else if constexpr (detail::has_emplace_front<T>::value)
                obj.emplace_front(std::move(e));
            else
                obj.emplace(std::move(e));
        }
    }
    // else: leave default-initialized (covers unknown POD structs)
}


// ----------------------------------------------------------------
// C-array definitions (call back into primary dispatcher)
// ----------------------------------------------------------------
template <class T, size_t N>
void fill(T (&arr)[N], size_t lower, size_t upper, size_t min, size_t max) {
    for (auto& elem : arr) fill(elem, lower, upper, min, max);
}

template <class T, size_t N, size_t M>
void fill(T (&arr)[N][M], size_t lower, size_t upper, size_t min, size_t max) {
    for (auto& row : arr) fill(row, lower, upper, min, max);
}

template <class T, size_t N, size_t M, size_t P>
void fill(T (&arr)[N][M][P], size_t lower, size_t upper, size_t min, size_t max) {
    for (auto& plane : arr) fill(plane, lower, upper, min, max);
}


// ----------------------------------------------------------------
// std::array<V,N> — explicit overload so fill(array, ...) calls
// the primary dispatcher on each element (not the C-array path)
// ----------------------------------------------------------------
template <class V, size_t N>
void fill(std::array<V,N>& arr, size_t lower, size_t upper, size_t min = 5, size_t max = 10) {
    for (auto& elem : arr) fill(elem, lower, upper, min, max);
}


// ----------------------------------------------------------------
// std::pair<K,V> — map value_type is pair<const K, V>; const_cast
// on key is safe because the pair itself is non-const stack storage
// ----------------------------------------------------------------
template <class K, class V>
void fill(std::pair<K,V>& p, size_t lower, size_t upper, size_t min = 5, size_t max = 10) {
    fill(const_cast<std::remove_const_t<K>&>(p.first),  lower, upper, min, max);
    fill(p.second,                                       lower, upper, min, max);
}

} // namespace h5::test
