/* SPDX-License-Identifier: MIT
 * This file is part of H5CPP.
 * Copyright (c) 2018-2024 Varga Consulting, Toronto, ON, Canada.
 * Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. */

#pragma once

#include <hdf5.h>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

namespace h5::compat {
    template<std::size_t...> struct index_sequence {};
    template<std::size_t n, std::size_t... next>
    struct index_sequence_impl_t : index_sequence_impl_t<n - 1U, n - 1U, next...> {};
    template<std::size_t... next> struct index_sequence_impl_t<0U, next...> {
        using type = index_sequence<next...>;
    };
    template<std::size_t N> using make_index_sequence = typename index_sequence_impl_t<N>::type;
    template<class T, class tuple_t, std::size_t... i>
    constexpr herr_t apply_impl(T&& f, tuple_t&& t, index_sequence<i...>) {
        return std::forward<T>(f)(std::get<i>(std::forward<tuple_t>(t))...);
    }

    template<class T, class tuple_t>
    constexpr herr_t apply(T&& f, tuple_t&& t) {
        using tp_t = typename std::decay<tuple_t>::type;
        return apply_impl(
            std::forward<T>(f),
            std::forward<tuple_t>(t),
            make_index_sequence<std::tuple_size<tp_t>::value>{}
        );
    }
    struct nonesuch {
        nonesuch() = delete;
        ~nonesuch() = delete;
        nonesuch(nonesuch const&) = delete;
        void operator=(nonesuch const&) = delete;
    };
    namespace detail {
        template<class default_t, class always_void_t, template<class...> class op_t, class... args_t>
        struct detector_t {
            using value_t = std::false_type;
            using type = default_t;
        };

        template<class default_t, template<class...> class op_t, class... args_t>
        struct detector_t<default_t, std::void_t<op_t<args_t...>>, op_t, args_t...> {
            using value_t = std::true_type;
            using type = op_t<args_t...>;
        };

    } // namespace detail

    template<template<class...> class op_t, class... args_t> using is_detected = typename detail::detector_t<nonesuch, void, op_t, args_t...>::value_t;
    template<template<class...> class op_t, class... args_t> using detected_t = typename detail::detector_t<nonesuch, void, op_t, args_t...>::type;
    template<class default_t, template<class...> class op_t, class... args_t> using detected_or = detail::detector_t<default_t, void, op_t, args_t...>;
    template<template<class...> class op_t, class... args_t> constexpr bool is_detected_v = is_detected<op_t, args_t...>::value;
    template<class default_t, template<class...> class op_t, class... args_t> using detected_or_t = typename detected_or<default_t, op_t, args_t...>::type;
    template<class expected_t, template<class...> class op_t, class... args_t> using is_detected_exact = std::is_same<expected_t, detected_t<op_t, args_t...>>;
    template<class expected_t, template<class...> class op_t, class... args_t> constexpr bool is_detected_exact_v = is_detected_exact<expected_t, op_t, args_t...>::value;
    template<class to_t, template<class...> class op_t, class... args_t> using is_detected_convertible = std::is_convertible<detected_t<op_t, args_t...>, to_t>;
    template<class to_t, template<class...> class op_t, class... args_t> constexpr bool is_detected_convertible_v = is_detected_convertible<to_t, op_t, args_t...>::value;

    template<class T> using is_pod = std::is_pod<typename std::remove_cv<typename std::remove_reference<T>::type>::type>;
    template<class T> constexpr bool is_pod_v = is_pod<T>::value;      
} // namespace h5::compat

namespace h5::meta::compat {
    using h5::compat::detected_or;
    using h5::compat::detected_or_t;
    using h5::compat::detected_t;
    using h5::compat::is_detected;
    using h5::compat::is_detected_convertible;
    using h5::compat::is_detected_convertible_v;
    using h5::compat::is_detected_exact;
    using h5::compat::is_detected_exact_v;
    using h5::compat::is_detected_v;
    using h5::compat::is_pod;
    using h5::compat::is_pod_v;    
    using h5::compat::nonesuch;
}