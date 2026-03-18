#pragma once


#include <type_traits>

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
}