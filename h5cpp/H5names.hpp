/* SPDX-License-Identifier: MIT
 * This file is part of H5CPP.
 * Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. */

#pragma once
#include <cstddef>
#include <complex> 
#include "H5meta.hpp"

namespace h5 {
    template <class T, int N> struct name<T[N]> {
        static constexpr const char* value = meta::extent_to_string<T[N]>::value;
    };
}
// half float support: 
// TODO: factor out in a separate file
#ifdef HALF_HALF_HPP
    H5CPP_REGISTER_HALF_FLOAT(half_float::half)
#endif
// Open XDR doesn-t define namespace or 
#ifdef WITH_OPENEXR_HALF 
    H5CPP_REGISTER_HALF_FLOAT(OPENEXR_NAMESPACE::half)
#endif

namespace h5 {

    //template <> struct name<char const*> { static constexpr char const * value = "char const*"; };
    //template <> struct name<char**> { static constexpr char const * value = "char**"; };
    template <> struct name<short**> { static constexpr char const * value = "short**"; };
    template <> struct name<int**> { static constexpr char const * value = "int**"; };
    template <> struct name<long**> { static constexpr char const * value = "long**"; };
    template <> struct name<long long **> { static constexpr char const * value = "long long**"; };
    template <> struct name<float**> { static constexpr char const * value = "float**"; };
    template <> struct name<double**> { static constexpr char const * value = "double**"; };
    template <> struct name<void**> { static constexpr char const * value = "void**"; };
    template <class T> struct name<std::complex<T>> { static constexpr char const * value = "complex<T>"; };
    template <> struct name<std::complex<float>> { static constexpr char const * value = "complex<float>"; };
    template <> struct name<std::complex<double>> { static constexpr char const * value = "complex<double>"; };
    template <> struct name<std::complex<short>> { static constexpr char const * value = "complex<short>"; };
    template <> struct name<std::complex<int>> { static constexpr char const * value = "complex<int>"; };
}