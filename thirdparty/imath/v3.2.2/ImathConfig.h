// SPDX-License-Identifier: BSD-3-Clause
// Copyright Contributors to the OpenEXR Project.

// This file is a pre-filled static version of ImathConfig.h.in for
// vendored header-only use in h5cpp. Generated cmake variables have
// been substituted with sensible cross-platform defaults.
// Upstream: https://github.com/AcademySoftwareFoundation/Imath  tag v3.2.2

#ifndef INCLUDED_IMATH_CONFIG_H
#define INCLUDED_IMATH_CONFIG_H 1

#pragma once

// Use lookup table for half-to-float conversion by default.
// The compiler can override with F16C flags or IMATH_HALF_NO_LOOKUP_TABLE.
/* #undef IMATH_HALF_USE_LOOKUP_TABLE */

// Large stack support — conservative default (disabled).
/* #undef IMATH_HAVE_LARGE_STACK */

//////////////////////
// C++ namespace configuration

#define IMATH_INTERNAL_NAMESPACE_CUSTOM 0
#define IMATH_INTERNAL_NAMESPACE Imath

#define IMATH_NAMESPACE_CUSTOM 0
#define IMATH_NAMESPACE Imath

// Version information
#define IMATH_VERSION_STRING "3.2.2"
#define IMATH_PACKAGE_STRING "Imath 3.2.2"

#define IMATH_VERSION_MAJOR 3
#define IMATH_VERSION_MINOR 2
#define IMATH_VERSION_PATCH 2
#define IMATH_VERSION_RELEASE_TYPE ""

#define IMATH_VERSION_HEX                                                      \
    ((uint32_t(IMATH_VERSION_MAJOR) << 24) |                                   \
     (uint32_t(IMATH_VERSION_MINOR) << 16) |                                   \
     (uint32_t(IMATH_VERSION_PATCH) << 8))

// Library API version (SOVERSION): 30.2.2
#define IMATH_LIB_VERSION_STRING "30.2.2"

// Use noexcept (enabled by default for modern C++17+).
#define IMATH_USE_NOEXCEPT 1
#if IMATH_USE_NOEXCEPT
#    define IMATH_NOEXCEPT noexcept
#else
#    define IMATH_NOEXCEPT
#endif

// Enable interoperability constructors/assignments by default.
#ifndef IMATH_FOREIGN_VECTOR_INTEROP
#    if defined(__GNUC__) && __GNUC__ == 4 && !defined(__clang__)
#        define IMATH_FOREIGN_VECTOR_INTEROP 0
#    else
#        define IMATH_FOREIGN_VECTOR_INTEROP 1
#    endif
#endif

// GPU (CUDA / HIP) decorator.
#if defined(__CUDACC__) || defined(__HIP__)
#    define IMATH_HOSTDEVICE __host__ __device__
#else
#    define IMATH_HOSTDEVICE
#endif

// Branch-prediction hints.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#    ifdef __cplusplus
#        define IMATH_LIKELY(x) (__builtin_expect(static_cast<bool>(x), true))
#        define IMATH_UNLIKELY(x) (__builtin_expect(static_cast<bool>(x), false))
#    else
#        define IMATH_LIKELY(x) (__builtin_expect((x), 1))
#        define IMATH_UNLIKELY(x) (__builtin_expect((x), 0))
#    endif
#else
#    define IMATH_LIKELY(x) (x)
#    define IMATH_UNLIKELY(x) (x)
#endif

#ifndef __has_attribute
#    define __has_attribute(x) 0
#endif

// Deprecation decorator.
#if defined(_MSC_VER)
#    define IMATH_DEPRECATED(msg) __declspec(deprecated(msg))
#elif defined(__cplusplus) && __cplusplus >= 201402L
#    define IMATH_DEPRECATED(msg) [[deprecated(msg)]]
#elif defined(__GNUC__) || defined(__clang__)
#    define IMATH_DEPRECATED(msg) __attribute__((deprecated(msg)))
#else
#    define IMATH_DEPRECATED(msg)
#endif

// Symbol visibility — disabled for vendored header-only use.
/* #undef IMATH_ENABLE_API_VISIBILITY */

#endif // INCLUDED_IMATH_CONFIG_H
