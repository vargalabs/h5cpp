# SPDX-License-Identifier: MIT
# Copyright (c) 2018-2026 vargaconsulting, Toronto,ON Canada

# Imath v3.2.2 — half-precision floating-point from the OpenEXR project
# https://github.com/AcademySoftwareFoundation/Imath  License: BSD-3-Clause
#
# NOTE: half.cpp provides stream operators (operator<< / operator>>)
# and printBits() that are declared in half.h as IMATH_EXPORT symbols.
# A static convenience library is compiled so consumers do not need to
# compile half.cpp themselves.

set(name imath)
project(lib${name} LANGUAGES CXX)

set(current_version_dir "v3.2.2")
set(imath_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})

add_library(lib${name} STATIC ${imath_source_dir}/half.cpp)
add_library(h5cpp::${name} ALIAS lib${name})

target_include_directories(lib${name} PUBLIC
    $<BUILD_INTERFACE:${imath_source_dir}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(lib${name} PUBLIC cxx_std_17)

# The vendored ImathConfig.h does not enable API visibility tagging;
# no special compile definitions are needed.
