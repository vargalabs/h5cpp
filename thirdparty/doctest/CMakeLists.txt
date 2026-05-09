# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. 

set(name doctest)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v2.4.11")
add_library(${libname} INTERFACE)
target_include_directories(${libname} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}/include)

