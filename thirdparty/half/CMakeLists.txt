# SPDX-License-Identifier: MIT
# Copyright (c) 2018-2026 vargaconsulting, Toronto,ON Canada

# half v2.2.0 — Christian Rau's half-precision floating-point library
# https://half.sourceforge.net/  License: MIT

set(name half)
set(libname h5cpp::${name})
project(lib${name} LANGUAGES CXX)

set(current_version_dir "v2.2.0")
add_library(lib${name} INTERFACE)
add_library(${libname} ALIAS lib${name})
target_include_directories(lib${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
