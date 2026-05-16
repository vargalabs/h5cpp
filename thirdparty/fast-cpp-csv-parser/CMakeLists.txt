# SPDX-License-Identifier: MIT
# Copyright (c) 2018-2026 vargaconsulting, Toronto,ON Canada

# fast-cpp-csv-parser — Ben Strasser's header-only CSV parser
# https://github.com/ben-strasser/fast-cpp-csv-parser  License: BSD-3-Clause
# Vendored at commit 574a9fe4 (master, 2025-era)

set(name csv)
project(libh5cpp-${name} LANGUAGES CXX)

set(current_version_dir "574a9fe4")
add_library(libh5cpp-${name} INTERFACE)
add_library(h5cpp::${name} ALIAS libh5cpp-${name})
target_include_directories(libh5cpp-${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
