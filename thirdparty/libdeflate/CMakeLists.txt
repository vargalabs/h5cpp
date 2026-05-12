# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

add_library(libh5cpp-libdeflate INTERFACE)

if(H5CPP_LIBDEFLATE_FOUND)
    target_include_directories(libh5cpp-libdeflate INTERFACE ${LIBDEFLATE_INCLUDE_DIR})
    target_link_libraries(libh5cpp-libdeflate INTERFACE ${LIBDEFLATE_LIBRARY})
    target_compile_definitions(libh5cpp-libdeflate INTERFACE H5CPP_HAS_LIBDEFLATE=1)
else()
    target_compile_definitions(libh5cpp-libdeflate INTERFACE H5CPP_DISABLE_LIBDEFLATE=1)
endif()
