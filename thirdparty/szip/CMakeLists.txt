# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

set(name szip)
set(libname libh5cpp-${name})
project(${libname} LANGUAGES C)

set(current_version_dir "v2.2.0")
set(szip_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}/src)

add_library(libh5cpp-szip INTERFACE)

if(H5CPP_USE_SZIP AND EXISTS ${szip_source_dir}/szlib.h)
    set(szip_sources
        ${szip_source_dir}/encoding.c
        ${szip_source_dir}/rice.c
        ${szip_source_dir}/sz_api.c
    )

    add_library(libh5cpp-szip-static STATIC ${szip_sources})

    # SZconfig.h lives one level above the versioned src dir (our static replacement
    # for the autoconf-generated header).  List it before the src dir so it is found
    # in preference to any SZconfig.h.in-derived artefact.
    target_include_directories(libh5cpp-szip-static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${szip_source_dir}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_features(libh5cpp-szip-static PUBLIC c_std_99)
    target_compile_definitions(libh5cpp-szip-static PUBLIC H5CPP_HAS_SZIP=1)

    # Suppress harmless warnings from the vendored C89-era source.
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(libh5cpp-szip-static PRIVATE
            -Wno-implicit-function-declaration
            -Wno-unused-variable
            -Wno-sign-compare)
    endif()

    target_link_libraries(libh5cpp-szip INTERFACE libh5cpp-szip-static)
    target_compile_definitions(libh5cpp-szip INTERFACE H5CPP_HAS_SZIP=1)

    set(H5CPP_SZIP_FOUND TRUE CACHE INTERNAL "H5CPP vendored szip was configured" FORCE)
    set(H5CPP_SZIP_INCLUDE_DIR ${szip_source_dir} CACHE INTERNAL "H5CPP vendored szip include directory" FORCE)
    set(H5CPP_SZIP_TARGET libh5cpp-szip CACHE INTERNAL "H5CPP vendored szip target" FORCE)
else()
    set(H5CPP_SZIP_FOUND FALSE CACHE INTERNAL "H5CPP vendored szip was configured" FORCE)
endif()
