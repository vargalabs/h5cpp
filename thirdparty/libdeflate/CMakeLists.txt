# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

set(name libdeflate)
set(libname libh5cpp-${name})
project(${libname} LANGUAGES C)

set(current_version_dir "v1.25.0")
set(libdeflate_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})

add_library(libh5cpp-libdeflate INTERFACE)

if(H5CPP_USE_LIBDEFLATE AND EXISTS ${libdeflate_source_dir}/libdeflate.h)
    set(libdeflate_sources
        ${libdeflate_source_dir}/lib/adler32.c
        ${libdeflate_source_dir}/lib/crc32.c
        ${libdeflate_source_dir}/lib/deflate_compress.c
        ${libdeflate_source_dir}/lib/deflate_decompress.c
        ${libdeflate_source_dir}/lib/gzip_compress.c
        ${libdeflate_source_dir}/lib/gzip_decompress.c
        ${libdeflate_source_dir}/lib/utils.c
        ${libdeflate_source_dir}/lib/zlib_compress.c
        ${libdeflate_source_dir}/lib/zlib_decompress.c
    )

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$")
        list(APPEND libdeflate_sources ${libdeflate_source_dir}/lib/x86/cpu_features.c)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|armv[0-9].*)$")
        list(APPEND libdeflate_sources ${libdeflate_source_dir}/lib/arm/cpu_features.c)
    endif()

    add_library(libh5cpp-libdeflate-static STATIC ${libdeflate_sources})
    target_include_directories(libh5cpp-libdeflate-static PUBLIC
        $<BUILD_INTERFACE:${libdeflate_source_dir}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_features(libh5cpp-libdeflate-static PUBLIC c_std_99)
    target_compile_definitions(libh5cpp-libdeflate-static PUBLIC H5CPP_HAS_LIBDEFLATE=1)

    target_link_libraries(libh5cpp-libdeflate INTERFACE libh5cpp-libdeflate-static)
    target_compile_definitions(libh5cpp-libdeflate INTERFACE H5CPP_HAS_LIBDEFLATE=1)

    set(H5CPP_LIBDEFLATE_FOUND TRUE CACHE INTERNAL "H5CPP vendored libdeflate was configured" FORCE)
    set(H5CPP_LIBDEFLATE_INCLUDE_DIR ${libdeflate_source_dir} CACHE INTERNAL "H5CPP vendored libdeflate include directory" FORCE)
    set(H5CPP_LIBDEFLATE_TARGET libh5cpp-libdeflate CACHE INTERNAL "H5CPP vendored libdeflate target" FORCE)
elseif(H5CPP_REQUIRE_LIBDEFLATE)
    message(FATAL_ERROR "H5CPP_REQUIRE_LIBDEFLATE is ON, but ${libdeflate_source_dir} is missing libdeflate.h")
else()
    target_compile_definitions(libh5cpp-libdeflate INTERFACE H5CPP_DISABLE_LIBDEFLATE=1)
    set(H5CPP_LIBDEFLATE_FOUND FALSE CACHE INTERNAL "H5CPP vendored libdeflate was configured" FORCE)
endif()
