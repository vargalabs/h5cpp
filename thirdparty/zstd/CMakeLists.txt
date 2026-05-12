# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

set(name zstd)
set(libname libh5cpp-${name})
project(${libname} LANGUAGES C)

set(current_version_dir "v1.5.7")
set(zstd_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
set(zstd_lib_dir    ${zstd_source_dir}/lib)

add_library(libh5cpp-zstd INTERFACE)

if(H5CPP_USE_ZSTD AND EXISTS ${zstd_lib_dir}/zstd.h)
    set(zstd_sources
        # common
        ${zstd_lib_dir}/common/debug.c
        ${zstd_lib_dir}/common/entropy_common.c
        ${zstd_lib_dir}/common/error_private.c
        ${zstd_lib_dir}/common/fse_decompress.c
        ${zstd_lib_dir}/common/pool.c
        ${zstd_lib_dir}/common/threading.c
        ${zstd_lib_dir}/common/xxhash.c
        ${zstd_lib_dir}/common/zstd_common.c
        # compress
        ${zstd_lib_dir}/compress/fse_compress.c
        ${zstd_lib_dir}/compress/hist.c
        ${zstd_lib_dir}/compress/huf_compress.c
        ${zstd_lib_dir}/compress/zstd_compress.c
        ${zstd_lib_dir}/compress/zstd_compress_literals.c
        ${zstd_lib_dir}/compress/zstd_compress_sequences.c
        ${zstd_lib_dir}/compress/zstd_compress_superblock.c
        ${zstd_lib_dir}/compress/zstd_double_fast.c
        ${zstd_lib_dir}/compress/zstd_fast.c
        ${zstd_lib_dir}/compress/zstd_lazy.c
        ${zstd_lib_dir}/compress/zstd_ldm.c
        ${zstd_lib_dir}/compress/zstd_opt.c
        ${zstd_lib_dir}/compress/zstd_preSplit.c
        ${zstd_lib_dir}/compress/zstdmt_compress.c
        # decompress
        ${zstd_lib_dir}/decompress/huf_decompress.c
        ${zstd_lib_dir}/decompress/zstd_ddict.c
        ${zstd_lib_dir}/decompress/zstd_decompress.c
        ${zstd_lib_dir}/decompress/zstd_decompress_block.c
    )

    # x86-64 optimised decompressor (assembly; optional, falls back gracefully)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
        if(EXISTS ${zstd_lib_dir}/decompress/huf_decompress_amd64.S)
            enable_language(ASM)
            list(APPEND zstd_sources ${zstd_lib_dir}/decompress/huf_decompress_amd64.S)
        endif()
    endif()

    add_library(libh5cpp-zstd-static STATIC ${zstd_sources})
    target_include_directories(libh5cpp-zstd-static PUBLIC
        $<BUILD_INTERFACE:${zstd_lib_dir}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_features(libh5cpp-zstd-static PUBLIC c_std_99)
    target_compile_definitions(libh5cpp-zstd-static PUBLIC H5CPP_HAS_ZSTD=1)
    # Keep zstd single-threaded within our library; users can enable MT externally
    target_compile_definitions(libh5cpp-zstd-static PRIVATE ZSTD_MULTITHREAD=0)

    target_link_libraries(libh5cpp-zstd INTERFACE libh5cpp-zstd-static)
    target_compile_definitions(libh5cpp-zstd INTERFACE H5CPP_HAS_ZSTD=1)

    set(H5CPP_ZSTD_FOUND TRUE CACHE INTERNAL "H5CPP vendored zstd was configured" FORCE)
    set(H5CPP_ZSTD_INCLUDE_DIR ${zstd_lib_dir} CACHE INTERNAL "H5CPP vendored zstd include directory" FORCE)
    set(H5CPP_ZSTD_TARGET libh5cpp-zstd CACHE INTERNAL "H5CPP vendored zstd target" FORCE)
else()
    set(H5CPP_ZSTD_FOUND FALSE CACHE INTERNAL "H5CPP vendored zstd was configured" FORCE)
endif()
