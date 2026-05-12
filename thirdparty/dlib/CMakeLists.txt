#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name dlib)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v19.24")
set(dlib_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}/dlib)
set(dlib_build_dir ${CMAKE_CURRENT_BINARY_DIR}/${current_version_dir}-build)

# Platform-dependent library path (Unix .a vs Windows .lib, with multi-config Release/ subdir)
if(WIN32)
    # dlib on MSVC adds a release postfix (e.g. dlib19_24_release_64bit_msvc1939.lib);
    # start with the plain name and fall back to a glob after the build.
    set(dlib_lib_path "${dlib_build_dir}/Release/dlib.lib")
else()
    set(dlib_lib_path "${dlib_build_dir}/libdlib.a")
endif()

if(NOT EXISTS ${dlib_lib_path})
    file(MAKE_DIRECTORY ${dlib_build_dir})

    # Determine warning suppression for newer Clang versions (dlib v19.24 uses
    # '::template go(...)' without template args, which Clang 17+ errors on)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(dlib_extra_cxx_flags "-Wno-missing-template-arg-list-after-template-kw")
    else()
        set(dlib_extra_cxx_flags "")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND}
            -S ${dlib_source_dir}
            -B ${dlib_build_dir}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_STANDARD=14
            -DCMAKE_CXX_FLAGS=${dlib_extra_cxx_flags}
            -DDLIB_ISO_CPP_ONLY=OFF
            -DDLIB_NO_GUI_SUPPORT=ON
            -DDLIB_JPEG_SUPPORT=OFF
            -DDLIB_PNG_SUPPORT=OFF
            -DDLIB_GIF_SUPPORT=OFF
            -DDLIB_WEBP_SUPPORT=OFF
            -DDLIB_LINK_WITH_SQLITE3=OFF
            -DDLIB_USE_BLAS=OFF
            -DDLIB_USE_LAPACK=OFF
            -DDLIB_USE_CUDA=OFF
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        RESULT_VARIABLE dlib_configure_result
        OUTPUT_VARIABLE dlib_configure_output
        ERROR_VARIABLE dlib_configure_error
    )
    if(NOT dlib_configure_result EQUAL 0)
        message(FATAL_ERROR "dlib configure failed:\n${dlib_configure_error}")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${dlib_build_dir} --config Release
        RESULT_VARIABLE dlib_build_result
        OUTPUT_VARIABLE dlib_build_output
        ERROR_VARIABLE dlib_build_error
    )
    if(NOT dlib_build_result EQUAL 0)
        message(FATAL_ERROR "dlib build failed:\n${dlib_build_error}")
    endif()

    # Re-discover library path after build (Windows MSVC adds a postfix)
    if(WIN32)
        if(NOT EXISTS ${dlib_lib_path})
            file(GLOB dlib_lib_candidates "${dlib_build_dir}/Release/dlib*.lib")
            if(dlib_lib_candidates)
                list(GET dlib_lib_candidates 0 dlib_lib_path)
            endif()
        endif()
    endif()
endif()

add_library(${libname} STATIC IMPORTED GLOBAL)
set_target_properties(${libname} PROPERTIES
    IMPORTED_LOCATION ${dlib_lib_path}
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}
    INTERFACE_COMPILE_DEFINITIONS "DLIB_NO_GUI_SUPPORT;DLIB_DISABLE_ASSERTS;DLIB__CMAKE_GENERATED_A_CONFIG_H_FILE;DLIB_CHECK_FOR_VERSION_MISMATCH=DLIB_VERSION_MISMATCH_CHECK__EXPECTED_VERSION_19_24_0"
)
