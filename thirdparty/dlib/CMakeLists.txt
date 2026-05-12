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

if(NOT EXISTS ${dlib_build_dir}/libdlib.a)
    file(MAKE_DIRECTORY ${dlib_build_dir})

    execute_process(
        COMMAND ${CMAKE_COMMAND}
            -S ${dlib_source_dir}
            -B ${dlib_build_dir}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_STANDARD=14
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
endif()

add_library(${libname} STATIC IMPORTED GLOBAL)
set_target_properties(${libname} PROPERTIES
    IMPORTED_LOCATION ${dlib_build_dir}/libdlib.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}
    INTERFACE_COMPILE_DEFINITIONS "DLIB_NO_GUI_SUPPORT;DLIB_DISABLE_ASSERTS;DLIB__CMAKE_GENERATED_A_CONFIG_H_FILE;DLIB_CHECK_FOR_VERSION_MISMATCH=DLIB_VERSION_MISMATCH_CHECK__EXPECTED_VERSION_19_24_0"
)
