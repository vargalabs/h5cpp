#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name itpp)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v4.3.1")
set(itpp_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
set(itpp_build_dir ${CMAKE_CURRENT_BINARY_DIR}/${current_version_dir}-build)

if(NOT EXISTS ${itpp_build_dir}/itpp/libitpp_static.a)
    file(MAKE_DIRECTORY ${itpp_build_dir})

    execute_process(
        COMMAND ${CMAKE_COMMAND}
            -S ${itpp_source_dir}
            -B ${itpp_build_dir}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_STANDARD=14
            -DCMAKE_CXX_FLAGS="-Wno-register"
            -DITPP_SHARED_LIB=OFF
            -DOLD_TESTS=OFF
            -DHTML_DOCS=OFF
        RESULT_VARIABLE itpp_configure_result
        OUTPUT_VARIABLE itpp_configure_output
        ERROR_VARIABLE itpp_configure_error
    )
    if(NOT itpp_configure_result EQUAL 0)
        message(FATAL_ERROR "itpp configure failed:\n${itpp_configure_error}")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${itpp_build_dir} --config Release --parallel 1
        RESULT_VARIABLE itpp_build_result
        OUTPUT_VARIABLE itpp_build_output
        ERROR_VARIABLE itpp_build_error
    )
    if(NOT itpp_build_result EQUAL 0)
        message(FATAL_ERROR "itpp build failed:\n${itpp_build_error}")
    endif()
endif()

add_library(${libname} STATIC IMPORTED GLOBAL)
set_target_properties(${libname} PROPERTIES
    IMPORTED_LOCATION ${itpp_build_dir}/itpp/libitpp_static.a
    INTERFACE_INCLUDE_DIRECTORIES "${itpp_source_dir};${itpp_build_dir}"
)
