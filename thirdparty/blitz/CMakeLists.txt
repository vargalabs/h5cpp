#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name blitz)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v1.0.2")
set(blitz_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
set(blitz_build_dir ${CMAKE_CURRENT_BINARY_DIR}/${current_version_dir}-build)

if(NOT EXISTS ${blitz_build_dir}/src/libblitz.a)
    file(MAKE_DIRECTORY ${blitz_build_dir})

    execute_process(
        COMMAND ${CMAKE_COMMAND}
            -S ${blitz_source_dir}
            -B ${blitz_build_dir}
            -DCMAKE_BUILD_TYPE=Release
            -DBUILD_TESTING=OFF
            -DBUILD_DOC=OFF
        RESULT_VARIABLE blitz_configure_result
        OUTPUT_VARIABLE blitz_configure_output
        ERROR_VARIABLE blitz_configure_error
    )
    if(NOT blitz_configure_result EQUAL 0)
        message(FATAL_ERROR "Blitz++ configure failed:\n${blitz_configure_error}")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} --build ${blitz_build_dir} --config Release
        RESULT_VARIABLE blitz_build_result
        OUTPUT_VARIABLE blitz_build_output
        ERROR_VARIABLE blitz_build_error
    )
    if(NOT blitz_build_result EQUAL 0)
        message(FATAL_ERROR "Blitz++ build failed:\n${blitz_build_error}")
    endif()
endif()

add_library(${libname} STATIC IMPORTED GLOBAL)
set_target_properties(${libname} PROPERTIES
    IMPORTED_LOCATION ${blitz_build_dir}/src/libblitz.a
    INTERFACE_INCLUDE_DIRECTORIES "${blitz_build_dir};${blitz_source_dir}"
)
