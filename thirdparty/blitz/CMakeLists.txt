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

# Platform-dependent library path (Unix .a vs Windows import .lib, with multi-config Release/ subdir)
if(WIN32)
    # On Windows blitz only builds SHARED; the import library is blitz.lib
    set(blitz_lib_path "${blitz_build_dir}/src/Release/blitz.lib")
else()
    set(blitz_lib_path "${blitz_build_dir}/src/libblitz.a")
endif()

if(NOT EXISTS ${blitz_lib_path})
    file(MAKE_DIRECTORY ${blitz_build_dir})

    # Windows-only: force absolute output dirs (default IMPLIB resolution under
    # WINDOWS_EXPORT_ALL_SYMBOLS lands at drive root) and enable symbol export
    # so MSVC actually produces the import .lib alongside the .dll.
    set(blitz_extra_args)
    if(WIN32)
        list(APPEND blitz_extra_args
            -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE=${blitz_build_dir}/src/Release
            -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE=${blitz_build_dir}/src/Release
            -DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE=${blitz_build_dir}/src/Release
        )
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND}
            -S ${blitz_source_dir}
            -B ${blitz_build_dir}
            -DCMAKE_BUILD_TYPE=Release
            -DBUILD_TESTING=OFF
            -DBUILD_DOC=OFF
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5
            ${blitz_extra_args}
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

    # Verify the import/static lib was produced; fail loudly with diagnostics if not.
    if(NOT EXISTS "${blitz_lib_path}")
        message(FATAL_ERROR "Blitz++ library not found after build at ${blitz_lib_path}.\n"
            "Build output:\n${blitz_build_output}\n${blitz_build_error}")
    endif()
endif()

add_library(${libname} STATIC IMPORTED GLOBAL)
set_target_properties(${libname} PROPERTIES
    IMPORTED_LOCATION ${blitz_lib_path}
    INTERFACE_INCLUDE_DIRECTORIES "${blitz_build_dir};${blitz_source_dir}"
)
