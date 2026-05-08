# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

enable_testing()
find_package(ZLIB REQUIRED)

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${HDF5_C_INCLUDE_DIRS}
)

set(LIBS
    ${HDF5_LIBRARIES}
    ZLIB::ZLIB
    libdoctest
)

## begin thirdparty library dependencies
  add_subdirectory(thirdparty EXCLUDE_FROM_ALL)
## end library dependencies

# usage:
#   add_test_case(H5cmake)
#   add_test_case(H5Iall)
#   add_test_case(H5Tmeta LIBRARIES some_extra_target)
function(add_test_case file_name)
    set(options)
    set(one_value_args)
    set(multi_value_args LIBRARIES)
    cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    string(TOLOWER "${file_name}" file_name_lc)
    set(test_name test-${file_name_lc})

    add_executable(${test_name} ${file_name}.cpp)

    set_target_properties(${test_name} PROPERTIES
        CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    target_include_directories(${test_name} PRIVATE ${INCLUDE_DIRS})
    target_link_libraries(${test_name} PRIVATE ${LIBS})

    foreach(lib IN LISTS ARG_LIBRARIES)
        target_link_libraries(${test_name} PRIVATE ${lib})
    endforeach()

    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# smoke test — validates harness builds and runs
add_test_case(example)

# test content lands in issue/94
add_test_case(H5compat)
add_test_case(H5Dgather)
add_test_case(H5meta)
add_test_case(H5Tmeta)
# add_test_case(H5Iall)
