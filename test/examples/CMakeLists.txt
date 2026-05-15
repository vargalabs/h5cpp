# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Example-derived I/O round-trip tests

# Each file here is a standalone doctest executable converted from an example
function(add_example_test file_name)
    set(options)
    set(one_value_args)
    set(multi_value_args LIBRARIES)
    cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    set(test_name example-test-${file_name})
    add_executable(${test_name} ${file_name}.cpp)
    set_target_properties(${test_name} PROPERTIES
        CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    target_include_directories(${test_name} PRIVATE ${INCLUDE_DIRS})
    target_link_libraries(${test_name} PRIVATE ${LIBS})

    foreach(lib IN LISTS ARG_LIBRARIES)
        target_link_libraries(${test_name} PRIVATE ${lib})
    endforeach()

    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_example_test(test_string_io)
add_example_test(test_stl_vector_io)
add_example_test(test_basics_io)
add_example_test(test_raw_memory_io)
add_example_test(test_csv_io)
add_example_test(test_compound_struct_io)
add_example_test(test_attributes_io)
add_example_test(test_two_bit_io)
add_example_test(test_n_bit_io)

# Category C: Linear Algebra (requires thirdparty libs)
add_example_test(test_arma_io     LIBRARIES libarmadillo)
# NOTE: eigen3 test disabled — including <Eigen/Dense> before <h5cpp/all>
# triggers an HDF5 "not a datatype" error even in the original example.
# add_example_test(test_eigen3_io   LIBRARIES libeigen3)
add_example_test(test_blaze_io    LIBRARIES libblaze)
add_example_test(test_dlib_io     LIBRARIES libdlib)

# ublas v1.74 uses std::allocator::construct which was removed in C++20.
# Skip the ublas test under C++20 builds.
if(CMAKE_CXX_STANDARD LESS 20)
    find_path(H5CPP_BOOST_CONFIG_INCLUDE NAMES boost/config.hpp)
    if(H5CPP_BOOST_CONFIG_INCLUDE)
        add_example_test(test_ublas_io LIBRARIES libublas)
    else()
        message(STATUS "Boost config headers not found — skipping test_ublas_io")
    endif()
else()
    message(STATUS "C++20 build — skipping test_ublas_io (ublas v1.74 requires std::allocator::construct, removed in C++20)")
endif()

add_example_test(test_blitz_io    LIBRARIES libblitz)

# Category D: Advanced Features
add_example_test(test_packet_table_io LIBRARIES libarmadillo)
add_example_test(test_optimized_io    LIBRARIES libarmadillo)
add_example_test(test_transform_io    LIBRARIES libarmadillo)
add_example_test(test_utf_io          LIBRARIES libarmadillo)

# Category F: Multi-TU / Compiler-generated
add_example_test(test_multi_tu        LIBRARIES libarmadillo)
