# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Example-derived I/O round-trip tests

# Each file here is a standalone doctest executable converted from an example
function(add_example_test file_name)
    set(test_name example-test-${file_name})
    add_executable(${test_name} ${file_name}.cpp)
    set_target_properties(${test_name} PROPERTIES
        CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    target_include_directories(${test_name} PRIVATE ${INCLUDE_DIRS})
    target_link_libraries(${test_name} PRIVATE ${LIBS})
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_example_test(test_string_io)
add_example_test(test_stl_vector_io)
add_example_test(test_basics_io)
add_example_test(test_raw_memory_io)
add_example_test(test_csv_io)
