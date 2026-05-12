# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

cmake_minimum_required(VERSION 3.22.0)
project(libh5cpp-dev VERSION 1.10.4.6 LANGUAGES CXX C)

# ─── Standard Settings ────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(WARNING "In-source builds are discouraged.")
else()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_BINARY_DIR}/compile_commands.json
            ${CMAKE_SOURCE_DIR}/compile_commands.json
        OUTPUT_QUIET ERROR_QUIET
    )
endif()

# ─── Git Version Injection ────────────────────────────────────────────────────────
execute_process(COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
execute_process(COMMAND git describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_TAG OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

set(H5CPP_BASE_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
if(GIT_TAG)
    set(H5CPP_VERSION "${GIT_TAG}")
else()
    set(H5CPP_VERSION "${H5CPP_BASE_VERSION}-${GIT_COMMIT_HASH}")
endif()

# ─── Options ──────────────────────────────────────────────────────────────────────
option(H5CPP_BUILD_EXAMPLES "Build examples" OFF)
option(H5CPP_BUILD_TESTS    "Build tests"    OFF)

# ─── System Dependencies ──────────────────────────────────────────────────────────
include(GNUInstallDirs)

find_package(Threads REQUIRED QUIET)
find_package(Python3 COMPONENTS Interpreter QUIET)

# match hdf5 versioning x.x.x.h5cpp-version
set(H5CPP_VERSION
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
)

if(WIN32)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "C:/h5cpp" CACHE PATH "default install path" FORCE)
    endif()
elseif(UNIX)
    set(HDF5_ROOT "/usr/local" "/usr")
    find_package(PkgConfig REQUIRED QUIET)
endif()

if(APPLE)
    message(STATUS "H5CPP macOS support: CI runs on macos-13 with Homebrew HDF5 and LLVM clang-17/18")
endif()

if(ANDROID)
    message(STATUS "***WARNING: h5cpp is not fully tested on android ***")
endif()

find_package(Threads REQUIRED QUIET)

# ─── HDF5 ─────────────────────────────────────────────────────────────────────────
find_package(HDF5 REQUIRED COMPONENTS C)

if(HDF5_VERSION VERSION_LESS ${H5CPP_BASE_VERSION})
    message(FATAL_ERROR
        "-- !!! H5CPP examples require HDF5 v${H5CPP_BASE_VERSION} or greater !!!"
    )
else()
    message(STATUS
        "H5CPP ${PROJECT_VERSION} matches with minimum required HDF5 v${HDF5_VERSION}"
    )
endif()

if(HDF5_VERSION VERSION_GREATER 1.12.0)
    message(STATUS "H5CPP KITA support is available, see webpage for details...")
endif()

# ─── MPI (optional) ───────────────────────────────────────────────────────────────
find_package(MPI QUIET COMPONENTS C CXX)
if(MPI_FOUND AND HDF5_IS_PARALLEL)
    message(STATUS "MPI and PHDF5 found: Parallel H5CPP enabled")
endif()

# ─── Include Directories & Libraries ──────────────────────────────────────────────
list(APPEND H5CPP_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${HDF5_INCLUDE_DIRS}
)

list(APPEND H5CPP_LIBS
    Threads::Threads
    ${HDF5_LIBRARIES}
)

# ─── Thirdparty ───────────────────────────────────────────────────────────────────
add_subdirectory(thirdparty EXCLUDE_FROM_ALL)

# ─── Header-only Library Target ───────────────────────────────────────────────────
add_library(h5cpp INTERFACE)
add_library(h5cpp::h5cpp ALIAS h5cpp)

target_include_directories(h5cpp INTERFACE
    $<BUILD_INTERFACE:${H5CPP_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(h5cpp INTERFACE ${H5CPP_LIBS})
target_compile_features(h5cpp INTERFACE cxx_std_17)

# ─── Subdirectories ───────────────────────────────────────────────────────────────
if(H5CPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

enable_testing()
if(H5CPP_BUILD_TESTS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt)
    add_subdirectory(test)
endif()

# ─── Install Rules ────────────────────────────────────────────────────────────────
if(WIN32)
    set(INSTALL_DIR_CMAKE    "${CMAKE_INSTALL_PREFIX}/share/h5cpp/cmake")
    set(INSTALL_DIR_INCLUDE  "${CMAKE_INSTALL_PREFIX}/include")
    set(INSTALL_DIR_EXAMPLES "${CMAKE_INSTALL_PREFIX}/examples")
elseif(UNIX)
    set(INSTALL_DIR_CMAKE    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/h5cpp")
    set(INSTALL_DIR_INCLUDE  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
    set(INSTALL_DIR_EXAMPLES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/h5cpp")
endif()

message(STATUS "Install directory: ${INSTALL_DIR_INCLUDE}")

# Generate single-header amalgamation during build
set(H5CPP_AMALGAMATED_HEADER ${CMAKE_BINARY_DIR}/h5cpp.hpp)
add_custom_command(
    OUTPUT ${H5CPP_AMALGAMATED_HEADER}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/amalgamate.py
            ${CMAKE_SOURCE_DIR}/h5cpp ${H5CPP_AMALGAMATED_HEADER}
    DEPENDS ${CMAKE_SOURCE_DIR}/scripts/amalgamate.py
            ${CMAKE_SOURCE_DIR}/h5cpp/core
            ${CMAKE_SOURCE_DIR}/h5cpp/io
    COMMENT "Amalgamating h5cpp headers into h5cpp.hpp"
)
add_custom_target(h5cpp-amalgamate ALL DEPENDS ${H5CPP_AMALGAMATED_HEADER})

# Install header-only library (directory + single-header)
install(DIRECTORY h5cpp DESTINATION ${INSTALL_DIR_INCLUDE})
install(FILES ${H5CPP_AMALGAMATED_HEADER} DESTINATION ${INSTALL_DIR_INCLUDE})
install(TARGETS h5cpp EXPORT h5cpp-targets)

# Install examples
install(DIRECTORY examples/attributes   DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/basics       DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY thirdparty/           DESTINATION "${INSTALL_DIR_EXAMPLES}/thirdparty")
install(DIRECTORY examples/compound     DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/csv          DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/datasets     DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/datatypes    DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/half-float   DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/multi-tu     DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/optimized    DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/packet-table DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/raw_memory   DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/stl          DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/string       DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/transform    DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/utf          DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(FILES     examples/CMakeLists.txt DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(FILES     examples/Makefile       DESTINATION "${INSTALL_DIR_EXAMPLES}")

# Package config
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/config.cmake.in
    ${CMAKE_BINARY_DIR}/h5cpp-config.cmake
    INSTALL_DESTINATION ${INSTALL_DIR_CMAKE}
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/h5cpp-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(PACKAGE h5cpp)
export(EXPORT h5cpp-targets)

install(EXPORT h5cpp-targets DESTINATION ${INSTALL_DIR_CMAKE})
install(FILES
    ${CMAKE_BINARY_DIR}/h5cpp-config.cmake
    ${CMAKE_BINARY_DIR}/h5cpp-config-version.cmake
    DESTINATION ${INSTALL_DIR_CMAKE}
)
