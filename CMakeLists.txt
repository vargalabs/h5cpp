# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. 

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

if(CMAKE_VERSION VERSION_GREATER "3.17.0")
    cmake_policy(SET CMP0074 NEW)
endif()

project(libh5cpp-dev VERSION 1.10.4.6 LANGUAGES CXX C)


option(H5CPP_BUILD_EXAMPLES "Build examples" OFF)
option(H5CPP_BUILD_TESTS    "Build tests"    OFF)

# Generate compile_commands.json for clangd / VS Code / tooling
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

find_package(Threads REQUIRED QUIET)

# match hdf5 versioning x.x.x.h5cpp-version
set(H5CPP_VERSION
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
)

if(WIN32)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "C:/h5cpp" CACHE PATH "default install path" FORCE)
    endif()
elseif(UNIX)
    # Prefer system locations where distro HDF5 packages are commonly installed
    set(HDF5_ROOT "/usr/local" "/usr")
    find_package(PkgConfig REQUIRED QUIET)
endif()

if(APPLE)
    message(STATUS "***WARNING: h5cpp is not fully tested on apple ***")
endif()

if(ANDROID)
    message(STATUS "***WARNING: h5cpp is not fully tested on android ***")
endif()

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)

if(HDF5_VERSION VERSION_LESS ${H5CPP_VERSION})
    message(FATAL_ERROR
        "-- !!! H5CPP examples require HDF5 v${H5CPP_VERSION} or greater !!!"
    )
else()
    message(STATUS
        "H5CPP ${PROJECT_VERSION} matches with minimum required HDF5 v${HDF5_VERSION}"
    )
endif()

if(HDF5_VERSION VERSION_GREATER 1.12.0)
    message(STATUS "H5CPP KITA support is available, see webpage for details...")
endif()

# MPI is optional
find_package(MPI QUIET COMPONENTS C CXX)
if(MPI_FOUND AND HDF5_IS_PARALLEL)
    message(STATUS "MPI and PHDF5 found: Parallel H5CPP enabled")
endif()

option(H5CPP_BUILD_TESTS "Build tests" OFF)

# Header-only target
add_library(h5cpp INTERFACE)
add_library(h5cpp::h5cpp ALIAS h5cpp)

target_include_directories(h5cpp INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(h5cpp INTERFACE
    Threads::Threads
    ${HDF5_LIBRARIES}
)

target_compile_features(h5cpp INTERFACE cxx_std_17)

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

# Install header-only library
install(DIRECTORY h5cpp DESTINATION ${INSTALL_DIR_INCLUDE})
install(TARGETS h5cpp EXPORT h5cpp-targets)

# Install examples
install(DIRECTORY examples/attributes   DESTINATION "${INSTALL_DIR_EXAMPLES}")
install(DIRECTORY examples/basics       DESTINATION "${INSTALL_DIR_EXAMPLES}")
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


if(H5CPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

enable_testing()
if(H5CPP_BUILD_TESTS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt)
    add_subdirectory(test)
endif()
