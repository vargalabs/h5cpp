# Install script for directory: /home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/random" TYPE FILE FILES
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/F.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/beta.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/chisquare.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/default.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/discrete-uniform.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/exponential.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/gamma.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/mt.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/mtparam.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/normal.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/random/uniform.h"
    )
endif()

