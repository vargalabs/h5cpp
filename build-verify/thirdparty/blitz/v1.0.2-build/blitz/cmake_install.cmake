# Install script for directory: /home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/blitz" TYPE FILE FILES
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/array-impl.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/array.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/array.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/bench.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/bench.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/benchext.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/benchext.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/blitz.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/bounds.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/bzdebug.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/bzconfig.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/compiler.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/constpointerstack.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/etbase.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/et-forward.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/funcs.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/globeval.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/indexexpr.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/indexmap-forward.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/levicivita.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/limits-hack.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/listinit.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/memblock.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/memblock.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/minmax.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/numinquire.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/numtrait.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/ops.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/prettyprint.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/promote.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/range.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/range.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/ranks.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/reduce.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/shapecheck.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/simdtypes.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tau.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/timer.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinymat2.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinymat2.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinymat2io.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinyvec2.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinyvec2.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tinyvec2io.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tm2fastiter.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tmevaluate.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tv2fastiter.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tvevaluate.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/traversal.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/traversal.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tuning.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tvcross.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/tvecglobs.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/update.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/thirdparty/blitz/v1.0.2/blitz/wrap-climits.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/vecbops.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/vecuops.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/vecwhere.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/vecbfn.cc"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/matbops.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/matuops.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/mathfunc.h"
    "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/promote-old.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/blitz" TYPE FILE FILES "/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz//bzconfig.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/generate/cmake_install.cmake")
  include("/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/meta/cmake_install.cmake")
  include("/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/170-fix-include-correctness/build-verify/thirdparty/blitz/v1.0.2-build/blitz/array/cmake_install.cmake")

endif()

