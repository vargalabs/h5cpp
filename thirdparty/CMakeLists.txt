# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.

# Linear algebra third-party dependencies
add_subdirectory(eigen3)
add_subdirectory(blaze)
add_subdirectory(ublas)
add_subdirectory(armadillo)
add_subdirectory(dlib)
add_subdirectory(blitz)
# NOTE: itpp disabled — upstream v4.3.1 uses the 'register' keyword extensively,
# which was removed in C++17. Building with GCC 14 / -std=c++17 fails.
# This is an ITPP upstream compatibility issue, not an h5cpp bug.
# add_subdirectory(itpp)

# Test framework
add_subdirectory(doctest)
add_subdirectory(libdeflate)
