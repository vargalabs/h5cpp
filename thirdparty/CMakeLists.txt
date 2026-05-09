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
add_subdirectory(itpp)

# Test framework
add_subdirectory(doctest)
