// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/linalg/dlib.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <dlib/matrix.h>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] dlib matrix round-trip") {
    const char* filename = "test_dlib_io.h5";
    std::filesystem::remove(filename);

    using Matrix = dlib::matrix<double>;

    // BUILD and WRITE matrix
    Matrix M(2, 3);
    M(0, 0) = 1.0; M(0, 1) = 2.0; M(0, 2) = 3.0;
    M(1, 0) = 4.0; M(1, 1) = 5.0; M(1, 2) = 6.0;

    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write(fd, "matrix", M);
    }

    // READ BACK matrix
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<Matrix>(fd, "matrix");

        CHECK(readback.nr() == M.nr());
        CHECK(readback.nc() == M.nc());
        for (long i = 0; i < M.nr(); ++i) {
            for (long j = 0; j < M.nc(); ++j) {
                CHECK(readback(i, j) == doctest::Approx(M(i, j)));
            }
        }
    }

    std::filesystem::remove(filename);
}
