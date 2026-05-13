// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/linalg/blitz.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <blitz/array.h>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] blitz Array round-trip") {
    const char* filename = "test_blitz_io.h5";
    std::filesystem::remove(filename);

    using Matrix = blitz::Array<double, 2>;
    using Vector = blitz::Array<double, 1>;

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

        CHECK(readback.rows() == M.rows());
        CHECK(readback.cols() == M.cols());
        for (int i = 0; i < M.rows(); ++i) {
            for (int j = 0; j < M.cols(); ++j) {
                CHECK(readback(i, j) == doctest::Approx(M(i, j)));
            }
        }
    }

    // BUILD and WRITE vector
    Vector V(4);
    for (blitz::sizeType i = 0; i < V.size(); ++i)
        V(i) = static_cast<double>(i * 2);

    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDWR);
        h5::write(fd, "vector", V);
    }

    // READ BACK vector
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<Vector>(fd, "vector");

        CHECK(readback.size() == V.size());
        for (blitz::sizeType i = 0; i < V.size(); ++i) {
            CHECK(readback(i) == doctest::Approx(V(i)));
        }
    }

    std::filesystem::remove(filename);
}
