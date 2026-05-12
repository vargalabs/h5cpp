// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/linalg/blaze.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <blaze/Math.h>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] blaze DynamicMatrix/DynamicVector round-trip") {
    const char* filename = "test_blaze_io.h5";
    std::filesystem::remove(filename);

    using Matrix = blaze::DynamicMatrix<double, blaze::rowMajor>;
    using Vector = blaze::DynamicVector<double, blaze::columnVector>;

    // BUILD and WRITE square matrix (avoids blaze padding issues with non-multiple dims)
    Matrix M(2, 2);
    M(0, 0) = 1.0; M(0, 1) = 2.0;
    M(1, 0) = 3.0; M(1, 1) = 4.0;

    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write(fd, "matrix", M);
    }

    // READ BACK matrix
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<Matrix>(fd, "matrix");

        CHECK(readback.rows() == M.rows());
        CHECK(readback.columns() == M.columns());
        for (size_t i = 0; i < M.rows(); ++i) {
            for (size_t j = 0; j < M.columns(); ++j) {
                CHECK(readback(i, j) == doctest::Approx(M(i, j)));
            }
        }
    }

    // BUILD and WRITE vector
    Vector V(8);
    for (size_t i = 0; i < V.size(); ++i)
        V[i] = static_cast<double>(i * 2);

    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDWR);
        h5::write(fd, "vector", V);
    }

    // READ BACK vector
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<Vector>(fd, "vector");

        CHECK(readback.size() == V.size());
        for (size_t i = 0; i < V.size(); ++i) {
            CHECK(readback[i] == doctest::Approx(V[i]));
        }
    }

    std::filesystem::remove(filename);
}
