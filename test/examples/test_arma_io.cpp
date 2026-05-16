// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/linalg/arma.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <armadillo>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] armadillo mat/vec round-trip") {
    const char* filename = "test_arma_io.h5";
    std::filesystem::remove(filename);

    // BUILD and WRITE non-square matrix — shape (3 rows x 5 cols) catches row/col transposition bugs
    arma::mat M = arma::linspace<arma::rowvec>(1.0, 15.0, 15);
    M.reshape(3, 5);

    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write(fd, "matrix", M);
    }

    // READ BACK matrix
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<arma::mat>(fd, "matrix");

        CHECK(readback.n_rows == M.n_rows);
        CHECK(readback.n_cols == M.n_cols);
        for (arma::uword i = 0; i < M.n_elem; ++i) {
            CHECK(readback(i) == doctest::Approx(M(i)));
        }
    }

    // BUILD and WRITE vector
    arma::vec V(8);
    for (arma::uword i = 0; i < V.n_elem; ++i)
        V(i) = static_cast<double>(i * 2);

    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDWR);
        h5::write(fd, "vector", V);
    }

    // READ BACK vector
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<arma::vec>(fd, "vector");

        CHECK(readback.n_elem == V.n_elem);
        for (arma::uword i = 0; i < V.n_elem; ++i) {
            CHECK(readback(i) == doctest::Approx(V(i)));
        }
    }

    std::filesystem::remove(filename);
}
