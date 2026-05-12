// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/transform/transform.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <armadillo>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] data transform round-trip") {
    const char* filename = "test_transform_io.h5";
    std::filesystem::remove(filename);

    arma::mat M(5, 5);
    M.ones();

    // WRITE with data transform: 2*x + 5
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write(fd, "transform", M, h5::data_transform{"2*x+5"});
    }

    // READ BACK with inverse transform: (x - 5) / 2
    {
        auto readback = h5::read<arma::mat>(filename, "transform", h5::data_transform{"(x-5)/2"});

        CHECK(readback.n_rows == M.n_rows);
        CHECK(readback.n_cols == M.n_cols);
        for (arma::uword i = 0; i < M.n_elem; ++i) {
            CHECK(readback(i) == doctest::Approx(M(i)));
        }
    }

    std::filesystem::remove(filename);
}
