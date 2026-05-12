// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/optimized/optimized.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <armadillo>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] optimized offset write/read round-trip") {
    const char* filename = "test_optimized_io.h5";
    std::filesystem::remove(filename);

    arma::mat M(10, 1);
    M.zeros();

    // CREATE dataset with chunked/gzipped properties
    h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
    h5::ds_t ds = h5::create<short>(fd, "dataset",
        h5::current_dims{10, 5},
        h5::max_dims{10, H5S_UNLIMITED},
        h5::chunk{10, 1} | h5::gzip{9} | h5::fill_value<short>{0});

    // WRITE with offsets
    for (hsize_t i = 0; i < 5; ++i) {
        M.zeros();
        M(0, 0) = static_cast<double>(i);
        h5::write(ds, M, h5::count{10, 1}, h5::offset{0, i});
    }

    // READ BACK and verify each offset
    for (hsize_t i = 0; i < 5; ++i) {
        auto R = h5::read<arma::mat>(ds, h5::count{10, 1}, h5::offset{0, i});
        CHECK(R(0, 0) == doctest::Approx(static_cast<double>(i)));
    }

    std::filesystem::remove(filename);
}
