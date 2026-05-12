// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/utf/utf.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <armadillo>
#include <h5cpp/all>
#include <filesystem>
#include <string>
#include <vector>

TEST_CASE("[example] UTF-8 dataset path round-trip") {
    const char* filename = "test_utf_io.h5";
    std::filesystem::remove(filename);

    std::vector<std::string> paths = {
        "hello world",
        "مرحبا بالعالم",
        "你好，世界",
        "こんにちは世界",
        "안녕 세상"
    };

    // WRITE to UTF-8 named datasets
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC,
            h5::default_fcpl, h5::libver_bounds({H5F_LIBVER_V18, H5F_LIBVER_V18}));
        for (const auto& path : paths) {
            arma::mat M = arma::ones(3, 3);
            h5::write(fd, path, M);
        }
    }

    // READ BACK and verify
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        for (const auto& path : paths) {
            auto readback = h5::read<arma::mat>(fd, path);
            CHECK(readback.n_rows == 3);
            CHECK(readback.n_cols == 3);
            CHECK(readback(0, 0) == doctest::Approx(1.0));
        }
    }

    std::filesystem::remove(filename);
}
