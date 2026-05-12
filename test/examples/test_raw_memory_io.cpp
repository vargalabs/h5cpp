// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/raw_memory/raw.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <cstddef>
#include <cstring>
#include <filesystem>

TEST_CASE("[example] raw memory pointer round-trip") {
    const char* filename = "test_raw_memory_io.h5";
    std::filesystem::remove(filename);

    double ad[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double* ptr = static_cast<double*>(std::calloc(10, sizeof(double)));

    // WRITE: 1D memory into 2D file space
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write<double>(fd, "dataset", ad, h5::count{1, 10});
    }

    // READ BACK: 2D file space -> 1D mem space with offset
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        h5::read<double>(fd, "dataset", ptr, h5::count{1, 8}, h5::offset{0, 2});

        // ptr[0..7] should contain ad[2..9]
        for (int i = 0; i < 8; ++i) {
            CHECK(ptr[i] == doctest::Approx(ad[i + 2]));
        }
    }

    // READ BACK: different count with memory offset
    {
        std::memset(ptr, 0, 10 * sizeof(double));
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        h5::read<double>(fd, "dataset", ptr + 4, h5::count{1, 3});

        // ptr[4..6] should contain ad[0..2]
        for (int i = 0; i < 3; ++i) {
            CHECK(ptr[4 + i] == doctest::Approx(ad[i]));
        }
    }

    std::free(ptr);
    std::filesystem::remove(filename);
}
