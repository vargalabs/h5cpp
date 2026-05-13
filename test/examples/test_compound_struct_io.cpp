// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/compound/struct.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <filesystem>
#include <vector>

#include "examples/compound/struct.h"
#include "examples/compound/generated.h"
#include "examples/compound/utils.hpp"

TEST_CASE("[example] compound struct round-trip") {
    const char* filename = "test_compound_struct_io.h5";
    std::filesystem::remove(filename);

    // BUILD test data
    std::vector<sn::example::Record> original = h5::utils::get_test_data<sn::example::Record>(20);

    // WRITE
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::write(fd, "orm/partial/vector one_shot", original);

        // Also test dataset creation with compound type and custom properties
        h5::create<sn::example::Record>(fd, "/orm/chunked_2D",
            h5::current_dims{4, 5}, h5::chunk{1, 5} | h5::gzip{8});
    }

    // READ BACK
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<std::vector<sn::example::Record>>(fd, "orm/partial/vector one_shot");

        CHECK(readback.size() == original.size());
        for (size_t i = 0; i < original.size(); ++i) {
            CHECK(readback[i].idx == original[i].idx);
        }
    }

    std::filesystem::remove(filename);
}
