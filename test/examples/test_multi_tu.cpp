// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Compile-time verification of compiler-generated multi-TU headers.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <filesystem>

#include "examples/multi-tu/struct.h"
#include <h5cpp/core>
#include "examples/multi-tu/tu-01.h"
#include <h5cpp/io>
#include "examples/multi-tu/utils.hpp"

TEST_CASE("[example] multi-TU compiler-generated header compiles and runs") {
    const char* filename = "test_multi_tu.h5";
    std::filesystem::remove(filename);

    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);

        // Verify sn::example::Record works with h5::create (uses generated tu-01.h)
        auto ds = h5::create<sn::example::Record>(fd, "/orm/record",
            h5::current_dims{5, 3}, h5::chunk{1, 3} | h5::gzip{8});
        CHECK(ds > 0);

        // Verify sn::typecheck::Record works with h5::create
        auto ds_tc = h5::create<sn::typecheck::Record>(fd, "/orm/typecheck",
            h5::max_dims{H5S_UNLIMITED});
        CHECK(ds_tc > 0);

        // Write and round-trip sn::example::Record using utils helper
        std::vector<sn::example::Record> original = h5::utils::get_test_data<sn::example::Record>(10);
        h5::write(fd, "/data/vector", original);

        auto readback = h5::read<std::vector<sn::example::Record>>(fd, "/data/vector");
        CHECK(readback.size() == original.size());
        for (size_t i = 0; i < original.size(); ++i) {
            CHECK(readback[i].idx == original[i].idx);
        }
    }

    std::filesystem::remove(filename);
}
