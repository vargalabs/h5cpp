// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/csv/csv2hdf5.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <filesystem>
#include <vector>
#include <cstring>

// Reuse the CSV example struct definition and HDF5 compound type registration
#include "examples/csv/struct.h"
#include "examples/csv/generated.h"

TEST_CASE("[example] compound struct CSV-like round-trip") {
    const char* filename = "test_csv_io.h5";
    std::filesystem::remove(filename);
    auto set_location = [](char (&dst)[STR_ARRAY_SIZE], const char* src) {
        std::memset(dst, 0, STR_ARRAY_SIZE);
        std::strncpy(dst, src, STR_ARRAY_SIZE - 1);
    };

    // BUILD test data (simulating parsed CSV rows)
    std::vector<input_t> original(3);
    original[0] = input_t{902363382, 0, 39.15920668, -86.52587356, {}};
    set_location(original[0].ReportedLocation, "1ST & FESS");
    original[1] = input_t{902364268, 1500, 39.16144, -86.534848, {}};
    set_location(original[1].ReportedLocation, "2ND & COLLEGE");
    original[2] = input_t{902364412, 2300, 39.14978027, -86.56889006, {}};
    set_location(original[2].ReportedLocation, "BASSWOOD & BLOOMFIELD");

    // WRITE using packet table append (same pattern as csv example)
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<input_t>(fd, "simple approach/dataset.csv",
            h5::current_dims{0}, h5::max_dims{H5S_UNLIMITED}, h5::chunk{3} | h5::gzip{9});
        h5::pt_t pt = ds;

        ds["data set"] = "monroe-county-crash-data2003-to-2015.csv";
        ds["csv parser"] = "https://github.com/ben-strasser/fast-cpp-csv-parser";

        for (const auto& row : original) {
            h5::append(pt, row);
        }
    }

    // READ BACK and verify
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<std::vector<input_t>>(fd, "simple approach/dataset.csv");

        CHECK(readback.size() == original.size());
        for (size_t i = 0; i < original.size(); ++i) {
            CHECK(readback[i].MasterRecordNumber == original[i].MasterRecordNumber);
            CHECK(readback[i].Hour == original[i].Hour);
            CHECK(readback[i].Latitude == doctest::Approx(original[i].Latitude));
            CHECK(readback[i].Longitude == doctest::Approx(original[i].Longitude));
            CHECK(std::strncmp(readback[i].ReportedLocation, original[i].ReportedLocation, STR_ARRAY_SIZE) == 0);
        }

        // Verify attributes
        h5::ds_t ds = h5::open(fd, "simple approach/dataset.csv");
        std::string attr_dataset = h5::aread<std::string>(ds, "data set");
        CHECK(attr_dataset == "monroe-county-crash-data2003-to-2015.csv");

        std::string attr_parser = h5::aread<std::string>(ds, "csv parser");
        CHECK(attr_parser == "https://github.com/ben-strasser/fast-cpp-csv-parser");
    }

    std::filesystem::remove(filename);
}
