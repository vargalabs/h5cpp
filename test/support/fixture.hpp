/*
 * Copyright (c) 2026 vargaconsulting, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 *
 * Test fixture helpers for h5cpp I/O tests.
 * Issue #114.
 */
#pragma once
#include <hdf5.h>
#include <h5cpp/H5Fcreate.hpp>
#include <string>

namespace h5::test {
    struct file_fixture_t {
        h5::fd_t fd;
        explicit file_fixture_t(const std::string& path = "test-io.h5") {
            fd = h5::create(path, H5F_ACC_TRUNC);
        }
    };
}
