/* Copyright (c) 2018-2026 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargalabs.com> */

// Network-dependent ROS3 smoke test.
// Requires HDF5 built with ROS3 support and outbound S3 connectivity.
// Enable with: cmake -DH5CPP_NETWORK_TESTS=ON
//
// Test file: s3://rhdf5-public/h5ex_t_array.h5  (eu-central-1, public read)
//   /DS1  —  H5T_ARRAY{[3][5] int64},  4 elements

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Pall.hpp>
#include <hdf5.h>

#ifdef H5_HAVE_ROS3_VFD

static constexpr const char* BUCKET  = "s3://rhdf5-public/h5ex_t_array.h5";
static constexpr const char* REGION  = "eu-central-1";
static constexpr const char* DATASET = "/DS1";

TEST_CASE("ros3 unauthenticated open — rhdf5-public bucket") {
    h5::ros3 fapl{false, REGION, "", ""};
    h5::fd_t fd = h5::open(BUCKET, H5F_ACC_RDONLY, fapl);
    CHECK(H5Iis_valid(static_cast<hid_t>(fd)));
}

TEST_CASE("ros3 dataset existence — /DS1 present") {
    h5::ros3 fapl{false, REGION, "", ""};
    h5::fd_t fd = h5::open(BUCKET, H5F_ACC_RDONLY, fapl);
    htri_t exists = H5Lexists(static_cast<hid_t>(fd), DATASET, H5P_DEFAULT);
    CHECK(exists > 0);
}

TEST_CASE("ros3 dataset dimensions — DS1 has 4 elements") {
    h5::ros3 fapl{false, REGION, "", ""};
    h5::fd_t fd = h5::open(BUCKET, H5F_ACC_RDONLY, fapl);

    hid_t ds   = H5Dopen(static_cast<hid_t>(fd), DATASET, H5P_DEFAULT);
    hid_t sp   = H5Dget_space(ds);
    hsize_t dims[1] = {};
    int ndims = H5Sget_simple_extent_dims(sp, dims, nullptr);

    CHECK(ndims == 1);
    CHECK(dims[0] == 4);

    H5Sclose(sp);
    H5Dclose(ds);
}

#else

TEST_CASE("ros3 network tests skipped — H5_HAVE_ROS3_VFD not defined") {
    MESSAGE("ROS3 VFD not compiled into this HDF5 build; skipping.");
    CHECK(true);
}

#endif
