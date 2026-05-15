/* Copyright (c) 2018-2026 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargalabs.com> */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Pall.hpp>

#ifdef H5_HAVE_ROS3_VFD
#include <H5FDros3.h>

TEST_CASE("fapl_ros3 default ctor — unauthenticated") {
    h5::fapl_ros3 p;
    CHECK(static_cast<hid_t>(p) > 0);

    H5FD_ros3_fapl_t out{};
    REQUIRE(H5Pget_fapl_ros3(static_cast<hid_t>(p), &out) >= 0);
    CHECK(out.authenticate == false);
}

TEST_CASE("fapl_ros3 v1 ctor — authenticated round-trip") {
    h5::fapl_ros3 p{true, "us-east-1", "AKIAIOSFODNN7EXAMPLE", "wJalrXUtnFEMI"};
    CHECK(static_cast<hid_t>(p) > 0);

    H5FD_ros3_fapl_t out{};
    REQUIRE(H5Pget_fapl_ros3(static_cast<hid_t>(p), &out) >= 0);
    CHECK(out.authenticate == true);
    CHECK(std::string(out.aws_region) == "us-east-1");
    CHECK(std::string(out.secret_id)  == "AKIAIOSFODNN7EXAMPLE");
    CHECK(std::string(out.secret_key) == "wJalrXUtnFEMI");
}

TEST_CASE("fapl_ros3 alias ros3 is identical type") {
    h5::ros3 p{false, "", "", ""};
    CHECK(static_cast<hid_t>(p) > 0);
}

TEST_CASE("h5::have_ros3_vfd constexpr is true") {
    CHECK(h5::have_ros3_vfd == true);
}

#if H5_VERSION_GE(1,12,1)
TEST_CASE("fapl_ros3 v2 ctor — session token round-trip") {
    h5::fapl_ros3 p{true, "eu-west-1", "KEY_ID", "SECRET", "SESSION_TOKEN"};
    CHECK(static_cast<hid_t>(p) > 0);

    H5FD_ros3_fapl_t out{};
    REQUIRE(H5Pget_fapl_ros3(static_cast<hid_t>(p), &out) >= 0);
    CHECK(std::string(out.session_token) == "SESSION_TOKEN");
}
#endif

#else  // H5_HAVE_ROS3_VFD not defined

TEST_CASE("h5::have_ros3_vfd constexpr is false when driver absent") {
    CHECK(h5::have_ros3_vfd == false);
}

#endif
