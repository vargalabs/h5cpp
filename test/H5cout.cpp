#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5cout.hpp>
#include <sstream>
#include "support/fixture.hpp"

TEST_CASE("operator<< for dxpl_t prints handle") {
    std::ostringstream oss;
    oss << h5::default_dxpl;
    CHECK(oss.str().find("handle") != std::string::npos);
}

TEST_CASE("operator<< for sp_t prints rank and dims") {
    h5::test::file_fixture_t f("test-cout-sp.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{3, 4});
    h5::sp_t sp = h5::get_space(ds);
    std::ostringstream oss;
    oss << sp;
    CHECK(oss.str().find("rank") != std::string::npos);
}

TEST_CASE("operator<< for impl::array prints dims") {
    h5::current_dims_t dims{3, 4, 5};
    std::ostringstream oss;
    oss << dims;
    CHECK(oss.str().find("3") != std::string::npos);
    CHECK(oss.str().find("4") != std::string::npos);
    CHECK(oss.str().find("5") != std::string::npos);
}

TEST_CASE("operator<< for impl::array with inf max_dims") {
    h5::max_dims_t dims{10, H5S_UNLIMITED};
    std::ostringstream oss;
    oss << dims;
    CHECK(oss.str().find("inf") != std::string::npos);
}

TEST_CASE("operator<< for std::vector prints elements") {
    std::vector<int> vec = {1, 2, 3};
    std::ostringstream oss;
    oss << vec;
    CHECK(oss.str().find("1") != std::string::npos);
}
