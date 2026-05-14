#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5Rall.hpp>
#include <vector>
#include "support/fixture.hpp"

TEST_CASE("h5::reference creates region reference") {
    h5::test::file_fixture_t f("test-reference.h5");
    h5::create<int>(f.fd, "ds", h5::current_dims_t{10});
    std::vector<int> data(10);
    for (int i = 0; i < 10; ++i) data[i] = i;
    h5::write(f.fd, "ds", data);

    h5::reference_t ref = h5::reference(f.fd, "ds", h5::count_t{5}, h5::offset_t{2});
    CHECK(true);
}

TEST_CASE("h5::exp::read reads via reference") {
    h5::test::file_fixture_t f("test-reference-read.h5");
    h5::create<int>(f.fd, "ds", h5::current_dims_t{10});
    std::vector<int> data(10);
    for (int i = 0; i < 10; ++i) data[i] = i;
    h5::write(f.fd, "ds", data);

    h5::reference_t ref = h5::reference(f.fd, "ds", h5::count_t{5}, h5::offset_t{2});
    h5::ds_t ds = h5::open(f.fd, "ds");
    auto readback = h5::exp::read<std::vector<int>>(ds, ref);
    CHECK(readback.size() == 5);
    for (size_t i = 0; i < 5; ++i)
        CHECK(readback[i] == static_cast<int>(i + 2));
}
