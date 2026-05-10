#include <hdf5.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Fcreate.hpp>
#include <h5cpp/H5Fopen.hpp>
#include <h5cpp/H5Dcreate.hpp>
#include <h5cpp/H5Dopen.hpp>
#include <h5cpp/H5Acreate.hpp>
#include <h5cpp/H5Aopen.hpp>
#include <h5cpp/H5Awrite.hpp>
#include <h5cpp/H5Adelete.hpp>
#include "support/fixture.hpp"

TEST_CASE("h5::adelete removes an existing attribute") {
    h5::test::file_fixture_t f("test-adelete.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{1});

    ds["attr_to_delete"] = 42;
    CHECK(H5Aexists(static_cast<hid_t>(ds), "attr_to_delete") > 0);

    h5::adelete(ds, "attr_to_delete");
    CHECK(H5Aexists(static_cast<hid_t>(ds), "attr_to_delete") == 0);
}

TEST_CASE("h5::adelete throws on non-existent attribute") {
    h5::test::file_fixture_t f("test-adelete-missing.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{1});

    bool threw = false;
    try {
        h5::adelete(ds, "missing_attr");
    } catch (const h5::error::io::attribute::delete_&) {
        threw = true;
    }
    CHECK(threw);
}
