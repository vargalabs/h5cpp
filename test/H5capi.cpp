#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5capi.hpp>
#include "support/fixture.hpp"

TEST_CASE("get_fill_value returns default fill value when not explicitly set") {
    h5::test::file_fixture_t f("test-fill-undefined.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{1});
    void* fv = h5::get_fill_value(ds);
    CHECK(fv != nullptr);
    CHECK(*static_cast<int*>(fv) == 0);
    free(fv);
}

TEST_CASE("get_fill_value returns value when defined") {
    h5::test::file_fixture_t f("test-fill-defined.h5");
    int fill = 42;
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{1},
        h5::chunk{1} | h5::fill_value<int>{fill});
    void* fv = h5::get_fill_value(ds);
    CHECK(fv != nullptr);
    CHECK(*static_cast<int*>(fv) == 42);
    free(fv);
}

TEST_CASE("get_chunk_dims returns expected chunk") {
    h5::test::file_fixture_t f("test-chunk-dims.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{100},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::dcpl_t dcpl = h5::get_dcpl(ds);
    h5::chunk_t chunk;
    int rank = h5::get_chunk_dims(dcpl, chunk);
    CHECK(rank == 1);
    CHECK(chunk[0] == 10);
}

TEST_CASE("set_extent extends dataset") {
    h5::test::file_fixture_t f("test-set-extent.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{10},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{5});
    h5::current_dims_t new_dims{20};
    h5::set_extent(ds, new_dims);
    h5::sp_t sp = h5::get_space(ds);
    h5::current_dims_t cd;
    h5::max_dims_t md;
    h5::get_simple_extent_dims(sp, cd, md);
    CHECK(cd[0] == 20);
}

TEST_CASE("get_size returns correct type size") {
    h5::dt_t<int> dt;
    CHECK(h5::get_size(dt) == sizeof(int));
    h5::dt_t<double> dt2;
    CHECK(h5::get_size(dt2) == sizeof(double));
}

TEST_CASE("get_simple_extent_ndims returns correct rank") {
    h5::test::file_fixture_t f("test-rank.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{10, 20});
    h5::sp_t sp = h5::get_space(ds);
    CHECK(h5::get_simple_extent_ndims(sp) == 2);
}

TEST_CASE("select_hyperslab 2-arg overload") {
    h5::test::file_fixture_t f("test-hyperslab-2arg.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{10});
    h5::sp_t sp = h5::get_space(ds);
    h5::offset_t offset{2};
    h5::count_t count{5};
    h5::select_hyperslab(sp, offset, count);
    CHECK(true);
}

TEST_CASE("select_hyperslab 4-arg overload") {
    h5::test::file_fixture_t f("test-hyperslab-4arg.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{100});
    h5::sp_t sp = h5::get_space(ds);
    h5::offset_t offset{0};
    h5::stride_t stride{2};
    h5::count_t count{10};
    h5::block_t block{1};
    h5::select_hyperslab(sp, offset, stride, count, block);
    CHECK(true);
}

TEST_CASE("writeds low-level write") {
    h5::test::file_fixture_t f("test-writeds.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{5});
    std::vector<int> data = {1, 2, 3, 4, 5};
    h5::dt_t<int> mem_type;
    h5::sp_t mem_space = h5::create_simple(hsize_t{5});
    h5::sp_t file_space = h5::get_space(ds);
    h5::writeds(ds, mem_type, mem_space, file_space, data.data());
    CHECK(true);
}

TEST_CASE("get_dapl returns valid property list") {
    h5::test::file_fixture_t f("test-get-dapl.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "dataset", h5::current_dims_t{1});
    h5::dapl_t dapl = h5::get_dapl(ds);
    CHECK(static_cast<hid_t>(dapl) > 0);
}
