#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Sall.hpp>

TEST_CASE("get_current_dims with count only") {
    h5::count_t count{10, 20};
    auto dims = h5::impl::get_current_dims(count);
    CHECK(dims.rank == 2);
    CHECK(dims[0] == 10);
    CHECK(dims[1] == 20);
}

TEST_CASE("get_current_dims with count and stride") {
    h5::count_t count{5, 3};
    h5::stride_t stride{2, 4};
    auto dims = h5::impl::get_current_dims(count, stride);
    CHECK(dims.rank == 2);
    CHECK(dims[0] == 10);
    CHECK(dims[1] == 12);
}

TEST_CASE("get_current_dims with count, stride, and block") {
    h5::count_t count{5, 3};
    h5::stride_t stride{3, 5};
    h5::block_t block{2, 3};
    auto dims = h5::impl::get_current_dims(count, stride, block);
    CHECK(dims.rank == 2);
    CHECK(dims[0] == 10);
    CHECK(dims[1] == 9);
}

TEST_CASE("get_current_dims with count, stride, block, and offset") {
    h5::count_t count{5, 3};
    h5::stride_t stride{3, 5};
    h5::block_t block{2, 3};
    h5::offset_t offset{1, 2};
    auto dims = h5::impl::get_current_dims(count, stride, block, offset);
    CHECK(dims.rank == 2);
    CHECK(dims[0] == 11);
    CHECK(dims[1] == 11);
}

TEST_CASE("get_current_dims with count and offset") {
    h5::count_t count{10};
    h5::offset_t offset{5};
    auto dims = h5::impl::get_current_dims(count, offset);
    CHECK(dims.rank == 1);
    CHECK(dims[0] == 15);
}
