#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5Zpipeline.hpp>
#include <h5cpp/H5Zpipeline_basic.hpp>
#include <vector>
#include "support/fixture.hpp"

TEST_CASE("basic_pipeline_t set_cache with chunked gzip dataset") {
    h5::test::file_fixture_t f("test-pipeline-cache.h5");
    h5::ds_t ds = h5::create<double>(f.fd, "ds", h5::current_dims_t{100},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10} | h5::gzip{6});
    h5::dcpl_t dcpl = h5::get_dcpl(ds);
    h5::impl::basic_pipeline_t pipeline;
    pipeline.set_cache(dcpl, sizeof(double));
    CHECK(pipeline.rank == 1);
    CHECK(pipeline.n == 10);
}

TEST_CASE("basic_pipeline_t write/read round-trip with gzip filter") {
    h5::test::file_fixture_t f("test-pipeline-roundtrip.h5");
    h5::ds_t ds = h5::create<double>(f.fd, "ds", h5::current_dims_t{100},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10} | h5::gzip{6});
    h5::dcpl_t dcpl = h5::get_dcpl(ds);
    h5::impl::basic_pipeline_t pipeline;
    pipeline.set_cache(dcpl, sizeof(double));

    std::vector<double> data(100);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = static_cast<double>(i);

    h5::offset_t offset{0};
    h5::stride_t stride{1};
    h5::block_t block{1};
    h5::count_t count{100};

    pipeline.write(ds, offset, stride, block, count, h5::default_dxpl, data.data());

    std::vector<double> readback(100);
    pipeline.read(ds, offset, stride, block, count, h5::default_dxpl, readback.data());

    for (size_t i = 0; i < data.size(); ++i)
        CHECK(readback[i] == data[i]);
}

TEST_CASE("basic_pipeline_t write/read with no filters") {
    h5::test::file_fixture_t f("test-pipeline-nofilter.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{50},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::dcpl_t dcpl = h5::get_dcpl(ds);
    h5::impl::basic_pipeline_t pipeline;
    pipeline.set_cache(dcpl, sizeof(int));

    std::vector<int> data(50);
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = static_cast<int>(i * i);

    h5::offset_t offset{0};
    h5::stride_t stride{1};
    h5::block_t block{1};
    h5::count_t count{50};

    pipeline.write(ds, offset, stride, block, count, h5::default_dxpl, data.data());

    std::vector<int> readback(50);
    pipeline.read(ds, offset, stride, block, count, h5::default_dxpl, readback.data());

    for (size_t i = 0; i < data.size(); ++i)
        CHECK(readback[i] == data[i]);
}

TEST_CASE("high_throughput flag is constructible") {
    auto ht = h5::high_throughput;
    CHECK(true);
}
