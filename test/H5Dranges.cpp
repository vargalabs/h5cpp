/* This file is part of the H5CPP project and is licensed under the MIT License.
 *
 * Copyright © 2025–2026 Varga Consulting, Toronto, ON, Canada 🇨🇦
 *
 * Tests for h5::view<T> — C++20 ranges streaming view (issue #181).
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>

#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5Dranges.hpp>

#include <numeric>
#include <vector>

#include "support/fixture.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static h5::ds_t make_chunked_dataset(h5::fd_t& fd, const char* path,
                                      std::size_t n_elem, std::size_t chunk,
                                      bool compress) {
    if (compress)
        return h5::create<double>(fd, path,
            h5::current_dims_t{n_elem}, h5::max_dims_t{H5S_UNLIMITED},
            h5::chunk{chunk} | h5::gzip{6});
    else
        return h5::create<double>(fd, path,
            h5::current_dims_t{n_elem}, h5::max_dims_t{H5S_UNLIMITED},
            h5::chunk{chunk});
}

static std::vector<double> iota_vec(std::size_t n) {
    std::vector<double> v(n);
    std::iota(v.begin(), v.end(), 0.0);
    return v;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("h5::view iterates uncompressed chunked dataset") {
    h5::test::file_fixture_t f("test-h5dranges-plain.h5");
    constexpr std::size_t N = 40, CHUNK = 10;

    auto ds   = make_chunked_dataset(f.fd, "plain", N, CHUNK, false);
    auto data = iota_vec(N);
    h5::write(ds, data);

    std::vector<double> got;
    got.reserve(N);
    for (double v : h5::view<double>(ds))
        got.push_back(v);

    REQUIRE(got.size() == N);
    for (std::size_t i = 0; i < N; ++i)
        CHECK(got[i] == data[i]);
}

TEST_CASE("h5::view iterates gzip-compressed chunked dataset") {
    h5::test::file_fixture_t f("test-h5dranges-gzip.h5");
    constexpr std::size_t N = 100, CHUNK = 20;

    auto ds   = make_chunked_dataset(f.fd, "gzip", N, CHUNK, true);
    auto data = iota_vec(N);
    h5::write(ds, data);

    std::vector<double> got;
    got.reserve(N);
    for (double v : h5::view<double>(ds))
        got.push_back(v);

    REQUIRE(got.size() == N);
    for (std::size_t i = 0; i < N; ++i)
        CHECK(got[i] == data[i]);
}

TEST_CASE("h5::view works when element count is an exact multiple of chunk size") {
    h5::test::file_fixture_t f("test-h5dranges-exact.h5");
    constexpr std::size_t N = 30, CHUNK = 10;

    auto ds   = make_chunked_dataset(f.fd, "exact", N, CHUNK, false);
    auto data = iota_vec(N);
    h5::write(ds, data);

    std::size_t count = 0;
    for (double v : h5::view<double>(ds)) {
        CHECK(v == data[count]);
        ++count;
    }
    CHECK(count == N);
}

TEST_CASE("h5::view works with int element type") {
    h5::test::file_fixture_t f("test-h5dranges-int.h5");
    constexpr std::size_t N = 50, CHUNK = 7;

    h5::ds_t ds = h5::create<int>(f.fd, "ints",
        h5::current_dims_t{N}, h5::max_dims_t{H5S_UNLIMITED},
        h5::chunk{CHUNK});

    std::vector<int> data(N);
    std::iota(data.begin(), data.end(), 0);
    h5::write(ds, data);

    std::vector<int> got;
    got.reserve(N);
    for (int v : h5::view<int>(ds))
        got.push_back(v);

    REQUIRE(got.size() == N);
    for (std::size_t i = 0; i < N; ++i)
        CHECK(got[i] == data[i]);
}

TEST_CASE("h5::view on rank-2 dataset throws") {
    h5::test::file_fixture_t f("test-h5dranges-rank2.h5");
    h5::ds_t ds = h5::create<double>(f.fd, "mat",
        h5::current_dims_t{4, 4}, h5::max_dims_t{H5S_UNLIMITED, H5S_UNLIMITED},
        h5::chunk{2, 2});
    // nodiscard result intentionally discarded — we only care about the throw
    CHECK_THROWS_AS((void)h5::view<double>(ds), std::runtime_error);
}
