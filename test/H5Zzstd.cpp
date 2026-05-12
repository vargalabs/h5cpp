#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/H5Zall.hpp>

#include <cstring>
#include <vector>
#include <numeric>

namespace {
std::vector<double> make_double_sequence(size_t count) {
    std::vector<double> v(count);
    for (size_t i = 0; i < count; ++i)
        v[i] = static_cast<double>(i) * 2.718281828;
    return v;
}
}

TEST_CASE("zstd: library availability check") {
#if defined(H5CPP_HAS_ZSTD)
    CHECK(true);
    MESSAGE("Zstd is active");
#else
    WARN("Zstd not configured (install libzstd-dev); filter falls back to passthrough");
#endif
}

TEST_CASE("zstd: encode + decode round-trip") {
    const auto input = make_double_sequence(512);
    const size_t nbytes = input.size() * sizeof(double);

    std::vector<uint8_t> compressed(h5::impl::filter::zstd_bound(nbytes));
    std::vector<uint8_t> decompressed(nbytes);

    const unsigned params[] = {3u, static_cast<unsigned>(nbytes)};

    const size_t comp_size = h5::impl::filter::zstd(
        compressed.data(), input.data(), nbytes, 0, 2, params);

    REQUIRE(comp_size > 0);

    const size_t decomp_size = h5::impl::filter::zstd(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 2, params);

    REQUIRE(decomp_size == nbytes);
    CHECK(std::memcmp(decompressed.data(), input.data(), nbytes) == 0);
}

TEST_CASE("zstd: compressed output is smaller than input for compressible data") {
    const auto input = make_double_sequence(1024);
    const size_t nbytes = input.size() * sizeof(double);
    std::vector<uint8_t> compressed(h5::impl::filter::zstd_bound(nbytes));

    const unsigned params[] = {3u, static_cast<unsigned>(nbytes)};
    const size_t comp_size = h5::impl::filter::zstd(
        compressed.data(), input.data(), nbytes, 0, 2, params);

    REQUIRE(comp_size > 0);
#if defined(H5CPP_HAS_ZSTD)
    CHECK(comp_size < nbytes);
#endif
}

TEST_CASE("zstd: decode with only level param (external-file compat)") {
    // External HDF5 Zstd plugin stores only compression level in cd_values[0];
    // no output size is stored.  set_cache injects block_size into params[1].
    const auto input = make_double_sequence(512);
    const size_t nbytes = input.size() * sizeof(double);

    std::vector<uint8_t> compressed(h5::impl::filter::zstd_bound(nbytes));
    std::vector<uint8_t> decompressed(nbytes);

    const unsigned enc_params[] = {3u, static_cast<unsigned>(nbytes)};
    const size_t comp_size = h5::impl::filter::zstd(
        compressed.data(), input.data(), nbytes, 0, 2, enc_params);
    REQUIRE(comp_size > 0);

    // Decode with params[1] = block_size (as injected by set_cache)
    const unsigned dec_params[] = {3u, static_cast<unsigned>(nbytes)};
    const size_t decomp_size = h5::impl::filter::zstd(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 2, dec_params);
    REQUIRE(decomp_size == nbytes);
    CHECK(std::memcmp(decompressed.data(), input.data(), nbytes) == 0);
}

TEST_CASE("zstd: corrupt compressed input returns 0") {
    std::vector<uint8_t> corrupt(64, 0xFF);
    std::vector<uint8_t> out(8192);
    const unsigned params[] = {3u, 8192u};

    const size_t result = h5::impl::filter::zstd(
        out.data(), corrupt.data(), corrupt.size(), H5Z_FLAG_REVERSE, 2, params);

#if defined(H5CPP_HAS_ZSTD)
    CHECK(result == 0);
#else
    CHECK(result == corrupt.size());
#endif
}
