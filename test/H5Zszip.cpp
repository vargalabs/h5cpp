#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/H5Zall.hpp>

#include <cstring>
#include <vector>
#include <numeric>

namespace {
// 16-bit integer data — typical sensor payload for szip.
std::vector<uint16_t> make_int16_payload(size_t count) {
    std::vector<uint16_t> v(count);
    for (size_t i = 0; i < count; ++i)
        v[i] = static_cast<uint16_t>(i % 4096); // 12-bit range
    return v;
}
}

TEST_CASE("szip: library availability check") {
#if defined(H5CPP_HAS_SZIP)
    CHECK(true);
    MESSAGE("SZIP is active (vendored v2.2.0)");
#else
    WARN("SZIP not configured; filter falls back to passthrough");
#endif
}

TEST_CASE("szip: encode + decode round-trip") {
    // HDF5 params layout: options_mask, bits_per_pixel, pixels_per_block, pixels_per_scanline
    const size_t count = 512;
    const auto input = make_int16_payload(count);
    const size_t nbytes = input.size() * sizeof(uint16_t);

    // Worst-case output buffer: szip cannot expand beyond ~1.5x + overhead.
    std::vector<uint8_t> compressed(nbytes + nbytes / 2 + 128);
    std::vector<uint8_t> decompressed(nbytes);

    // SZ_NN_OPTION_MASK | SZ_MSB_OPTION_MASK, 16 bpp, 32 ppb, count pixels
    const unsigned params[] = {
        32u | 16u,  // SZ_NN_OPTION_MASK=32 | SZ_MSB_OPTION_MASK=16
        16u,        // bits_per_pixel = 16
        32u,        // pixels_per_block = 32
        static_cast<unsigned>(count)
    };

    const size_t comp_size = h5::impl::filter::szip(
        compressed.data(), input.data(), nbytes, 0, 4, params);

#if defined(H5CPP_HAS_SZIP)
    REQUIRE(comp_size > 0);

    const size_t decomp_size = h5::impl::filter::szip(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 4, params);

    REQUIRE(decomp_size == nbytes);
    CHECK(std::memcmp(decompressed.data(), input.data(), nbytes) == 0);
#else
    // Passthrough: encode returns original size, data is unchanged.
    CHECK(comp_size == nbytes);
#endif
}

TEST_CASE("szip: compressed output is smaller than input for compressible data") {
    const size_t count = 1024;
    const auto input = make_int16_payload(count);
    const size_t nbytes = input.size() * sizeof(uint16_t);

    std::vector<uint8_t> compressed(nbytes + nbytes / 2 + 128);

    const unsigned params[] = {32u | 16u, 16u, 32u, static_cast<unsigned>(count)};
    const size_t comp_size = h5::impl::filter::szip(
        compressed.data(), input.data(), nbytes, 0, 4, params);

    REQUIRE(comp_size > 0);
#if defined(H5CPP_HAS_SZIP)
    CHECK(comp_size < nbytes);
#endif
}

TEST_CASE("szip: returns 0 on insufficient params") {
    // n < 4 must return 0 rather than crash.
    std::vector<uint8_t> buf(64);
    const unsigned params[] = {32u, 16u}; // only n=2
    const size_t result = h5::impl::filter::szip(
        buf.data(), buf.data(), buf.size(), 0, 2, params);
#if defined(H5CPP_HAS_SZIP)
    CHECK(result == 0);
#else
    CHECK(result == buf.size()); // passthrough
#endif
}
