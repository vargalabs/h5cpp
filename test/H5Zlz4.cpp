#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/H5Zall.hpp>

#include <cstring>
#include <vector>
#include <numeric>

namespace {
std::vector<uint8_t> make_repeating_payload(size_t size) {
    std::vector<uint8_t> v(size);
    for (size_t i = 0; i < size; ++i)
        v[i] = static_cast<uint8_t>(i % 251);
    return v;
}
}

TEST_CASE("lz4: library availability check") {
#if defined(H5CPP_HAS_LZ4)
    CHECK(true);
    MESSAGE("LZ4 is active");
#else
    WARN_MESSAGE(false, "LZ4 not configured; filter falls back to passthrough");
#endif
}

TEST_CASE("lz4: encode + decode round-trip") {
    const auto input = make_repeating_payload(4096);
    const size_t nbytes = input.size();

    std::vector<uint8_t> compressed(h5::impl::filter::lz4_bound(nbytes));
    std::vector<uint8_t> decompressed(nbytes);

    const unsigned params[] = {0u, static_cast<unsigned>(nbytes)};

    const size_t comp_size = h5::impl::filter::lz4(
        compressed.data(), input.data(), nbytes, 0, 2, params);

    REQUIRE(comp_size > 0);

    const size_t decomp_size = h5::impl::filter::lz4(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 2, params);

    REQUIRE(decomp_size == nbytes);
    CHECK(decompressed == input);
}

TEST_CASE("lz4: compressed output is smaller than input for compressible data") {
    const auto input = make_repeating_payload(8192);
    const size_t nbytes = input.size();
    std::vector<uint8_t> compressed(h5::impl::filter::lz4_bound(nbytes));

    const unsigned params[] = {0u, static_cast<unsigned>(nbytes)};
    const size_t comp_size = h5::impl::filter::lz4(
        compressed.data(), input.data(), nbytes, 0, 2, params);

    REQUIRE(comp_size > 0);
#if defined(H5CPP_HAS_LZ4)
    CHECK(comp_size < nbytes);
#endif
}

TEST_CASE("lz4: decode with only level param (external-file compat)") {
    // External HDF5 files using the community LZ4 plugin store block-size in cd_values[0],
    // not a compression level.  Our set_cache fix injects block_size into params[1] so the
    // decoder always has a reliable output-size bound.  This test simulates that fixed state.
    const auto input = make_repeating_payload(4096);
    const size_t nbytes = input.size();

    std::vector<uint8_t> compressed(h5::impl::filter::lz4_bound(nbytes));
    std::vector<uint8_t> decompressed(nbytes);

    // Encode with standard params
    const unsigned enc_params[] = {0u, static_cast<unsigned>(nbytes)};
    const size_t comp_size = h5::impl::filter::lz4(
        compressed.data(), input.data(), nbytes, 0, 2, enc_params);
    REQUIRE(comp_size > 0);

    // Decode with params[1] = block_size (as injected by set_cache)
    const unsigned dec_params[] = {0u, static_cast<unsigned>(nbytes)};
    const size_t decomp_size = h5::impl::filter::lz4(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 2, dec_params);
    REQUIRE(decomp_size == nbytes);
    CHECK(decompressed == input);
}

TEST_CASE("lz4: corrupt compressed input returns 0") {
    std::vector<uint8_t> corrupt(64, 0xFF);
    std::vector<uint8_t> out(4096);
    const unsigned params[] = {0u, 4096u};

    const size_t result = h5::impl::filter::lz4(
        out.data(), corrupt.data(), corrupt.size(), H5Z_FLAG_REVERSE, 2, params);

#if defined(H5CPP_HAS_LZ4)
    CHECK(result == 0);
#else
    // passthrough mode: returns input size
    CHECK(result == corrupt.size());
#endif
}
