#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/H5Zall.hpp>

#include <cstring>
#include <vector>
#include <numeric>

namespace {
// Generate a buffer of count doubles with predictable values.
std::vector<double> make_doubles(size_t count) {
    std::vector<double> v(count);
    for (size_t i = 0; i < count; ++i)
        v[i] = static_cast<double>(i) * 1.23456789;
    return v;
}
}

TEST_CASE("shuffle: identity for type_size == 1") {
    const std::vector<uint8_t> input = {1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<uint8_t> out(input.size());
    const unsigned params[] = {1};
    const size_t result = h5::impl::filter::shuffle(
        out.data(), input.data(), input.size(), 0, 1, params);
    REQUIRE(result == input.size());
    CHECK(out == input);
}

TEST_CASE("shuffle: round-trip for 4-byte elements") {
    const size_t count = 64;
    std::vector<uint32_t> input(count);
    std::iota(input.begin(), input.end(), 0u);

    const size_t nbytes = count * sizeof(uint32_t);
    std::vector<uint8_t> shuffled(nbytes);
    std::vector<uint8_t> restored(nbytes);

    const unsigned params[] = {sizeof(uint32_t)};

    const size_t fwd = h5::impl::filter::shuffle(
        shuffled.data(), input.data(), nbytes, 0, 1, params);
    REQUIRE(fwd == nbytes);

    const size_t rev = h5::impl::filter::shuffle(
        restored.data(), shuffled.data(), nbytes, H5Z_FLAG_REVERSE, 1, params);
    REQUIRE(rev == nbytes);
    CHECK(std::memcmp(restored.data(), input.data(), nbytes) == 0);
}

TEST_CASE("shuffle: round-trip for 8-byte (double) elements") {
    const auto input = make_doubles(128);
    const size_t nbytes = input.size() * sizeof(double);

    std::vector<uint8_t> shuffled(nbytes);
    std::vector<double>  restored(input.size());

    const unsigned params[] = {sizeof(double)};

    const size_t fwd = h5::impl::filter::shuffle(
        shuffled.data(), input.data(), nbytes, 0, 1, params);
    REQUIRE(fwd == nbytes);

    CHECK(std::memcmp(shuffled.data(), input.data(), nbytes) != 0);

    const size_t rev = h5::impl::filter::shuffle(
        restored.data(), shuffled.data(), nbytes, H5Z_FLAG_REVERSE, 1, params);
    REQUIRE(rev == nbytes);
    CHECK(restored == input);
}

TEST_CASE("shuffle + deflate: chain encodes and decodes correctly") {
    const auto orig = make_doubles(256);
    const size_t nbytes = orig.size() * sizeof(double);

    std::vector<uint8_t> shuffled(nbytes);
    std::vector<uint8_t> compressed(h5::impl::filter::deflate_bound(nbytes));
    std::vector<uint8_t> decompressed(nbytes);
    std::vector<double>  restored(orig.size());

    const unsigned shuf_params[] = {sizeof(double)};
    const unsigned def_params[]  = {6u, static_cast<unsigned>(nbytes)};

    // Encode: shuffle → deflate
    REQUIRE(h5::impl::filter::shuffle(
        shuffled.data(), orig.data(), nbytes, 0, 1, shuf_params) == nbytes);

    const size_t comp_size = h5::impl::filter::deflate(
        compressed.data(), shuffled.data(), nbytes, 0, 2, def_params);
    REQUIRE(comp_size > 0);
    CHECK(comp_size < nbytes); // doubles with shuffle should compress well

    // Decode: inflate → unshuffle
    const size_t decomp_size = h5::impl::filter::deflate(
        decompressed.data(), compressed.data(), comp_size, H5Z_FLAG_REVERSE, 2, def_params);
    REQUIRE(decomp_size == nbytes);

    REQUIRE(h5::impl::filter::shuffle(
        restored.data(), decompressed.data(), nbytes, H5Z_FLAG_REVERSE, 1, shuf_params) == nbytes);

    CHECK(restored == orig);
}
