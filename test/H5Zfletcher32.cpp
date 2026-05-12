#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/H5Zall.hpp>

#include <cstring>
#include <vector>
#include <numeric>

namespace {
std::vector<uint8_t> make_payload(size_t size) {
    std::vector<uint8_t> v(size);
    std::iota(v.begin(), v.end(), static_cast<uint8_t>(0));
    return v;
}
}

TEST_CASE("fletcher32: encode grows output by 4 bytes") {
    const auto input = make_payload(64);
    std::vector<uint8_t> out(input.size() + 4);
    const size_t result = h5::impl::filter::fletcher32(
        out.data(), input.data(), input.size(), 0, 0, nullptr);
    CHECK(result == input.size() + 4);
}

TEST_CASE("fletcher32: encode + decode round-trip") {
    const auto input = make_payload(128);
    const size_t nbytes = input.size();

    std::vector<uint8_t> encoded(nbytes + 4);
    std::vector<uint8_t> decoded(nbytes);

    const size_t enc_size = h5::impl::filter::fletcher32(
        encoded.data(), input.data(), nbytes, 0, 0, nullptr);
    REQUIRE(enc_size == nbytes + 4);

    const size_t dec_size = h5::impl::filter::fletcher32(
        decoded.data(), encoded.data(), enc_size, H5Z_FLAG_REVERSE, 0, nullptr);
    REQUIRE(dec_size == nbytes);
    CHECK(decoded == input);
}

TEST_CASE("fletcher32: detects single-byte corruption") {
    const auto input = make_payload(64);
    const size_t nbytes = input.size();

    std::vector<uint8_t> encoded(nbytes + 4);
    h5::impl::filter::fletcher32(
        encoded.data(), input.data(), nbytes, 0, 0, nullptr);

    // Corrupt a data byte
    encoded[10] ^= 0xFF;

    std::vector<uint8_t> decoded(nbytes);
    const size_t dec_size = h5::impl::filter::fletcher32(
        decoded.data(), encoded.data(), encoded.size(), H5Z_FLAG_REVERSE, 0, nullptr);
    CHECK(dec_size == 0);
}

TEST_CASE("fletcher32: detects checksum corruption") {
    const auto input = make_payload(64);
    const size_t nbytes = input.size();

    std::vector<uint8_t> encoded(nbytes + 4);
    h5::impl::filter::fletcher32(
        encoded.data(), input.data(), nbytes, 0, 0, nullptr);

    // Corrupt checksum bytes
    encoded[nbytes] ^= 0xFF;

    std::vector<uint8_t> decoded(nbytes);
    const size_t dec_size = h5::impl::filter::fletcher32(
        decoded.data(), encoded.data(), encoded.size(), H5Z_FLAG_REVERSE, 0, nullptr);
    CHECK(dec_size == 0);
}

TEST_CASE("fletcher32: odd-size payload round-trips") {
    const auto input = make_payload(65); // odd size
    const size_t nbytes = input.size();

    std::vector<uint8_t> encoded(nbytes + 4);
    std::vector<uint8_t> decoded(nbytes);

    const size_t enc_size = h5::impl::filter::fletcher32(
        encoded.data(), input.data(), nbytes, 0, 0, nullptr);
    REQUIRE(enc_size == nbytes + 4);

    const size_t dec_size = h5::impl::filter::fletcher32(
        decoded.data(), encoded.data(), enc_size, H5Z_FLAG_REVERSE, 0, nullptr);
    REQUIRE(dec_size == nbytes);
    CHECK(decoded == input);
}
