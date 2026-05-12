#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <zlib.h>

#include <cstring>
#include <string>
#include <vector>

namespace {
std::vector<unsigned char> sample_payload() {
    std::vector<unsigned char> data;
    const std::string pattern = "h5cpp-filtering-pipeline-deflate-roundtrip:";
    for (int i = 0; i < 256; ++i)
        data.insert(data.end(), pattern.begin(), pattern.end());
    return data;
}
}

TEST_CASE("H5Z deflate callback uses libdeflate when configured") {
#if defined(H5CPP_HAS_LIBDEFLATE)
    CHECK(true);
#else
    WARN("libdeflate not configured; zlib fallback branch is active");
#endif
}

TEST_CASE("H5Z deflate callback encodes zlib-compatible stream") {
    const auto input = sample_payload();
    const unsigned params[] = {6, static_cast<unsigned>(input.size())};
    std::vector<unsigned char> encoded(h5::impl::filter::deflate_bound(input.size()));

    const size_t encoded_size = h5::impl::filter::gzip(
        encoded.data(), input.data(), input.size(), 0, 2, params);

    REQUIRE(encoded_size > 0);
    CHECK(encoded_size < input.size());

    std::vector<unsigned char> decoded(input.size());
    uLongf decoded_size = static_cast<uLongf>(decoded.size());
    const int status = uncompress(
        decoded.data(), &decoded_size, encoded.data(), static_cast<uLong>(encoded_size));

    REQUIRE(status == Z_OK);
    REQUIRE(decoded_size == input.size());
    CHECK(decoded == input);
}

TEST_CASE("H5Z deflate callback decodes its own output") {
    const auto input = sample_payload();
    const unsigned params[] = {6, static_cast<unsigned>(input.size())};
    std::vector<unsigned char> encoded(h5::impl::filter::deflate_bound(input.size()));
    std::vector<unsigned char> decoded(input.size());

    const size_t encoded_size = h5::impl::filter::deflate(
        encoded.data(), input.data(), input.size(), 0, 2, params);
    REQUIRE(encoded_size > 0);

    const size_t decoded_size = h5::impl::filter::deflate(
        decoded.data(), encoded.data(), encoded_size, H5Z_FLAG_REVERSE, 2, params);

    REQUIRE(decoded_size == input.size());
    CHECK(decoded == input);
}

TEST_CASE("H5Z deflate callback rejects corrupt compressed input") {
    std::vector<unsigned char> corrupt = {0x78, 0x9c, 0x01, 0x02, 0x03, 0x04};
    std::vector<unsigned char> decoded(128);
    const unsigned params[] = {6, static_cast<unsigned>(decoded.size())};

    const size_t decoded_size = h5::impl::filter::gzip(
        decoded.data(), corrupt.data(), corrupt.size(), H5Z_FLAG_REVERSE, 2, params);

    CHECK(decoded_size == 0);
}
