#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <vector>
#include <forward_list>
#include <sstream>
#include <string>
#include "support/fixture.hpp"

TEST_CASE("packet table append scalar values") {
    h5::test::file_fixture_t f("test-pt-append.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::pt_t pt(ds);
    for (int i = 0; i < 25; ++i)
        h5::append(pt, i);
    h5::flush(pt);
    CHECK(true);
}

TEST_CASE("packet table append std::vector chunk") {
    h5::test::file_fixture_t f("test-pt-vector.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0, 5},
        h5::max_dims_t{H5S_UNLIMITED, 5}, h5::chunk{2, 5});
    h5::pt_t pt(ds);
    std::vector<int> chunk = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    h5::append(pt, chunk);
    h5::flush(pt);
    CHECK(true);
}

TEST_CASE("packet table partial flush fills remainder with fill value") {
    h5::test::file_fixture_t f("test-pt-partial.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    {
        h5::pt_t pt(ds);
        for (int i = 0; i < 3; ++i)
            h5::append(pt, i + 1);
        h5::flush(pt);
    }
    auto readback = h5::read<std::vector<int>>(f.fd, "ds");
    REQUIRE(readback.size() == 10);  // partial chunk padded to full chunk size
    for (int i = 0; i < 3; ++i)
        CHECK(readback[i] == i + 1);
    for (int i = 3; i < 10; ++i)
        CHECK(readback[i] == 0);  // fill value pads remainder
}

TEST_CASE("packet table auto-flush on destruction") {
    h5::test::file_fixture_t f("test-pt-dtor.h5");
    {
        h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
            h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
        h5::pt_t pt(ds);
        for (int i = 0; i < 3; ++i)
            h5::append(pt, i + 1);
        // pt destructor flushes partial chunk, then ds closes
    }
    auto readback = h5::read<std::vector<int>>(f.fd, "ds");
    REQUIRE(readback.size() == 10);  // partial chunk padded to full chunk size
    for (int i = 0; i < 3; ++i)
        CHECK(readback[i] == i + 1);
    for (int i = 3; i < 10; ++i)
        CHECK(readback[i] == 0);  // fill value pads remainder
}

TEST_CASE("packet table string append and flush") {
    h5::test::file_fixture_t f("test-pt-string.h5");
    h5::mute();
    {
        h5::ds_t ds = h5::create<std::string>(f.fd, "ds", h5::current_dims_t{0},
            h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{5});
        h5::pt_t pt(ds);
        std::vector<std::string> strings = {"alpha", "beta", "gamma", "delta", "epsilon"};
        for (const auto& s : strings)
            h5::append(pt, s);
        h5::flush(pt);
    }
    h5::unmute();
    CHECK(true);
}

TEST_CASE("packet table const char* append") {
    h5::test::file_fixture_t f("test-pt-charptr.h5");
    h5::mute();
    {
        h5::ds_t ds = h5::create<std::string>(f.fd, "ds", h5::current_dims_t{0},
            h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{3});
        h5::pt_t pt(ds);
        h5::append(pt, "hello");
        h5::append(pt, "world");
        h5::append(pt, "test");
        h5::flush(pt);
    }
    h5::unmute();
    CHECK(true);
}

TEST_CASE("packet table copy constructor") {
    h5::test::file_fixture_t f("test-pt-copy.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::pt_t pt1(ds);
    for (int i = 0; i < 5; ++i)
        h5::append(pt1, i);

    h5::pt_t pt2(pt1); // copy ctor
    h5::flush(pt2);
    CHECK(true);
}

TEST_CASE("packet table move assignment") {
    h5::test::file_fixture_t f("test-pt-move.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::pt_t pt1(ds);
    for (int i = 0; i < 5; ++i)
        h5::append(pt1, i);

    h5::pt_t pt2;
    pt2 = std::move(pt1); // move assignment

    h5::mute();
    h5::flush(pt2);
    h5::unmute();
    CHECK(true);
}

TEST_CASE("packet table self move assignment") {
    h5::test::file_fixture_t f("test-pt-self-move.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::pt_t pt(ds);
    for (int i = 0; i < 5; ++i)
        h5::append(pt, i);
    pt = std::move(pt); // self move assignment
    h5::flush(pt);
    CHECK(true);
}

TEST_CASE("packet table output stream operator") {
    h5::test::file_fixture_t f("test-pt-ostream.h5");
    h5::ds_t ds = h5::create<int>(f.fd, "ds", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{10});
    h5::pt_t pt(ds);
    std::ostringstream oss;
    oss << pt;
    CHECK(oss.str().find("packet table:") != std::string::npos);
}

TEST_CASE("packet table invalid handle conversion ctor") {
    h5::pt_t pt(h5::ds_t{H5I_UNINIT});
    CHECK(true);
}

TEST_CASE("packet table output stream for invalid handle") {
    h5::pt_t pt;
    std::ostringstream oss;
    oss << pt;
    CHECK(oss.str().find("H5I_UNINIT") != std::string::npos);
}

TEST_CASE("[#232] std::forward_list<int> append streams elements into chunked dataset") {
    h5::test::file_fixture_t f("test-pt-fwdlist.h5");
    // forward_list is append/view only — h5::write/read intentionally unsupported.
    // Each element is streamed one-by-one into the packet table (partial chunk flushed explicitly).
    std::forward_list<int> src = {10, 20, 30, 40, 50};
    constexpr std::size_t N = 5;

    h5::ds_t ds = h5::create<int>(f.fd, "fwdlist", h5::current_dims_t{0},
        h5::max_dims_t{H5S_UNLIMITED}, h5::chunk{5});  // chunk == list size: auto-flush at boundary
    {
        h5::pt_t pt(ds);
        h5::append(pt, src);  // iterates element-by-element via iterator dispatch
        h5::flush(pt);        // flush partial chunk to file
    }

    auto readback = h5::read<std::vector<int>>(f.fd, "fwdlist");
    REQUIRE(readback.size() == N);
    const std::vector<int> expected(src.begin(), src.end());
    CHECK(readback == expected);
}
