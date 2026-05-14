// SPDX-License-Identifier: MIT
// This file is part of H5CPP.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <filesystem>
#include <vector>

TEST_CASE("H5R reference adapter preserves dataset-region references") {
    const char* filename = "test-h5r-reference.h5";
    std::filesystem::remove(filename);

    h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
    std::vector<int> data(16, 0);
    h5::write(fd, "data", data, h5::current_dims{4,4});

    h5::reference_t ref = h5::reference(fd, "data", h5::offset{1,1}, h5::count{2,2});
    h5::reference_t ref_copy;
    CHECK(h5::impl::reference::copy(ref_copy, ref) >= 0);
    CHECK(h5::impl::reference::equal(ref, ref_copy) > 0);

    std::vector<h5::reference_t> refs_data{ref_copy};
    h5::write(fd, "refs", refs_data, h5::current_dims{1});
    h5::ds_t refs = h5::open(fd, "refs");
    auto stored = h5::read<std::vector<h5::reference_t>>(refs);
    REQUIRE(stored.size() == 1);

    int patch[] = {1,2,3,4};
    h5::exp::write(refs, stored[0], patch);

    auto out = h5::read<std::vector<int>>(fd, "data");
    REQUIRE(out.size() == 16);
    CHECK(out[5] == 1);
    CHECK(out[6] == 2);
    CHECK(out[9] == 3);
    CHECK(out[10] == 4);

    CHECK(h5::impl::reference::destroy(ref) >= 0);
    CHECK(h5::impl::reference::destroy(ref_copy) >= 0);
    for(auto& item: stored)
        CHECK(h5::impl::reference::destroy(item) >= 0);
    std::filesystem::remove(filename);
}
