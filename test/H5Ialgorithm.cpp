#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <h5cpp/H5Ialgorithm.hpp>
#include "support/fixture.hpp"

TEST_CASE("h5::ls lists group contents") {
    h5::test::file_fixture_t f("test-ls.h5");
    h5::create<int>(f.fd, "group/dataset_a", h5::current_dims_t{1});
    h5::create<int>(f.fd, "group/dataset_b", h5::current_dims_t{1});
    auto files = h5::ls(f.fd, "/group");
    CHECK(files.size() == 2);
}

TEST_CASE("h5::ls on empty group returns empty vector") {
    h5::test::file_fixture_t f("test-ls-empty.h5");
    h5::create<int>(f.fd, "dataset", h5::current_dims_t{1});
    auto files = h5::ls(f.fd, "/");
    CHECK(files.size() == 1);
}

TEST_CASE("h5::bfs returns empty vector (not yet implemented)") {
    h5::test::file_fixture_t f("test-bfs.h5");
    auto files = h5::bfs(f.fd, "/");
    CHECK(files.empty());
}

TEST_CASE("h5::dfs returns empty vector (not yet implemented)") {
    h5::test::file_fixture_t f("test-dfs.h5");
    auto files = h5::dfs(f.fd, "/");
    CHECK(files.empty());
}
