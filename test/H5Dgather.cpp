#include <hdf5.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/H5Dgather.hpp>
#include <algorithm>
#include <deque>
#include <forward_list>
#include <list>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

TEST_CASE("H5Dgather gathers list int in iteration order") {
    std::list<int> values = {1, 2, 3};
    std::vector<int> elements;
    const int* ptr = h5::gather(values, elements);

    CHECK(ptr == elements.data());
    CHECK(elements == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather gathers forward list int in iteration order") {
    std::forward_list<int> values = {1, 2, 3};
    std::vector<int> elements;
    const int* ptr = h5::gather(values, elements);

    CHECK(ptr == elements.data());
    CHECK(elements == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather gathers deque int in logical order") {
    std::deque<int> values = {1, 2, 3};
    std::vector<int> elements;
    const int* ptr = h5::gather(values, elements);

    CHECK(ptr == elements.data());
    CHECK(elements == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather gathers set int in sorted order") {
    std::set<int> values = {3, 1, 2};
    std::vector<int> elements;
    const int* ptr = h5::gather(values, elements);

    CHECK(ptr == elements.data());
    CHECK(elements == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather gathers unordered set int without sequence assumption") {
    std::unordered_set<int> values = {3, 1, 2};
    std::vector<int> elements;
    const int* ptr = h5::gather(values, elements);
    std::set<int> gathered(elements.begin(), elements.end());

    CHECK(ptr == elements.data());
    CHECK(gathered == std::set<int>{1, 2, 3});
}

TEST_CASE("H5Dgather materializes list int from flat buffer") {
    const int values[] = {1, 2, 3};
    std::list<int> out;

    h5::materialize(out, values, 3);

    CHECK(std::vector<int>(out.begin(), out.end()) == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather materializes forward list int from flat buffer") {
    const int values[] = {1, 2, 3};
    std::forward_list<int> out;

    h5::materialize(out, values, 3);

    CHECK(std::vector<int>(out.begin(), out.end()) == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather materializes deque int from flat buffer") {
    const int values[] = {1, 2, 3};
    std::deque<int> out;

    h5::materialize(out, values, 3);

    CHECK(std::vector<int>(out.begin(), out.end()) == std::vector<int>{1, 2, 3});
}

TEST_CASE("H5Dgather materializes set int from flat buffer") {
    const int values[] = {3, 1, 2};
    std::set<int> out;

    h5::materialize(out, values, 3);

    CHECK(out == std::set<int>{1, 2, 3});
}

TEST_CASE("H5Dgather materializes unordered set int from flat buffer") {
    const int values[] = {3, 1, 2};
    std::unordered_set<int> out;

    h5::materialize(out, values, 3);

    CHECK(out == std::unordered_set<int>{1, 2, 3});
}

TEST_CASE("H5Dgather preserves vector string pointer table gather") {
    std::vector<std::string> values = {"alpha", "beta", "gamma"};
    std::vector<const char*> ptrs;
    const char** table = h5::gather(values, ptrs);

    CHECK(table == ptrs.data());
    REQUIRE(ptrs.size() == values.size());
    CHECK(std::string(table[0]) == "alpha");
    CHECK(std::string(table[1]) == "beta");
    CHECK(std::string(table[2]) == "gamma");
    CHECK(table[0] == values[0].data());
    CHECK(table[1] == values[1].data());
    CHECK(table[2] == values[2].data());
}
