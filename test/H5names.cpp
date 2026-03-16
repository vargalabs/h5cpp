
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>

#include "support/common.hpp"
#include "h5cpp/H5names.hpp"

#include <complex>
#include <string>
#include <type_traits>

/*
namespace h5 {
    template <>
    struct name<ns::pod_t> {
        static constexpr char const* value = "ns::pod_t";
    };
}
*/

TEST_CASE("H5names primitive specializations") {
    CHECK(std::string(h5::name<char const*>::value) == "const char*");
    CHECK(std::string(h5::name<std::complex<float>>::value) == "complex<float>");
    CHECK(std::string(h5::name<std::complex<double>>::value) == "complex<double>");
    //CHECK(std::string(h5::name<ns::pod_t>::value) == "ns::pod_t");
}

TEST_CASE("H5names generic complex fallback") {
    CHECK(std::string(h5::name<std::complex<int>>::value) == "complex<int>");
    CHECK(std::string(h5::name<std::complex<short>>::value) == "complex<short>");
}

TEST_CASE("H5names array extent formatting is non-empty") {
    using array_t = int[7];

    CHECK(std::string(h5::name<array_t>::value).size() > 0);
}

TEST_CASE("H5names array extent formatting is non-empty") {
    using array_t = int[7];

    CHECK(std::string(h5::name<array_t>::value).size() > 0);
}
