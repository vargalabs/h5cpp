#include <hdf5.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "support/common.hpp"
#include <h5cpp/compat.hpp>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace test {

    struct with_value_type_t {
        using value_type = int;
    };

    struct with_size_type_t {
        using size_type = std::size_t;
    };

    struct without_members_t {};

    struct plain_pod_t {
        int a;
        double b;
    };

    struct non_pod_t {
        non_pod_t() {}
        int x;
    };

    template<class T>
    using value_type_t = typename T::value_type;

    template<class T>
    using size_type_t = typename T::size_type;

    template<class T>
    using begin_expr_t = decltype(std::declval<T&>().begin());

}

TEST_CASE("compat make_index_sequence builds expected sequence") {
    CHECK((std::is_same_v<
        h5::compat::make_index_sequence<0>,
        h5::compat::index_sequence<>
    >));

    CHECK((std::is_same_v<
        h5::compat::make_index_sequence<1>,
        h5::compat::index_sequence<0>
    >));

    CHECK((std::is_same_v<
        h5::compat::make_index_sequence<4>,
        h5::compat::index_sequence<0,1,2,3>
    >));
}

TEST_CASE("compat apply unpacks tuple into callable") {
    auto fn = [](int a, double b, char c) -> herr_t {
        CHECK(a == 7);
        CHECK(b == doctest::Approx(2.5));
        CHECK(c == 'x');
        return 42;
    };

    auto args = std::make_tuple(7, 2.5, 'x');
    CHECK(h5::compat::apply(fn, args) == 42);
}

TEST_CASE("compat detection idiom detects nested types") {
    CHECK(h5::meta::compat::is_detected_v<test::value_type_t, test::with_value_type_t>);
    CHECK_FALSE(h5::meta::compat::is_detected_v<test::value_type_t, test::without_members_t>);

    CHECK(h5::meta::compat::is_detected_v<test::size_type_t, test::with_size_type_t>);
    CHECK_FALSE(h5::meta::compat::is_detected_v<test::size_type_t, test::without_members_t>);
}

TEST_CASE("compat detected_t and fallback types behave correctly") {
    CHECK((std::is_same_v<
        h5::meta::compat::detected_t<test::value_type_t, test::with_value_type_t>,
        int
    >));

    CHECK((std::is_same_v<
        h5::meta::compat::detected_t<test::value_type_t, test::without_members_t>,
        h5::meta::compat::nonesuch
    >));

    CHECK((std::is_same_v<
        h5::meta::compat::detected_or_t<double, test::value_type_t, test::without_members_t>,
        double
    >));
}

TEST_CASE("compat exact and convertible detection work") {
    CHECK(h5::meta::compat::is_detected_exact_v<int, test::value_type_t, test::with_value_type_t>);
    CHECK_FALSE(h5::meta::compat::is_detected_exact_v<double, test::value_type_t, test::with_value_type_t>);

    CHECK(h5::meta::compat::is_detected_convertible_v<long, test::value_type_t, test::with_value_type_t>);
    CHECK_FALSE(h5::meta::compat::is_detected_convertible_v<std::vector<int>, test::value_type_t, test::with_value_type_t>);
}

TEST_CASE("compat is_pod strips cv and reference qualifiers") {
}

TEST_CASE("meta compat aliases mirror compat namespace") {
    CHECK((std::is_same_v<h5::meta::compat::nonesuch, h5::meta::compat::nonesuch>));

    CHECK(h5::meta::compat::is_detected_v<test::value_type_t, test::with_value_type_t>);
    CHECK((std::is_same_v<
        h5::meta::compat::detected_t<test::value_type_t, test::with_value_type_t>,
        int
    >));
}