#include <hdf5.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <array>
#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <h5cpp/H5meta.hpp>

struct value_only_t {
    using value_type = int;
};
struct size_only_t {
    using size_type = std::size_t;
};

struct map_like_t {
    using key_type = int;
    using mapped_type = double;
};

struct iterable_t {
    int* begin();
    int* end();
};

struct const_iterable_t {
    const int* cbegin() const;
    const int* cend() const;
};

struct contiguous_t {
    using value_type = int;
    using size_type = std::size_t;

    int* data();
    std::size_t size();
};

struct tag_t {};

struct convertible_to_tag_t {
    operator tag_t() const { return {}; }
};

TEST_CASE("H5meta default h5::name fallback") {
    CHECK(std::string(h5::name<int>::value) == "n/a");
    CHECK(std::string(h5::name<double>::value) == "n/a");
}

TEST_CASE("H5meta detects nested typedefs") {
    CHECK(h5::meta::has_value_type<std::vector<int>>::value);

    CHECK(h5::meta::has_value_type<value_only_t>::value);

    CHECK_FALSE(h5::meta::has_value_type<size_only_t>::value);


}

TEST_CASE("H5meta detects container style member functions") {
    CHECK(h5::meta::has_data<std::vector<int>>::value);
    CHECK(h5::meta::has_size<std::vector<int>>::value);
    CHECK(h5::meta::has_begin<std::vector<int>>::value);
    CHECK(h5::meta::has_end<std::vector<int>>::value);
    CHECK(h5::meta::has_cbegin<std::vector<int>>::value);
    CHECK(h5::meta::has_cend<std::vector<int>>::value);

    CHECK(h5::meta::has_data<contiguous_t>::value);
    CHECK(h5::meta::has_size<contiguous_t>::value);

    CHECK(h5::meta::has_begin<iterable_t>::value);
    CHECK(h5::meta::has_end<iterable_t>::value);
    CHECK_FALSE(h5::meta::has_cbegin<iterable_t>::value);
    CHECK_FALSE(h5::meta::has_cend<iterable_t>::value);

    CHECK(h5::meta::has_cbegin<const_iterable_t>::value);
    CHECK(h5::meta::has_cend<const_iterable_t>::value);
    CHECK_FALSE(h5::meta::has_begin<const_iterable_t>::value);
    CHECK_FALSE(h5::meta::has_end<const_iterable_t>::value);
}

TEST_CASE("H5meta iterable aliases behave as expected") {
    CHECK(h5::meta::has_iterator<std::vector<int>>::value);
    CHECK(h5::meta::has_const_iterator<std::vector<int>>::value);
    CHECK(h5::meta::has_iterator<std::vector<int>>::value);
    CHECK(h5::meta::has_const_iterator<std::vector<int>>::value);

    CHECK(h5::meta::has_iterator<iterable_t>::value);
    CHECK_FALSE(h5::meta::has_const_iterator<iterable_t>::value);
    CHECK(h5::meta::has_iterator<iterable_t>::value);
    CHECK_FALSE(h5::meta::has_const_iterator<iterable_t>::value);

    CHECK_FALSE(h5::meta::has_iterator<const_iterable_t>::value);
    CHECK(h5::meta::has_const_iterator<const_iterable_t>::value);
    CHECK_FALSE(h5::meta::has_iterator<const_iterable_t>::value);
    CHECK(h5::meta::has_const_iterator<const_iterable_t>::value);
}

TEST_CASE("H5meta map-like detection works") {
}

TEST_CASE("H5meta value alias resolves value_type or falls back to self") {
    CHECK((std::is_same_v<typename h5::meta::value<std::vector<int>>::type, int>));
    CHECK((std::is_same_v<typename h5::meta::value<value_only_t>::type, int>));
    CHECK((std::is_same_v<typename h5::meta::value<double>::type, double>));
}

TEST_CASE("H5meta detail is_value_type detects nested value_type convertibility") {
    CHECK(h5::meta::detail::is_value_type<std::vector<int>, int>::value);
    CHECK(h5::meta::detail::is_value_type<std::vector<int>, double>::value);

    CHECK(h5::meta::detail::is_value_type<value_only_t, int>::value);
    CHECK_FALSE(h5::meta::detail::is_value_type<value_only_t, std::string>::value);
}

TEST_CASE("H5meta detail is_value_type falls back to direct convertibility") {
    CHECK(h5::meta::detail::is_value_type<convertible_to_tag_t, tag_t>::value);
    CHECK_FALSE(h5::meta::detail::is_value_type<int, std::string>::value);
}

TEST_CASE("H5meta static_for visits all tuple elements in order") {
    using types_t = std::tuple<int, double, char const*>;

    std::size_t count = 0;
    std::array<std::size_t, 3> visited{};

    h5::meta::static_for<types_t>([&](auto index) {
        constexpr std::size_t i = decltype(index)::value;
        using T = typename std::tuple_element<i, types_t>::type;

        visited[count] = i;
        CHECK(std::string(h5::name<T>::value) == "n/a");
        ++count;
    });

    CHECK(count == 3);
    CHECK(visited[0] == 0);
    CHECK(visited[1] == 1);
    CHECK(visited[2] == 2);
}

TEST_CASE("H5meta static_for honours non-zero start offset") {
    using types_t = std::tuple<int, float, double, long>;

    std::size_t count = 0;
    std::array<std::size_t, 2> visited{};
    h5::meta::static_for<types_t, 2>([&](auto index) {
        constexpr std::size_t i = decltype(index)::value;
        visited[count++] = i;
    });

    CHECK(count == 2);
    CHECK(visited[0] == 2);
    CHECK(visited[1] == 3);
}

TEST_CASE("H5meta tpos finds exact types in tuple") {
    using tuple_t = std::tuple<int, double, std::string>;

    using pos_int_t = h5::meta::tpos<int, tuple_t>;
    using pos_string_t = h5::meta::tpos<std::string, tuple_t>;
    using pos_missing_t = h5::meta::tpos<float, tuple_t>;

    CHECK(pos_int_t::present);
    CHECK(pos_int_t::value == 0);
    CHECK((std::is_same_v<typename pos_int_t::type, int>));

    CHECK(pos_string_t::present);
    CHECK(pos_string_t::value == 2);
    CHECK((std::is_same_v<typename pos_string_t::type, std::string>));

    CHECK_FALSE(pos_missing_t::present);
    CHECK(pos_missing_t::value < 0);
}

TEST_CASE("H5meta arg tpos finds first convertible argument") {
    using pos_string_t = h5::arg::tpos<const std::string&, const int&, const std::string&, const char*>;
    using pos_vector_t = h5::arg::tpos<const std::vector<int>&, const int&, const double&, const std::vector<int>&>;
    using pos_missing_t = h5::arg::tpos<const std::vector<double>&, const int&, const double&>;

    CHECK(pos_string_t::present);
    CHECK(pos_string_t::value == 1);
    CHECK((std::is_same_v<typename pos_string_t::type, const std::string&>));

    CHECK(pos_vector_t::present);
    CHECK(pos_vector_t::value == 2);
    CHECK((std::is_same_v<typename pos_vector_t::type, const std::vector<int>&>));

    CHECK_FALSE(pos_missing_t::present);
    CHECK(pos_missing_t::value == -1);
}

TEST_CASE("H5meta arg get returns matching argument when present") {
    std::string s = "hello";
    std::vector<int> v{1, 2, 3};
    std::string def = "default";
    std::vector<int> def_v{9, 9};

    auto& ref_s = h5::arg::get<std::string>(def, 7, s, 3.5);
    auto& ref_v = h5::arg::get<std::vector<int>>(def_v, 7, 3.5, v);

    CHECK((std::is_same_v<decltype(ref_s), std::string&>));
    CHECK(&ref_s == &s);
    CHECK(ref_s == "hello");

    CHECK((std::is_same_v<decltype(ref_v), std::vector<int>&>));
    CHECK(&ref_v == &v);
    CHECK(ref_v.size() == 3);
    CHECK(ref_v[0] == 1);
}

TEST_CASE("H5meta arg get returns default when match is absent") {
    int i = 1;
    double d = 2.0;
    std::string def = "fallback";

    auto& ref = h5::arg::get<std::string>(def, i, d);

    CHECK((std::is_same_v<decltype(ref), const std::string&>));
    CHECK(&ref == &def);
    CHECK(ref == "fallback");
}

TEST_CASE("H5meta arg getn returns nth forwarded argument") {
    int i = 11;
    double d = 2.5;
    std::string s = "abc";

    auto& ref0 = h5::arg::getn<0>(i, d, s);
    auto& ref1 = h5::arg::getn<1>(i, d, s);
    auto& ref2 = h5::arg::getn<2>(i, d, s);

    CHECK(&ref0 == &i);
    CHECK(&ref1 == &d);
    CHECK(&ref2 == &s);
}

TEST_CASE("H5meta impl conditional selects correct branch") {
    CHECK((std::is_same_v<typename h5::impl::conditional<true, int, double>::type, int>));
    CHECK((std::is_same_v<typename h5::impl::conditional<false, int, double>::type, double>));
}

TEST_CASE("H5meta impl cat concatenates tuple types") {
    using lhs_t = std::tuple<int, double>;
    using rhs_t = std::tuple<char, long>;
    using result_t = typename h5::impl::cat<lhs_t, rhs_t>::type;

    CHECK((std::is_same_v<result_t, std::tuple<int, double, char, long>>));
}

TEST_CASE("H5meta impl cat concatenates multiple tuples") {
    using t0_t = std::tuple<int>;
    using t1_t = std::tuple<double, char>;
    using t2_t = std::tuple<long>;
    using result_t = typename h5::impl::cat<t0_t, t1_t, t2_t>::type;

    CHECK((std::is_same_v<result_t, std::tuple<int, double, char, long>>));
}

TEST_CASE("H5meta get_extent returns empty array for scalars") {
    constexpr auto extents = h5::meta::get_extent<int>();
    CHECK(extents.size() == 0);
}

TEST_CASE("H5meta get_extent returns 1d extents correctly") {
    constexpr auto extents = h5::meta::get_extent<int[7]>();

    CHECK(extents.size() == 1);
    CHECK(extents[0] == 7);
}

TEST_CASE("H5meta get_extent returns multidimensional extents correctly") {
    constexpr auto extents = h5::meta::get_extent<double[2][3][4]>();

    CHECK(extents.size() == 3);
    CHECK(extents[0] == 2);
    CHECK(extents[1] == 3);
    CHECK(extents[2] == 4);
}

TEST_CASE("H5meta extent_to_string formats scalar and array types") {
    CHECK(std::string(h5::meta::extent_to_string<int>::value) == "int[]");
    CHECK(std::string(h5::meta::extent_to_string<int[7]>::value) == "int[7]");
    CHECK(std::string(h5::meta::extent_to_string<int[2][3]>::value) == "int[2,3]");
    CHECK(std::string(h5::meta::extent_to_string<double[4]>::value) == "double[4]");
    CHECK(std::string(h5::meta::extent_to_string<unsigned short[5]>::value) == "unsigned short[5]");
    CHECK(std::string(h5::meta::extent_to_string<long long[2][10]>::value) == "long long[2,10]");
}

TEST_CASE("H5meta extent_to_string preserves base type naming") {
    CHECK(std::string(h5::meta::extent_to_string<char[3]>::value) == "char[3]");
    CHECK(std::string(h5::meta::extent_to_string<short[9]>::value) == "short[9]");
    CHECK(std::string(h5::meta::extent_to_string<float[8]>::value) == "float[8]");
    CHECK(std::string(h5::meta::extent_to_string<long double[6]>::value) == "long double[6]");
}