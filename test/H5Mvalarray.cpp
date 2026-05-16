#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <valarray>
#include <filesystem>
#include <cmath>

namespace {
    struct temp_file {
        std::string path;
        explicit temp_file(const std::string& name)
            : path((std::filesystem::temp_directory_path() / name).string()) {}
        ~temp_file() {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    };
}

TEST_CASE("[#149] std::valarray<double> write/read round-trip") {
    temp_file tf("h5mvalarray_double.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::valarray<double> written(10);
    for (std::size_t i = 0; i < written.size(); ++i)
        written[i] = static_cast<double>(i) * 0.5;

    h5::write(fd, "valarray_double", written);
    auto read = h5::read<std::valarray<double>>(fd, "valarray_double");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i)
        CHECK(read[i] == doctest::Approx(written[i]));
}

TEST_CASE("[#149] std::valarray<float> write/read round-trip") {
    temp_file tf("h5mvalarray_float.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::valarray<float> written(8);
    for (std::size_t i = 0; i < written.size(); ++i)
        written[i] = static_cast<float>(i) * 1.5f;

    h5::write(fd, "valarray_float", written);
    auto read = h5::read<std::valarray<float>>(fd, "valarray_float");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i)
        CHECK(read[i] == doctest::Approx(written[i]));
}

TEST_CASE("[#149] std::valarray<int> write/read round-trip") {
    temp_file tf("h5mvalarray_int.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::valarray<int> written = {1, 2, 3, 4, 5, 6, 7, 8};
    h5::write(fd, "valarray_int", written);
    auto read = h5::read<std::valarray<int>>(fd, "valarray_int");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i)
        CHECK(read[i] == written[i]);
}

TEST_CASE("[#149] std::valarray<double> extent matches size()") {
    temp_file tf("h5mvalarray_extent.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    constexpr std::size_t N = 16;
    std::valarray<double> written(N);
    for (std::size_t i = 0; i < N; ++i) written[i] = static_cast<double>(i);

    h5::write(fd, "va", written);

    auto ds = h5::open(fd, "va", H5F_ACC_RDONLY);
    auto sp = h5::get_space(ds);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(static_cast<hid_t>(sp), dims, nullptr);
    CHECK(dims[0] == N);
}
