#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <complex>
#include <vector>
#include <filesystem>

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

TEST_CASE("[#43] std::complex<double> scalar dataset round-trip") {
    temp_file tf("h5tcomplex_scalar.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::vector<std::complex<double>> written = { {1.5, -2.5} };
    h5::write(fd, "scalar_cplx", written);
    auto read = h5::read<std::vector<std::complex<double>>>(fd, "scalar_cplx");

    REQUIRE(read.size() == 1);
    CHECK(read[0].real() == doctest::Approx(1.5));
    CHECK(read[0].imag() == doctest::Approx(-2.5));
}

TEST_CASE("[#43] std::vector<std::complex<float>> round-trip") {
    temp_file tf("h5tcomplex_float_vec.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::vector<std::complex<float>> written = {
        {1.0f, 2.0f}, {-3.5f, 4.25f}, {0.0f, -1.0f}
    };
    h5::write(fd, "float_vec", written);
    auto read = h5::read<std::vector<std::complex<float>>>(fd, "float_vec");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i) {
        CHECK(read[i].real() == doctest::Approx(written[i].real()));
        CHECK(read[i].imag() == doctest::Approx(written[i].imag()));
    }
}

TEST_CASE("[#43] std::vector<std::complex<double>> round-trip") {
    temp_file tf("h5tcomplex_double_vec.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::vector<std::complex<double>> written = {
        {1.0, 2.0}, {-3.5, 4.25}, {0.0, -1.0}, {1e10, -1e-10}
    };
    h5::write(fd, "double_vec", written);
    auto read = h5::read<std::vector<std::complex<double>>>(fd, "double_vec");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i) {
        CHECK(read[i].real() == doctest::Approx(written[i].real()));
        CHECK(read[i].imag() == doctest::Approx(written[i].imag()));
    }
}

TEST_CASE("[#43] std::complex<double> attribute round-trip") {
    temp_file tf("h5tcomplex_attr.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    auto ds = h5::create<double>(fd, "dataset", h5::current_dims_t{1});

    std::vector<std::complex<double>> written = { {3.14, -2.72} };
    h5::awrite(ds, "cplx_attr", written);
    auto read = h5::aread<std::vector<std::complex<double>>>(ds, "cplx_attr");

    REQUIRE(read.size() == 1);
    CHECK(read[0].real() == doctest::Approx(3.14));
    CHECK(read[0].imag() == doctest::Approx(-2.72));
}
