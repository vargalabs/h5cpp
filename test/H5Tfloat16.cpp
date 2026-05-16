#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>

#if __cplusplus >= 202302L && defined(__STDCPP_FLOAT16_T__)
#include <h5cpp/all>
#include <stdfloat>
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

TEST_CASE("[#224] std::float16_t scalar dataset create/write/read round-trip") {
    temp_file tf("h5tfloat16_scalar.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    // Use a single-element vector to exercise scalar-value write/read
    std::vector<std::float16_t> written = { static_cast<std::float16_t>(1.5f) };
    h5::write(fd, "scalar_f16", written);
    auto read = h5::read<std::vector<std::float16_t>>(fd, "scalar_f16");

    REQUIRE(read.size() == 1);
    CHECK(static_cast<float>(read[0]) == doctest::Approx(static_cast<float>(written[0])));
}

TEST_CASE("[#224] std::vector<std::float16_t> dataset round-trip (N=8)") {
    temp_file tf("h5tfloat16_vec.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);

    std::vector<std::float16_t> written = {
        static_cast<std::float16_t>(0.0f),
        static_cast<std::float16_t>(1.0f),
        static_cast<std::float16_t>(-1.0f),
        static_cast<std::float16_t>(0.5f),
        static_cast<std::float16_t>(2.0f),
        static_cast<std::float16_t>(-2.0f),
        static_cast<std::float16_t>(3.14f),
        static_cast<std::float16_t>(-0.125f)
    };
    h5::write(fd, "vec_f16", written);
    auto read = h5::read<std::vector<std::float16_t>>(fd, "vec_f16");

    REQUIRE(read.size() == written.size());
    for (std::size_t i = 0; i < written.size(); ++i) {
        CHECK(static_cast<float>(read[i]) == doctest::Approx(static_cast<float>(written[i])).epsilon(0.01));
    }
}

TEST_CASE("[#224] std::float16_t attribute write/read via h5::awrite/h5::aread") {
    temp_file tf("h5tfloat16_attr.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    auto ds = h5::create<float>(fd, "dataset", h5::current_dims_t{1});

    std::float16_t written = static_cast<std::float16_t>(3.14f);
    h5::awrite(ds, "f16_attr", written);
    auto read = h5::aread<std::float16_t>(ds, "f16_attr");

    CHECK(static_cast<float>(read) == doctest::Approx(static_cast<float>(written)).epsilon(0.01));
}

#endif
