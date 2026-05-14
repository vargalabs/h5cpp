#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <filesystem>
#include <vector>
#include <string>

namespace {
    struct temp_file {
        std::string path;
        explicit temp_file(const std::string& name) {
            path = (std::filesystem::temp_directory_path() / name).string();
        }
        ~temp_file() {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    };
}

// NOTE: True scalar dataset round-trip is currently blocked on issue #89.
// h5::write(fd,"ds",42) creates an H5S_SCALAR space, but h5::read<T> always
// invokes H5Sselect_hyperslab which fails for scalars.  Until that is fixed
// we use a 1-D single-element vector to exercise the write/read plumbing.
TEST_CASE("[#132] scalar-like int round-trip") {
    temp_file tf("h5dio_scalar_int.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    std::vector<int> data = {42};
    h5::write(fd, "dataset", data);
    auto read_vec = h5::read<std::vector<int>>(fd, "dataset");
    CHECK(read_vec == data);
}

TEST_CASE("[#132] 1D vector<double> round-trip") {
    temp_file tf("h5dio_vector_double.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    std::vector<double> data = {1.0, 2.0, 3.0};
    h5::write(fd, "vector", data);
    auto read_vec = h5::read<std::vector<double>>(fd, "vector");
    CHECK(read_vec == data);
}

// NOTE: h5::read<std::string> does not compile because impl::decay<std::string>
// yields const char*, routing into the generic (non-string) read path.
// The working API for string datasets is std::vector<std::string>.
TEST_CASE("[#132] string round-trip") {
    temp_file tf("h5dio_string.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    std::vector<std::string> data = {"hello"};
    h5::write(fd, "string", data);
    auto read_vec = h5::read<std::vector<std::string>>(fd, "string");
    CHECK(read_vec == data);
}

TEST_CASE("[#132] file create + open round-trip") {
    temp_file tf("h5dio_reopen.h5");
    {
        auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
        h5::write(fd, "dataset", std::vector<int>{123});
    }
    auto fd = h5::open(tf.path, H5F_ACC_RDWR);
    auto read_vec = h5::read<std::vector<int>>(fd, "dataset");
    REQUIRE(read_vec.size() == 1);
    CHECK(read_vec[0] == 123);
}

TEST_CASE("[#132] read non-existent dataset throws") {
    temp_file tf("h5dio_nonexistent.h5");
    auto fd = h5::create(tf.path, H5F_ACC_TRUNC);
    CHECK_THROWS_AS(h5::read<int>(fd, "nonexistent"), h5::error::io::dataset::open);
}

TEST_CASE("[#178] read non-existent file throws") {
    CHECK_THROWS_AS(h5::open("nonexistent_file_xyz.h5", H5F_ACC_RDONLY), h5::error::io::file::open);
}

TEST_CASE("[#178] invalid property list throws on file create") {
    h5::fcpl_t bad_fcpl{H5I_UNINIT};
    CHECK_THROWS_AS(h5::create("test_bad_fcpl.h5", H5F_ACC_TRUNC, bad_fcpl), h5::error::io::file::create);
    std::error_code ec;
    std::filesystem::remove("test_bad_fcpl.h5", ec);
}
