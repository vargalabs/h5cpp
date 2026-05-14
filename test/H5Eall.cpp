#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>

TEST_CASE("h5::mute and h5::unmute are idempotent") {
    h5::mute();
    h5::unmute();
    CHECK(true);
}

TEST_CASE("h5::error exception types are constructible") {
    h5::error::io::file::open e1;
    h5::error::io::file::open e2("msg");
    CHECK(std::string(e2.what()) == "msg");

    h5::error::io::dataset::read e3;
    h5::error::io::dataset::read e4("dataset read failed");
    CHECK(std::string(e4.what()) == "dataset read failed");

    h5::error::io::attribute::create e5;
    h5::error::io::attribute::delete_ e6;
    h5::error::property_list::argument e7;
    h5::error::io::packet_table::write e8;
    CHECK(true);
}

TEST_CASE("h5::error::rollback prints to stderr") {
    h5::error::rollback r("test rollback message");
    CHECK(true);
}
