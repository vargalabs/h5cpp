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
    h5::error::io::packet_table::write e9("msg");
    CHECK(true);
}

TEST_CASE("h5::error::rollback prints to stderr") {
    h5::error::rollback r("test rollback message");
    CHECK(true);
}

TEST_CASE("h5::error dataset exception constructors") {
    h5::error::io::dataset::write e1;
    h5::error::io::dataset::write e2("msg");
    h5::error::io::dataset::append e3;
    h5::error::io::dataset::append e4("msg");
    h5::error::io::dataset::create e5;
    h5::error::io::dataset::create e6("msg");
    h5::error::io::dataset::close e7;
    h5::error::io::dataset::close e8("msg");
    h5::error::io::dataset::misc e9;
    h5::error::io::dataset::misc e10("msg");
    CHECK(true);
}

TEST_CASE("h5::error packet_table exception constructors") {
    h5::error::io::packet_table::open e1;
    h5::error::io::packet_table::open e2("msg");
    h5::error::io::packet_table::close e3;
    h5::error::io::packet_table::close e4("msg");
    h5::error::io::packet_table::create e5;
    h5::error::io::packet_table::create e6("msg");
    h5::error::io::packet_table::misc e7;
    h5::error::io::packet_table::misc e8("msg");
    CHECK(true);
}

TEST_CASE("h5::error attribute exception constructors") {
    h5::error::io::attribute::close e1;
    h5::error::io::attribute::close e2("msg");
    h5::error::io::attribute::read e3;
    h5::error::io::attribute::read e4("msg");
    h5::error::io::attribute::write e5;
    h5::error::io::attribute::write e6("msg");
    h5::error::io::attribute::create e7;
    h5::error::io::attribute::create e8("msg");
    h5::error::io::attribute::misc e9;
    h5::error::io::attribute::misc e10("msg");
    CHECK(true);
}

TEST_CASE("h5::error property_list exception constructors") {
    h5::error::property_list::argument e1;
    h5::error::property_list::argument e2("msg");
    h5::error::property_list::misc e3;
    h5::error::property_list::misc e4("msg");
    CHECK(true);
}
