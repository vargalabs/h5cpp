#include <hdf5.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Tmeta.hpp>
#include "support/types.hpp"
#include "support/fixture.hpp"
#include "support/generator.hpp"

using sr_t = h5::meta::storage_representation_t;

static constexpr size_t LOWER = 2, UPPER = 100, MIN_SZ = 4, MAX_SZ = 16;

TEST_CASE_TEMPLATE_DEFINE("[#114] h5 round-trip", T, io_round_trip) {
    constexpr auto rep = h5::meta::storage_representation_v<T>;
    const std::string type_name{h5::test::name<T>::value};
    h5::test::file_fixture_t f("test-io.h5");

    if constexpr (rep == sr_t::scalar) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER);
        WARN_MESSAGE(false,
            "scalar round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::c_array) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER);
        WARN_MESSAGE(false,
            "c_array round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::linear_value_dataset) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER, MIN_SZ, MAX_SZ);
        WARN_MESSAGE(false,
            "linear_value_dataset round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::ragged_vlen_dataset) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER, MIN_SZ, MAX_SZ);
        WARN_MESSAGE(false,
            "ragged_vlen_dataset round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::fixed_inner_extent_dataset) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER, MIN_SZ, MAX_SZ);
        WARN_MESSAGE(false,
            "fixed_inner_extent_dataset round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::key_value_dataset) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER, MIN_SZ, MAX_SZ);
        WARN_MESSAGE(false,
            "key_value_dataset round-trip not yet implemented (needs #89): " << type_name);

    } else if constexpr (rep == sr_t::vlen_text_dataset) {
        T obj{};
        h5::test::fill(obj, LOWER, UPPER, MIN_SZ, MAX_SZ);
        WARN_MESSAGE(false,
            "vlen_text_dataset round-trip not yet implemented (needs #89): " << type_name);
    }
}
TEST_CASE_TEMPLATE_APPLY(io_round_trip, h5::test::io_types);
