#include <hdf5.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/H5Tmeta.hpp>
#include "support/types.hpp"
#include "support/fixture.hpp"

using sr_t = h5::meta::storage_representation_t;

TEST_CASE_TEMPLATE_DEFINE("[#114] h5 round-trip", T, io_round_trip) {
    constexpr auto rep  = h5::meta::storage_representation_v<T>;
    constexpr auto type_name = h5::test::name<T>::value;
    h5::test::file_fixture_t f("test-io.h5");

    if constexpr (rep == sr_t::scalar) {
        WARN_MESSAGE(false, "not yet implemented: scalar round-trip for " << type_name);
    } else if constexpr (rep == sr_t::c_array) {
        WARN_MESSAGE(false, "not yet implemented: c_array round-trip for " << type_name);
    } else if constexpr (rep == sr_t::linear_value_dataset) {
        WARN_MESSAGE(false, "not yet implemented: linear_value_dataset round-trip for " << type_name);
    } else if constexpr (rep == sr_t::ragged_vlen_dataset) {
        WARN_MESSAGE(false, "not yet implemented: ragged_vlen_dataset round-trip for " << type_name);
    } else if constexpr (rep == sr_t::fixed_inner_extent_dataset) {
        WARN_MESSAGE(false, "not yet implemented: fixed_inner_extent_dataset round-trip for " << type_name);
    } else if constexpr (rep == sr_t::key_value_dataset) {
        WARN_MESSAGE(false, "not yet implemented: key_value_dataset round-trip for " << type_name);
    } else if constexpr (rep == sr_t::vlen_text_dataset) {
        WARN_MESSAGE(false, "not yet implemented: vlen_text_dataset round-trip for " << type_name);
    }
}
TEST_CASE_TEMPLATE_APPLY(io_round_trip, h5::test::io_types);
