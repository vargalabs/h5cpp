#pragma once

#include <hdf5.h>
#include <doctest/all>
#include <h5cpp/core>

#include "names.hpp"
#include "types.hpp"

namespace h5::test {
    template <class T, class...>
    struct name {
        static constexpr char const* value = "T";
    };
}
