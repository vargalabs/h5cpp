#pragma once

#include <doctest/all>
#include <h5cpp/core>

//#include "pod_t.hpp"
#include "names.hpp"
#include "types.hpp"

namespace h5::test {
    template <class T, class...>
    struct name {
        static constexpr char const* value = "T";
    };
}

#include "types.hpp"
#include "names.hpp"
