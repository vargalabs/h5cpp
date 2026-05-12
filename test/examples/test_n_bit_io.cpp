// SPDX-License-Identifier: MIT
// This file is part of H5CPP.
// Converted from examples/datatypes/n-bit.cpp into an I/O round-trip test.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/core>
#include <h5cpp/io>
#include <filesystem>
#include <vector>

// Custom n-bit type definition (extracted from examples/datatypes/n-bit.hpp
// to avoid Eigen dependency in the operator<< overload at bottom of header)
namespace bitstring {
    struct n_bit {
        explicit operator unsigned char() const { return value; }
        explicit operator unsigned int() const { return value; }
        n_bit() = default;
        n_bit(unsigned char value_) : value(value_) {}
        unsigned char value;
    };
}

namespace h5::impl::detail {
    template <> struct hid_t<bitstring::n_bit, H5Tclose, true, true, hdf5::type>
        : public dt_p<bitstring::n_bit> {
        using parent = dt_p<bitstring::n_bit>;
        using dt_p<bitstring::n_bit>::hid_t;
        using hidtype = bitstring::n_bit;
        hid_t() : parent(H5Tcopy(H5T_NATIVE_UCHAR)) {
            H5Tset_precision(handle, 2);
        }
    };
}
namespace h5 {
    template <> struct name<bitstring::n_bit> {
        static constexpr char const* value = "bitstring::n_bit";
    };
}

TEST_CASE("[example] custom n-bit type round-trip") {
    const char* filename = "test_n_bit_io.h5";
    std::filesystem::remove(filename);

    namespace bs = bitstring;

    // BUILD test data using STL vector (method 3 from original example)
    std::vector<bs::n_bit> original(12 * 8);
    for (size_t i = 0; i < original.size(); ++i) {
        original[i] = bs::n_bit(static_cast<unsigned char>(i % 4));
    }

    // WRITE
    {
        h5::fd_t fd = h5::create(filename, H5F_ACC_TRUNC);
        h5::ds_t ds = h5::create<bs::n_bit>(fd, "stl",
            h5::current_dims{12, 8}, h5::max_dims{12, H5S_UNLIMITED},
            h5::chunk{3, 4} | h5::nbit);
        h5::write<bs::n_bit>(ds, original.data(), h5::count{12, 8});
    }

    // READ BACK
    {
        h5::fd_t fd = h5::open(filename, H5F_ACC_RDONLY);
        auto readback = h5::read<std::vector<bs::n_bit>>(fd, "stl");

        CHECK(readback.size() == original.size());
        for (size_t i = 0; i < original.size(); ++i) {
            CHECK(readback[i].value == original[i].value);
        }
    }

    std::filesystem::remove(filename);
}
