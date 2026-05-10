/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

#include <string>
#include <type_traits>

namespace h5 {
    // Delete an attribute from an HDF5 object by name
    template<class HID_T>
    inline std::enable_if_t<h5::impl::is_valid_attr<HID_T>::value, void>
    adelete(const HID_T& parent, const std::string& name) {
        H5CPP_CHECK_NZ(
            H5Adelete(static_cast<hid_t>(parent), name.c_str()),
            h5::error::io::attribute::delete_,
            "couldn't delete attribute");
    }
}

