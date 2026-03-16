/* SPDX-License-Identifier: MIT
 * This file is part of H5CPP.
 * Copyright (c) 2018-2024 Varga Consulting, Toronto, ON, Canada.
 * Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada. */

#ifndef  H5CPP_AOPEN_HPP
#define  H5CPP_AOPEN_HPP
namespace h5 {
	template<class HID_T>
	inline typename std::enable_if<h5::impl::is_valid_attr<HID_T>::value,
    h5::at_t>::type open(const  HID_T& parent, const std::string& path, const h5::acpl_t& acpl = h5::default_acpl ){

		H5CPP_CHECK_PROP( acpl, h5::error::io::attribute::open, "invalid attribute creation property" );
		hid_t attr = H5I_UNINIT;
	   	H5CPP_CHECK_NZ((
			attr = H5Aopen( static_cast<hid_t>(parent), path.c_str(), static_cast<hid_t>(acpl))),
				   							h5::error::io::attribute::open, "can't open attribute..." );
     	return  h5::at_t{attr};
    }
}
#endif

