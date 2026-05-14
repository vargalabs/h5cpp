/*
 * Copyright (c) 2018-2021 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once
#include <cstring>

namespace h5::impl::reference {
    inline herr_t create_region(h5::reference_t& ref, ::hid_t loc, const char* path, ::hid_t space, ::hid_t oapl=static_cast<::hid_t>(0)) {
#if H5_VERSION_GE(1,12,0)
        return H5Rcreate_region(loc, path, space, oapl, &ref.value);
#else
        (void)oapl;
        return H5Rcreate(&ref.value, loc, path, H5R_DATASET_REGION, space);
#endif
    }

    inline h5::ds_t open_object(const h5::ds_t& ds, h5::reference_t& ref, ::hid_t oapl=static_cast<::hid_t>(0)) {
#if H5_VERSION_GE(1,12,0)
        (void)ds;
        return h5::ds_t(H5Ropen_object(&ref.value, H5P_DEFAULT, oapl));
#else
        return h5::ds_t(H5Rdereference2(ds, oapl, H5R_DATASET_REGION, &ref.value));
#endif
    }

    inline h5::sp_t open_region(const h5::ds_t& ds, h5::reference_t& ref, ::hid_t oapl=static_cast<::hid_t>(0)) {
#if H5_VERSION_GE(1,12,0)
        (void)ds;
        return h5::sp_t(H5Ropen_region(&ref.value, H5P_DEFAULT, oapl));
#else
        h5::ds_t target = open_object(ds, ref, oapl);
        return h5::sp_t(H5Rget_region(target, H5R_DATASET_REGION, &ref.value));
#endif
    }

    inline herr_t copy(h5::reference_t& dst, const h5::reference_t& src) {
#if H5_VERSION_GE(1,12,0)
        return H5Rcopy(&src.value, &dst.value);
#else
        std::memcpy(&dst.value, &src.value, sizeof(dst.value));
        return 0;
#endif
    }

    inline htri_t equal(const h5::reference_t& lhs, const h5::reference_t& rhs) {
#if H5_VERSION_GE(1,12,0)
        return H5Requal(&lhs.value, &rhs.value);
#else
        return std::memcmp(&lhs.value, &rhs.value, sizeof(lhs.value)) == 0 ? 1 : 0;
#endif
    }

    inline herr_t destroy(h5::reference_t& ref) {
#if H5_VERSION_GE(1,12,0)
        return H5Rdestroy(&ref.value);
#else
        (void)ref;
        return 0;
#endif
    }

    inline herr_t reclaim(::hid_t type, ::hid_t space, ::hid_t plist, void* buf) {
#if H5_VERSION_GE(1,12,0)
        return H5Treclaim(type, space, plist, buf);
#else
        return H5Dvlen_reclaim(type, space, plist, buf);
#endif
    }
}
