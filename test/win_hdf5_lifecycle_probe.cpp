// SPDX-License-Identifier: MIT
// This file is part of H5CPP.

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>
#include <hdf5.h>

#include <h5cpp/all>

#ifndef H5CPP_WIN_PROBE_MODE
#define H5CPP_WIN_PROBE_MODE 0
#endif

namespace {
    void probe(const char* label) {
        std::fprintf(stderr, "[WINLIFE] %s\n", label);
        std::fflush(stderr);
    }

    std::string filename(const char* suffix) {
        return std::string("win_lifecycle_") + suffix + ".h5";
    }

    void remove_file(const std::string& path) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    int raw_minimal() {
        auto path = filename("raw_minimal");
        remove_file(path);
        probe("raw_minimal: H5Fcreate");
        hid_t fd = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hsize_t dims[1] = {3};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t ds = H5Dcreate2(fd, "data", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        double data[3] = {1.0, 2.0, 3.0};
        probe("raw_minimal: H5Dwrite");
        herr_t err = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(ds);
        H5Sclose(space);
        H5Fclose(fd);
        remove_file(path);
        probe("raw_minimal: return");
        return err < 0 ? 1 : 0;
    }

    int raw_property_lists() {
        auto path = filename("raw_property_lists");
        remove_file(path);
        hsize_t dims[1] = {3};
        hsize_t chunk[1] = {3};
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_create_intermediate_group(lcpl, 1);
        H5Pset_char_encoding(lcpl, H5T_CSET_UTF8);
        H5Pset_chunk(dcpl, 1, chunk);
        H5Pset_deflate(dcpl, 1);
        hid_t fd = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t space = H5Screate_simple(1, dims, nullptr);
        probe("raw_property_lists: H5Dcreate2");
        hid_t ds = H5Dcreate2(fd, "nested/data", H5T_NATIVE_DOUBLE, space, lcpl, dcpl, dapl);
        double data[3] = {1.0, 2.0, 3.0};
        probe("raw_property_lists: H5Dwrite");
        herr_t err = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl, data);
        H5Dclose(ds);
        H5Sclose(space);
        H5Pclose(dxpl);
        H5Pclose(dapl);
        H5Pclose(dcpl);
        H5Pclose(lcpl);
        H5Fclose(fd);
        remove_file(path);
        probe("raw_property_lists: return");
        return err < 0 ? 1 : 0;
    }

    int h5cpp_static_defaults() {
        probe("h5cpp_static_defaults: inspect default handles");
        std::fprintf(stderr, "[WINLIFE] default_lcpl=%lld valid=%d\n",
            (long long)static_cast<hid_t>(h5::default_lcpl), (int)H5Iis_valid(h5::default_lcpl));
        std::fprintf(stderr, "[WINLIFE] default_dapl=%lld valid=%d\n",
            (long long)static_cast<hid_t>(h5::default_dapl), (int)H5Iis_valid(h5::default_dapl));
        std::fprintf(stderr, "[WINLIFE] default_dxpl=%lld valid=%d\n",
            (long long)static_cast<hid_t>(h5::default_dxpl), (int)H5Iis_valid(h5::default_dxpl));
        {
            h5::lcpl_t lcpl = h5::create_path | h5::utf8;
            h5::dcpl_t dcpl = h5::chunk{3} | h5::gzip{1};
            h5::dapl_t dapl{H5Pcreate(H5P_DATASET_ACCESS)};
            std::fprintf(stderr, "[WINLIFE] local_lcpl=%lld local_dcpl=%lld local_dapl=%lld\n",
                (long long)static_cast<hid_t>(lcpl), (long long)static_cast<hid_t>(dcpl),
                (long long)static_cast<hid_t>(dapl));
        }
        probe("h5cpp_static_defaults: local property scope closed");
        return 0;
    }

    int h5cpp_create_only() {
        auto path = filename("h5cpp_create_only");
        remove_file(path);
        {
            probe("h5cpp_create_only: create file");
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            probe("h5cpp_create_only: create dataset");
            h5::ds_t ds = h5::create<double>(fd, "nested/data", h5::current_dims{3});
            std::fprintf(stderr, "[WINLIFE] ds=%lld\n", (long long)static_cast<hid_t>(ds));
        }
        remove_file(path);
        probe("h5cpp_create_only: return");
        return 0;
    }

    int h5cpp_high_level_write() {
        auto path = filename("h5cpp_high_level_write");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("h5cpp_high_level_write: before h5::write(fd,path,vector)");
            h5::write(fd, "vec", data);
            probe("h5cpp_high_level_write: after h5::write(fd,path,vector)");
        }
        remove_file(path);
        probe("h5cpp_high_level_write: return");
        return 0;
    }

    int h5cpp_high_level_write_exit() {
        auto path = filename("h5cpp_high_level_write_exit");
        remove_file(path);
        h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
        std::vector<double> data{1.0, 2.0, 3.0};
        probe("h5cpp_high_level_write_exit: before h5::write(fd,path,vector)");
        h5::write(fd, "vec", data);
        probe("h5cpp_high_level_write_exit: after h5::write, calling _Exit");
        std::fflush(stderr);
        std::_Exit(0);
    }

    int h5cpp_high_level_read() {
        auto path = filename("h5cpp_high_level_read");
        remove_file(path);
        {
            hid_t fd = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            hsize_t dims[1] = {3};
            hid_t space = H5Screate_simple(1, dims, nullptr);
            hid_t ds = H5Dcreate2(fd, "vec", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            double data[3] = {1.0, 2.0, 3.0};
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            H5Dclose(ds);
            H5Sclose(space);
            H5Fclose(fd);
        }
        {
            h5::fd_t fd = h5::open(path, H5F_ACC_RDONLY);
            probe("h5cpp_high_level_read: before h5::read<vector<double>>");
            auto data = h5::read<std::vector<double>>(fd, "vec");
            probe("h5cpp_high_level_read: after h5::read<vector<double>>");
            std::fprintf(stderr, "[WINLIFE] read_size=%zu\n", data.size());
        }
        remove_file(path);
        probe("h5cpp_high_level_read: return");
        return 0;
    }

    int h5cpp_attribute_write() {
        auto path = filename("h5cpp_attribute_write");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::ds_t ds = h5::create<double>(fd, "data", h5::current_dims{3});
            probe("h5cpp_attribute_write: before scalar attribute");
            ds["scalar"] = 42;
            probe("h5cpp_attribute_write: before vector attribute");
            ds["vector"] = std::vector<double>{1.0, 2.0, 3.0};
            probe("h5cpp_attribute_write: before string attribute");
            ds["string"] = std::string("value");
        }
        remove_file(path);
        probe("h5cpp_attribute_write: return");
        return 0;
    }
}

int main() {
#if H5CPP_WIN_PROBE_MODE == 1
    return raw_minimal();
#elif H5CPP_WIN_PROBE_MODE == 2
    return raw_property_lists();
#elif H5CPP_WIN_PROBE_MODE == 3
    return h5cpp_static_defaults();
#elif H5CPP_WIN_PROBE_MODE == 4
    return h5cpp_create_only();
#elif H5CPP_WIN_PROBE_MODE == 5
    return h5cpp_high_level_write();
#elif H5CPP_WIN_PROBE_MODE == 6
    return h5cpp_high_level_write_exit();
#elif H5CPP_WIN_PROBE_MODE == 7
    return h5cpp_high_level_read();
#elif H5CPP_WIN_PROBE_MODE == 8
    return h5cpp_attribute_write();
#else
    probe("unknown H5CPP_WIN_PROBE_MODE");
    return 2;
#endif
}
