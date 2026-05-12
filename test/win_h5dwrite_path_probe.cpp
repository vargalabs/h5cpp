// SPDX-License-Identifier: MIT
// This file is part of H5CPP.

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>
#include <hdf5.h>

#include <h5cpp/all>

#ifndef H5CPP_WIN_WRITE_PROBE_MODE
#define H5CPP_WIN_WRITE_PROBE_MODE 0
#endif

namespace {
    void probe(const char* label) {
        std::fprintf(stderr, "[WINWRITE] %s\n", label);
        std::fflush(stderr);
    }

    std::string filename(const char* suffix) {
        return std::string("win_write_path_") + suffix + ".h5";
    }

    void remove_file(const std::string& path) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    h5::ds_t make_dataset(const h5::fd_t& fd, const char* path) {
        probe("make_dataset: h5::create<double>");
        return h5::create<double>(fd, path, h5::current_dims{3});
    }

    int ds_ptr_count() {
        auto path = filename("ds_ptr_count");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::ds_t ds = make_dataset(fd, "vec");
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("ds_ptr_count: before h5::write(ds, ptr, count)");
            h5::write(ds, data.data(), h5::count{data.size()});
            probe("ds_ptr_count: after h5::write(ds, ptr, count)");
        }
        remove_file(path);
        return 0;
    }

    int ds_object() {
        auto path = filename("ds_object");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::ds_t ds = make_dataset(fd, "vec");
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("ds_object: before h5::write(ds, vector)");
            h5::write(ds, data);
            probe("ds_object: after h5::write(ds, vector)");
        }
        remove_file(path);
        return 0;
    }

    int fd_ptr_count() {
        auto path = filename("fd_ptr_count");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("fd_ptr_count: before h5::write(fd, path, ptr, count)");
            h5::write(fd, "vec", data.data(), h5::count{data.size()});
            probe("fd_ptr_count: after h5::write(fd, path, ptr, count)");
        }
        remove_file(path);
        return 0;
    }

    int fd_object_explicit_dims() {
        auto path = filename("fd_object_explicit_dims");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("fd_object_explicit_dims: before h5::write(fd, path, vector, current_dims)");
            h5::write(fd, "vec", data, h5::current_dims{3});
            probe("fd_object_explicit_dims: after h5::write(fd, path, vector, current_dims)");
        }
        remove_file(path);
        return 0;
    }

    int fd_object_explicit_plists() {
        auto path = filename("fd_object_explicit_plists");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::lcpl_t lcpl = h5::create_path | h5::utf8;
            h5::dcpl_t dcpl{H5Pcreate(H5P_DATASET_CREATE)};
            h5::dapl_t dapl{H5Pcreate(H5P_DATASET_ACCESS)};
            h5::dxpl_t dxpl{H5Pcreate(H5P_DATASET_XFER)};
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("fd_object_explicit_plists: before h5::write(fd, path, vector, plists)");
            h5::write(fd, "vec", data, h5::current_dims{3}, lcpl, dcpl, dapl, dxpl);
            probe("fd_object_explicit_plists: after h5::write(fd, path, vector, plists)");
        }
        remove_file(path);
        return 0;
    }

    int fd_object_existing_dataset() {
        auto path = filename("fd_object_existing_dataset");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            make_dataset(fd, "vec");
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("fd_object_existing_dataset: before h5::write(fd, existing path, vector)");
            h5::write(fd, "vec", data);
            probe("fd_object_existing_dataset: after h5::write(fd, existing path, vector)");
        }
        remove_file(path);
        return 0;
    }

    int low_level_explicit_spaces() {
        auto path = filename("low_level_explicit_spaces");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::ds_t ds = make_dataset(fd, "vec");
            h5::sp_t mem_space = h5::create_simple(static_cast<hsize_t>(3));
            h5::sp_t file_space{H5Dget_space(static_cast<hid_t>(ds))};
            h5::dxpl_t dxpl{H5P_DEFAULT};
            double data[3] = {1.0, 2.0, 3.0};
            const double* ptr = data;
            void (*write_fn)(const h5::ds_t&, const h5::sp_t&, const h5::sp_t&,
                const h5::dxpl_t&, const double*) = h5::write<double>;
            h5::select_all(mem_space);
            h5::select_all(file_space);
            probe("low_level_explicit_spaces: before h5::write(ds, spaces, dxpl, ptr)");
            write_fn(ds, mem_space, file_space, dxpl, ptr);
            probe("low_level_explicit_spaces: after h5::write(ds, spaces, dxpl, ptr)");
        }
        remove_file(path);
        return 0;
    }

    int fd_object_high_throughput_dapl() {
        auto path = filename("fd_object_high_throughput_dapl");
        remove_file(path);
        {
            h5::fd_t fd = h5::create(path, H5F_ACC_TRUNC);
            h5::dapl_t dapl = h5::high_throughput;
            std::vector<double> data{1.0, 2.0, 3.0};
            probe("fd_object_high_throughput_dapl: before h5::write(fd, path, vector, high_throughput)");
            h5::write(fd, "vec", data, h5::current_dims{3}, dapl);
            probe("fd_object_high_throughput_dapl: after h5::write(fd, path, vector, high_throughput)");
        }
        remove_file(path);
        return 0;
    }
}

int main() {
#if H5CPP_WIN_WRITE_PROBE_MODE == 1
    return ds_ptr_count();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 2
    return ds_object();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 3
    return fd_ptr_count();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 4
    return fd_object_explicit_dims();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 5
    return fd_object_explicit_plists();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 6
    return fd_object_existing_dataset();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 7
    return low_level_explicit_spaces();
#elif H5CPP_WIN_WRITE_PROBE_MODE == 8
    return fd_object_high_throughput_dapl();
#else
    probe("unknown H5CPP_WIN_WRITE_PROBE_MODE");
    return 2;
#endif
}
