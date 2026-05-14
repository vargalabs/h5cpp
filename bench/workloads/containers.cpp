// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
//
// Bench: container adapter overhead (raw HDF5 C API vs h5cpp high-level)

#define ANKERL_NANOBENCH_IMPLEMENT
#include "../harness/nanobench.h"
#include "../harness/config.hpp"
#include "../fixtures/synthetic.hpp"
#include <h5cpp/all>
#include <vector>
#include <cstdio>

static const char* const kFile = "/dev/shm/h5cpp_bench_containers.h5";
static void cleanup() { std::remove(kFile); }

int main() {
    bench::fixture::Synthetic gen(42);
    const std::size_t n = bench::default_payload_size();
    auto data = gen.doubles(n);
    const std::size_t bytes = n * sizeof(double);

    // h5cpp high-level write
    ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
        "container/h5cpp/write/vector_1m", [&] {
            cleanup();
            h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
            h5::write(fd, "data", data);
        });

    // h5cpp high-level read
    {
        h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
        h5::write(fd, "data", data);
    }
    std::vector<double> out(n);
    ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
        "container/h5cpp/read/vector_1m", [&] {
            h5::read(kFile, "data", out);
            ankerl::nanobench::doNotOptimizeAway(out.data());
        });

    // Raw HDF5 C API write (baseline)
    ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
        "container/raw_c/write/vector_1m", [&] {
            cleanup();
            hid_t fid = H5Fcreate(kFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            hsize_t dims[1] = {n};
            hid_t sid = H5Screate_simple(1, dims, nullptr);
            hid_t did = H5Dcreate2(fid, "data", H5T_NATIVE_DOUBLE, sid,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
            H5Dclose(did);
            H5Sclose(sid);
            H5Fclose(fid);
        });

    // Raw HDF5 C API read (baseline)
    ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
        "container/raw_c/read/vector_1m", [&] {
            hid_t fid = H5Fopen(kFile, H5F_ACC_RDONLY, H5P_DEFAULT);
            hid_t did = H5Dopen2(fid, "data", H5P_DEFAULT);
            H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
            H5Dclose(did);
            H5Fclose(fid);
            ankerl::nanobench::doNotOptimizeAway(out.data());
        });

    cleanup();
}
