// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
//
// Bench: h5::high_throughput direct-chunk vs standard HDF5 path
// NOTE: h5::high_throughput activation is currently broken on HDF5 1.10.x.
//       This bench will be enabled once the activation path is fixed.

#define ANKERL_NANOBENCH_IMPLEMENT
#include "../harness/nanobench.h"
#include "../harness/config.hpp"
#include "../fixtures/synthetic.hpp"
#include <h5cpp/all>
#include <vector>
#include <cstdio>

static const char* const kFile = "/dev/shm/h5cpp_bench_ht.h5";
static void cleanup() { std::remove(kFile); }

int main() {
    bench::fixture::Synthetic gen(42);
    const std::size_t n = bench::default_payload_size();
    auto data = gen.doubles(n);
    const std::size_t bytes = n * sizeof(double);

    // Standard path
    ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
        "high_throughput/standard/write/1m", [&] {
            cleanup();
            h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
            h5::write(fd, "data", data);
        });

    // high_throughput path (when activation works)
    // ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
    //     "high_throughput/direct/write/1m", [&] {
    //         cleanup();
    //         h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
    //         h5::ds_t ds = h5::create<double>(fd, "data",
    //             h5::current_dims{n}, h5::max_dims{n}, h5::chunk{1024});
    //         h5::write(ds, data, h5::high_throughput{1024});
    //     });

    cleanup();
}
