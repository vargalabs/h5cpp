// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
//
// Bench: sequential write/read throughput for 1D vector<double>

#define ANKERL_NANOBENCH_IMPLEMENT
#include "../harness/nanobench.h"
#include "../fixtures/synthetic.hpp"
#include <h5cpp/all>
#include <vector>
#include <cstdio>
#include <cstdlib>

static const char* const kFile = "/dev/shm/h5cpp_bench_sequential.h5";

static void cleanup() { std::remove(kFile); }

int main() {
    bench::fixture::Synthetic gen(42);
    const std::vector<std::size_t> sizes = {
        125'000'000,    // 1 GB
        1'250'000'000,  // 10 GB
    };

    for (std::size_t n : sizes) {
        auto data = gen.doubles(n);
        const std::size_t bytes = data.size() * sizeof(double);

        // Warm-up + ensure file exists
        cleanup();
        {
            h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
            h5::write(fd, "data", data);
        }

        ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
            "seq_write/vector_1d/" + std::to_string(n), [&] {
                cleanup();
                h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
                h5::write(fd, "data", data);
            });

        ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
            "seq_read/vector_1d/" + std::to_string(n), [&] {
                std::vector<double> out(n);
                h5::read(kFile, "data", out);
                ankerl::nanobench::doNotOptimizeAway(out.data());
            });

        cleanup();
    }
}
