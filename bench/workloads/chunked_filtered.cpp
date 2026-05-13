// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
//
// Bench: chunked write with filter chains (gzip, shuffle)
// NOTE: zstd/lz4 are available in the filter engine but lack DCPL property
// wrappers in the public API. Added when H5P wrappers are implemented.

#define ANKERL_NANOBENCH_IMPLEMENT
#include "../harness/nanobench.h"
#include "../fixtures/synthetic.hpp"
#include <h5cpp/all>
#include <vector>
#include <cstdio>

static const char* const kFile = "/dev/shm/h5cpp_bench_chunked.h5";
static void cleanup() { std::remove(kFile); }

int main() {
    bench::fixture::Synthetic gen(42);
    const std::size_t n = 125'000'000;  // 1 GB
    auto data = gen.doubles(n);
    const std::size_t bytes = n * sizeof(double);

    auto run_filter = [&](const char* name, const h5::dcpl_t& dcpl) {
        ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
            std::string("chunked_write/") + name + "/1m", [&] {
                cleanup();
                h5::fd_t fd = h5::create(kFile, H5F_ACC_TRUNC);
                h5::write(fd, "data", data, dcpl);
            });
    };

    run_filter("none",       h5::chunk{1024});
    run_filter("gzip",       h5::chunk{1024} | h5::gzip{6});
    run_filter("shuffle",    h5::chunk{1024} | h5::shuffle);
    run_filter("shuffle+gz", h5::chunk{1024} | h5::shuffle | h5::gzip{6});

    cleanup();
}
