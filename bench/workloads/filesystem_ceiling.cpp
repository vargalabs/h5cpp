// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
//
// Bench: raw POSIX fwrite/fread ceiling for bandwidth comparison

#define ANKERL_NANOBENCH_IMPLEMENT
#include "../harness/nanobench.h"
#include "../harness/config.hpp"
#include "../fixtures/synthetic.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

static const char* const kFile = "/dev/shm/h5cpp_bench_filesystem.raw";

static void cleanup() { std::remove(kFile); }

int main() {
    bench::fixture::Synthetic gen(42);
    const std::vector<std::size_t> sizes = bench::payload_sizes();

    for (std::size_t n : sizes) {
        auto data = gen.doubles(n);
        const std::size_t bytes = data.size() * sizeof(double);
        std::vector<double> out(n);

        // Warm-up
        cleanup();
        {
            FILE* fp = std::fopen(kFile, "wb");
            std::fwrite(data.data(), sizeof(double), n, fp);
            std::fclose(fp);
        }

        ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
            "fwrite/vector_1d/" + std::to_string(n), [&] {
                cleanup();
                FILE* fp = std::fopen(kFile, "wb");
                std::fwrite(data.data(), sizeof(double), n, fp);
                std::fflush(fp);
                std::fclose(fp);
            });

        ankerl::nanobench::Bench().unit("byte").batch(bytes).run(
            "fread/vector_1d/" + std::to_string(n), [&] {
                FILE* fp = std::fopen(kFile, "rb");
                std::fread(out.data(), sizeof(double), n, fp);
                std::fclose(fp);
                ankerl::nanobench::doNotOptimizeAway(out.data());
            });

        cleanup();
    }
}
