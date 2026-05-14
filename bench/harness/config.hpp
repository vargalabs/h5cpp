// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
#pragma once

#include <vector>
#include <cstdlib>
#include <cstddef>

namespace bench {

inline bool is_ci() {
    const char* env = std::getenv("H5CPP_BENCH_CI");
    return env && env[0] == '1';
}

// Returns tiny sizes for CI (smoke test only) or 1GB/10GB for real benchmarking.
inline std::vector<std::size_t> payload_sizes() {
    if (is_ci()) {
        return {1'000};  // 8 KB — enough to exercise the code path
    }
    return {
        125'000'000,    // 1 GB
        1'250'000'000,  // 10 GB
    };
}

// Single default payload size for benches that don't sweep sizes.
inline std::size_t default_payload_size() {
    return is_ci() ? 1'000 : 125'000'000;  // 8 KB vs 1 GB
}

} // namespace bench
