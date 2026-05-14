// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
#pragma once

#include <vector>
#include <random>
#include <cstddef>

namespace bench::fixture {

// Deterministic pseudo-random data generator for reproducible benchmarks.
// Same seed => same bytes => same compression ratios.
class Synthetic {
public:
    explicit Synthetic(std::uint64_t seed = 42) : rng_(seed) {}

    std::vector<double> doubles(std::size_t n) {
        std::vector<double> out(n);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (auto& v : out) v = dist(rng_);
        return out;
    }

    std::vector<std::uint8_t> bytes(std::size_t n) {
        std::vector<std::uint8_t> out(n);
        std::uniform_int_distribution<int> dist(0, 255);
        for (auto& v : out) v = static_cast<std::uint8_t>(dist(rng_));
        return out;
    }

private:
    std::mt19937_64 rng_;
};

} // namespace bench::fixture
