---
skill: io-read
tags: [io, read, stl, vector, partial, offset, count, string]
hdf5_min: 1.10
tested: true
source: test/H5Dio_roundtrip.cpp, examples/mpi/independent.cpp
---

# Reading datasets

## Full read — returns owned object

```cpp
#include <h5cpp/all>

auto fd = h5::open("data.h5", H5F_ACC_RDONLY);

auto v      = h5::read<std::vector<double>>(fd, "sensor/temperature");
auto labels = h5::read<std::vector<std::string>>(fd, "labels");
auto M      = h5::read<arma::mat>(fd, "results/matrix");
```

The template argument is the target C++ type. h5cpp handles all allocation.

## Read into pre-allocated buffer (high-performance loops)

```cpp
// Open dataset once, outside the loop
auto ds = h5::open(fd, "sensor/temperature");
std::vector<double> buf(1024);

for (...) {
    h5::read(ds, buf.data(), h5::offset{i * 1024}, h5::count{1024});
}
```

Pass a raw pointer + `h5::count` when the buffer is pre-allocated.
This avoids allocation on every iteration.

## Partial read — hyperslab selection

```cpp
auto ds  = h5::open(fd, "results/matrix");
auto row = h5::read<arma::rowvec>(ds, h5::offset{3, 0}, h5::count{1, 20});
```

## Read via file path (without opening fd separately)

```cpp
auto M = h5::read<arma::mat>("data.h5", "results/matrix");
```

## Type constraints

| Want | Use |
|---|---|
| `std::string` (single) | `h5::read<std::vector<std::string>>` — returns `vector` of size 1 |
| scalar `int` | `h5::read<std::vector<int>>` — round-trip scalars blocked on #89 |
| matrix | `h5::read<arma::mat>` / `h5::read<Eigen::MatrixXd>` etc. |

**`h5::read<std::string>` does not compile.** The decay path routes
`std::string` into the generic read (not the string specialization).
Always use `std::vector<std::string>`.

## Error on missing dataset

```cpp
// Throws h5::error::io::dataset::open
h5::read<std::vector<int>>(fd, "nonexistent");
```

Use `h5::mute()` / `h5::unmute()` around expected failures to suppress
HDF5's CAPI error output.
