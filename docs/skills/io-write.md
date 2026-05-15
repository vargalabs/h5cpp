---
skill: io-write
tags: [io, write, stl, vector, partial, offset, stride, string]
hdf5_min: 1.10
tested: true
source: test/H5Dio_roundtrip.cpp, examples/container/container.cpp, examples/mpi/independent.cpp
---

# Writing datasets

## One-shot write — creates dataset and writes in one call

```cpp
#include <h5cpp/all>

auto fd = h5::create("data.h5", H5F_ACC_TRUNC);

std::vector<double> v = {1.0, 2.0, 3.0};
h5::write(fd, "sensor/temperature", v);   // dataset created from v's size

std::vector<std::string> labels = {"alpha", "beta", "gamma"};
h5::write(fd, "labels", labels);
```

## Write to an existing dataset handle

```cpp
auto ds = h5::open(fd, "sensor/temperature");
h5::write(ds, v);
```

## Partial write — offset and stride into an existing dataset

```cpp
// Dataset must already exist with sufficient dimensions
arma::mat M(2, 3); M.ones();
h5::write(ds, M, h5::offset{2, 2}, h5::stride{1, 3});
```

## Write with explicit current dimensions (creates dataset)

```cpp
// Places a 1D vector as a row inside a 2D dataset
std::vector<double> v(8);
h5::write(fd, "matrix_slice", v,
    h5::current_dims{40, 50},
    h5::offset{5, 0},
    h5::count{1, 1},
    h5::stride{3, 5},
    h5::block{2, 4},
    h5::max_dims{40, H5S_UNLIMITED});
```

## Known constraints

- `h5::write(fd, "ds", 42)` — scalar write creates `H5S_SCALAR` space.
  Round-tripping via `h5::read` is blocked on #89. Use a 1-element vector.
- Strings: write as `std::vector<std::string>`, not `std::string`.
- `h5::offset`, `h5::stride`, `h5::count`, `h5::block` arguments may appear
  in any order after the data argument.
