---
skill: io-create
tags: [io, create, file, dataset, raii, dcpl, chunk]
hdf5_min: 1.10
tested: true
source: test/H5Dio_roundtrip.cpp, examples/basics/basics.cpp
---

# File and dataset creation

## File create

```cpp
#include <h5cpp/all>

// Create (truncate if exists) — returns RAII fd_t, closes on scope exit
auto fd = h5::create("data.h5", H5F_ACC_TRUNC);

// Open existing
auto fd = h5::open("data.h5", H5F_ACC_RDONLY);
auto fd = h5::open("data.h5", H5F_ACC_RDWR);
```

`h5::fd_t` is move-only. `H5Fclose` is called automatically on destruction.

## File create with property lists

```cpp
auto fd = h5::create("data.h5", H5F_ACC_TRUNC,
    h5::file_space_page_size{4096} | h5::userblock{512});
```

## Dataset create — fixed dimensions

```cpp
// Type is the template parameter. Dimensions via current_dims.
auto ds = h5::create<double>(fd, "/results/matrix",
    h5::current_dims{10, 20});
```

## Dataset create — chunked, extendable, compressed

```cpp
auto ds = h5::create<short>(fd, "/sensor/raw",
    h5::current_dims{10, 20},
    h5::max_dims{10, H5S_UNLIMITED},
    h5::chunk{2, 3} | h5::gzip{9} | h5::shuffle | h5::fletcher32);
```

- `h5::max_dims{..., H5S_UNLIMITED}` marks an extendable dimension.
- Chunk must be set whenever `H5S_UNLIMITED` appears in `max_dims`.
- `h5::shuffle` before `h5::gzip` improves compression ratio for numeric data.

## Dataset create — packet table (streaming)

```cpp
// current_dims{0} + max_dims{UNLIMITED} + chunk = packet table layout
auto ds = h5::create<int>(fd, "/stream/events",
    h5::current_dims{0},
    h5::max_dims{H5S_UNLIMITED},
    h5::chunk{1024});
```

See `io-append.md` for the corresponding `h5::append` usage.

## Create with path auto-creation (lcpl)

```cpp
// Creates intermediate groups automatically
auto ds = h5::create<float>(fd, "/a/b/c/dataset",
    h5::current_dims{100},
    h5::create_path | h5::utf8);   // lcpl
```

## Handle types

| Type | Meaning |
|---|---|
| `h5::fd_t` | file descriptor — wraps `hid_t` from `H5Fcreate` / `H5Fopen` |
| `h5::ds_t` | dataset descriptor — wraps `hid_t` from `H5Dcreate` / `H5Dopen` |
| `h5::pt_t` | packet table — wraps `ds_t` for streaming append |

All are RAII — do not call `H5Fclose` / `H5Dclose` manually on managed handles.
