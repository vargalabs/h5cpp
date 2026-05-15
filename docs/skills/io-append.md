---
skill: io-append
tags: [io, append, streaming, pt, packet-table, chunk, flush]
hdf5_min: 1.10
tested: true
source: test/H5Dappend.cpp
---

# Streaming append / packet table

`h5::pt_t` wraps a chunked extendable dataset for efficient streaming append.
Elements accumulate in an internal buffer; the buffer flushes to disk when full
or when `h5::flush(pt)` is called. The destructor flushes automatically.

## Setup — create the backing dataset

```cpp
#include <h5cpp/all>

auto fd = h5::create("stream.h5", H5F_ACC_TRUNC);

// current_dims{0} = starts empty; max_dims{UNLIMITED} = grows without bound
auto ds = h5::create<int>(fd, "events",
    h5::current_dims{0},
    h5::max_dims{H5S_UNLIMITED},
    h5::chunk{1024});          // chunk size = flush granularity

h5::pt_t pt(ds);
```

## Append scalars

```cpp
for (int i = 0; i < 25; ++i)
    h5::append(pt, i);

h5::flush(pt);   // write remaining buffer to disk
```

## Append a vector chunk

```cpp
std::vector<int> chunk = {1, 2, 3, 4, 5};
h5::append(pt, chunk);
h5::flush(pt);
```

## Append strings

```cpp
auto ds = h5::create<std::string>(fd, "labels",
    h5::current_dims{0},
    h5::max_dims{H5S_UNLIMITED},
    h5::chunk{128});
h5::pt_t pt(ds);

h5::append(pt, std::string{"alpha"});
h5::append(pt, "beta");    // const char* also accepted
h5::flush(pt);
```

## Flush on destruction (RAII)

```cpp
{
    h5::pt_t pt(ds);
    for (int i = 0; i < 7; ++i)
        h5::append(pt, i);
    // destructor calls flush() — partial chunk written with fill value
}
```

## 2D streaming (fixed inner extent)

```cpp
// Shape: rows grow, columns fixed at 5
auto ds = h5::create<int>(fd, "matrix_stream",
    h5::current_dims{0, 5},
    h5::max_dims{H5S_UNLIMITED, 5},
    h5::chunk{64, 5});
h5::pt_t pt(ds);

std::vector<int> row = {1, 2, 3, 4, 5};
h5::append(pt, row);     // appends one row
h5::flush(pt);
```

## Key invariants

- Chunk size sets the write granularity — too small → many small writes;
  too large → large unflushed buffer. Tune to expected event rate × element size.
- A partial final chunk is written with the dataset's fill value on flush.
- `h5::pt_t` is copyable (shares the underlying dataset handle).
- Do not mix `h5::append` and `h5::write` on the same dataset handle.
