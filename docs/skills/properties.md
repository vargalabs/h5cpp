---
skill: properties
tags: [properties, dcpl, fapl, lcpl, chunk, gzip, compression, shuffle, fletcher32, fill_value]
hdf5_min: 1.10
tested: true
source: examples/basics/basics.cpp, test/H5Pall.cpp
---

# Property list composition

h5cpp property lists are composed with `operator|`. All properties in a chain
must be the same list type (dcpl, fapl, lcpl, etc.) — mixing types is a
compile-time error.

## Dataset creation properties (dcpl)

```cpp
// Chunk + compression pipeline (most common)
h5::dcpl_t dcpl = h5::chunk{2, 3}
                | h5::shuffle           // reorder bytes for better compression
                | h5::gzip{9}           // deflate level 0–9
                | h5::fletcher32;       // checksum

// Fill value
h5::dcpl_t dcpl = h5::chunk{64} | h5::fill_value<short>{42};

// All filters in one expression (passed directly to create)
auto ds = h5::create<short>(fd, "dataset",
    h5::current_dims{10, 20},
    h5::max_dims{10, H5S_UNLIMITED},
    h5::chunk{2, 3} | h5::gzip{6} | h5::shuffle);
```

### Available dcpl properties

| Property | Effect |
|---|---|
| `h5::chunk{n}` / `h5::chunk{n,m}` | enable chunked layout, set chunk shape |
| `h5::gzip{level}` | deflate compression, level 0–9 |
| `h5::shuffle` | byte-shuffle filter (improves gzip ratio for numerics) |
| `h5::nbit` | n-bit filter for integer types |
| `h5::fletcher32` | Fletcher-32 checksum |
| `h5::fill_value<T>{v}` | default fill value for unwritten elements |
| `h5::alloc_time{H5D_ALLOC_TIME_EARLY}` | allocation timing |
| `h5::fill_time{H5D_FILL_TIME_ALLOC}` | fill timing |
| `h5::layout{H5D_CHUNKED}` | explicit layout (usually inferred from chunk) |

## File access properties (fapl)

```cpp
h5::fapl_t fapl = h5::fclose_degree_weak | h5::stdio;
```

## File creation properties (fcpl)

```cpp
h5::fcpl_t fcpl = h5::file_space_page_size{4096} | h5::userblock{512};
auto fd = h5::create("data.h5", H5F_ACC_TRUNC, fcpl);
```

## Link creation properties (lcpl)

```cpp
// Auto-create intermediate groups + UTF-8 link names
h5::lcpl_t lcpl = h5::create_path | h5::utf8;

// Pass lcpl to create
auto ds = h5::create<float>(fd, "/a/b/c/dataset",
    h5::current_dims{100}, lcpl);
```

`h5::create_path` and `h5::utf8` are pre-built `lcpl_t` instances — they are
or'd directly without constructing a named object:

```cpp
auto ds = h5::create<float>(fd, "/a/b/c/ds",
    h5::current_dims{100},
    h5::create_path | h5::utf8,   // lcpl inline
    h5::chunk{32} | h5::gzip{6}); // dcpl inline
```

## Operator|= (in-place composition)

```cpp
h5::dcpl_t dcpl = h5::chunk{16} | h5::gzip{4};
dcpl |= h5::fletcher32;   // add to existing list
```

## Default property lists

| Symbol | Meaning |
|---|---|
| `h5::default_lcpl` | `create_path \| utf8` |
| `h5::default_dapl` | `H5P_DEFAULT` dataset access |
| `h5::default_fapl` | `H5P_DEFAULT` file access |

Pass these explicitly when you need a placeholder:

```cpp
auto ds = h5::create<short>(fd, "ds",
    h5::current_dims{10},
    h5::default_lcpl,
    h5::chunk{4} | h5::gzip{3},
    h5::default_dapl);
```
