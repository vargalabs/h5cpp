---
hide:
  - toc
---

# Supported Types

H5CPP is built around a simple idea: if a C++ type has a sensible **memory layout**, **rank**, and **element type**, it should map cleanly to HDF5.

That makes the supported-type story much less mysterious than template metaprogramming sometimes suggests. In practice, H5CPP works well with two broad families of data:

- **structured records** — C/C++ POD or standard-layout structs
- **homogeneous numerical containers** — vectors, matrices, cubes, arrays, and related dense containers

This page summarizes what is supported, how H5CPP thinks about these types, and where the boundaries are.

## Mental model

At a high level, H5CPP classifies a type along a few practical axes:

- is it a **scalar**, **vector**, **matrix**, or higher-rank array?
- does it expose **contiguous memory**?
- does it have a meaningful **element type**?
- can its extents be obtained from the type or object?
- for structs, does it have a stable layout suitable for HDF5 compound mapping?

If the answer is "yes" in the right places, the type usually fits naturally.


## 1. Scalar types

These are the trivial building blocks. They map directly to native HDF5 atomic datatypes.

Supported scalar families include:

- signed and unsigned integers
- `float`, `double`
- `std::byte`
- `std::complex<T>` for supported floating-point element types

Typical examples:

```cpp
int
double
std::byte
std::complex<float>
std::complex<double>
```

Use scalar datasets for:

- configuration constants
- single-value attributes
- metadata counters
- timestamps and scalar measurements


## 2. Fixed-size C arrays

Built-in arrays work well when the shape is known at compile time.

Examples:

```cpp
int[N]
double[N][M]
float[N][M][K]
```

These are useful for:

- small dense tensors
- fixed-size records embedded inside structs
- compact in-memory layouts with no allocator overhead

The key advantage here is predictability: the type already carries the rank and extents.

## 3. Compound / record types

This is one of H5CPP's strongest areas.

Standard-layout record types can be mapped to HDF5 **compound datatypes**, which makes them ideal for:

- trades
- quotes
- fills
- market snapshots
- telemetry
- simulation events
- checkpoint metadata

Typical example:

```cpp
struct quote_t {
    uint64_t ts;
    double bid;
    double ask;
    uint32_t bid_size;
    uint32_t ask_size;
};
```

When combined with the LLVM-based H5CPP compiler-assisted reflection, even fairly rich structs collapse into a very short persistence call. That is exactly the kind of thing that should be boring, and finally becomes boring in the good way.

### Practical guidance for structs

Best results come from types that are:

- standard-layout
- trivially copyable or close to it
- free of raw ownership semantics
- composed of supported field types

A useful rule of thumb:

> if you would trust `memcpy` on the object layout, H5CPP is much more likely to like it too.


## 4. STL containers

### `std::vector<T>`

This is the workhorse dynamic container.

Supported and useful when:

- `T` is a supported scalar type
- `T` is a supported compound type
- data is stored contiguously

Examples:

```cpp
std::vector<int>
std::vector<double>
std::vector<quote_t>
```

This is a very natural fit for:

- time series
- event batches
- feature vectors
- row-major flattened tensors

### `std::vector<std::string>`

H5CPP also supports string vectors, mapping them to HDF5 variable-length string datasets.

That is handy for:

- symbols
- identifiers
- labels
- categorical metadata

### `std::array<T, N>`

Fixed-size STL arrays behave much like built-in arrays and fit well when extents are static.

Examples:

```cpp
std::array<int, N>
std::array<std::array<int, N>, M>
```

## 5. Dense linear algebra backends

H5CPP has broad support for dense numerical libraries. This is where the library becomes particularly useful for scientific computing, quantitative research, and trading analytics.

### Armadillo

Supported patterns include:

- `arma::Row<T>`
- `arma::Col<T>`
- `arma::Mat<T>`
- `arma::Cube<T>`

Use cases:

- return matrices
- factor panels
- rolling features
- covariance inputs
- signal tensors

### Eigen

Supported patterns include:

- dynamic vectors
- dynamic matrices
- dynamic arrays

Typical forms:

```cpp
Eigen::Matrix<T, Dynamic, Dynamic>
Eigen::Matrix<T, Dynamic, 1>
Eigen::Matrix<T, 1, Dynamic>
Eigen::Array<T, Dynamic, Dynamic>
Eigen::Array<T, Dynamic, 1>
Eigen::Array<T, 1, Dynamic>
```

### Blaze

Dense dynamic vectors and matrices are supported, including row/column vector distinctions and row-major / column-major matrix layouts.

### Boost uBLAS

Dense vectors and matrices are supported.

### Blitz, ITPP, dlib

These backends are also supported for the dense numerical cases represented in the examples.

The common thread across all of them is this:

- there is a well-defined element type
- rank is known or inferable
- extents can be obtained
- memory access is compatible with HDF5 transfer

## 6. Sparse and special cases

Sparse data is always trickier, because HDF5 is fundamentally happiest with dense regular layouts.

The repository contains sparse-related examples, which is a good sign, but sparse storage usually deserves a more explicit schema discussion:

- CSR / CSC as struct-of-arrays
- triplet storage
- custom compound representations

So the honest position is:

- **dense** support is first-class and natural
- **sparse** support is possible, but should be treated as a data-model design question, not magic

## 7. Strings and UTF data

String handling matters more than people admit, right until they have to persist a symbol universe or a million identifiers.

H5CPP supports:

- `std::string`
- `std::vector<std::string>`
- UTF-related examples and custom datatype paths

This makes it practical to mix:

- numerical datasets
- structured records
- human-readable metadata

inside the same HDF5 container.

## 8. What “supported” really means

For H5CPP, a supported type usually has these ingredients:

1. a sensible HDF5 datatype mapping
2. a rank that can be inferred or specified
3. stable memory access
4. extents that are known or discoverable
5. element types that themselves are supported

So “supported” is really shorthand for:

> H5CPP knows how to describe the type, size the dataspace, and move the bytes without lying to HDF5.

That is the whole game.

## Dense vs structured: the two big workflows

For most users, almost everything falls into one of these buckets.

### Structured datasets

These are record-oriented datasets:

- trades
- quotes
- fills
- order events
- risk snapshots
- telemetry
- experiment logs

Best represented as:

- standard-layout structs
- vectors of structs
- appendable packet-table-like datasets

### Homogeneous dense datasets

These are matrix- or tensor-like datasets:

- returns
- features
- factors
- covariance inputs
- price panels
- model inputs
- simulation state arrays

Best represented as:

- vectors
- matrices
- cubes
- dense arrays from supported numerical libraries

If you keep those two workflows in mind, H5CPP becomes much easier to reason about.

## Caveats

A few things deserve straight talk.

### Not every C++ type should be persisted directly

Just because a type exists does not mean it should be written as-is.

Types with:

- owning raw pointers
- virtual inheritance
- allocator-heavy internal state
- unstable layout
- opaque runtime-only invariants

are usually bad candidates for direct storage.

### HDF5 likes regularity

The closer your data looks to one of these:

- scalar
- dense array
- fixed-size record
- string
- vector of homogeneous things

the happier the mapping tends to be.

### Cross-language reuse matters

One of HDF5's biggest strengths is that the same container can later be read from:

- C++
- Python
- Julia
- MATLAB
- R

That is easiest when the stored representation is regular, explicit, and unsurprising.

## Reading the type map

The accompanying SVG organizes supported families visually, in roughly the same spirit as the original type-map sketch:

- **rows** are representative type families
- **columns** are properties H5CPP cares about
- green means **yes**
- orange-red means **no**
- gray means **not applicable**

This is not a formal proof system. It is a practical map for humans.

## Summary

H5CPP is strongest when your data is either:

- a **structured record** with stable layout, or
- a **homogeneous dense container** with a well-defined element type and shape

That covers a surprisingly large fraction of real workloads:

- scientific computing
- simulation
- machine learning pipelines
- quantitative research
- market-data storage
- backtesting inputs and outputs
- HPC batch workflows

If your type is regular, contiguous, and honest about what it contains, H5CPP will usually meet you more than halfway.