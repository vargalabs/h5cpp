---
hide:
  - toc
---

# Getting Started

H5CPP is a modern C++ interface for HDF5 designed to make persistence of structured and numerical data feel natural in C++. It is particularly useful when working with:

- :material-file-code-outline:{.icon} structured records such as trades, quotes, fills, and other event data
- :material-grid:{.icon} dense numerical datasets such as vectors, matrices, tensors, and feature tables
- :material-lan-connect:{.icon} workflows that begin on a laptop and later move to MPI-enabled environments

This page gets you from zero to a working example with the least possible drama.

## Prerequisites

Before building examples, make sure you have:

- :material-linux:{.icon} a Linux system or another Unix-like environment
- :material-language-cpp:{.icon} a modern C++ compiler
- :material-database:{.icon} the HDF5 development libraries installed
- :material-hammer-wrench:{.icon} `make` or `cmake`

Recommended compiler versions:

- **GCC** 11+
- **Clang** 13+

H5CPP is header-only in normal use, so once dependencies are available, integration is straightforward.

## Install HDF5

On Ubuntu or Debian:

```bash
sudo apt update
sudo apt install libhdf5-dev
```

For MPI-enabled parallel HDF5:

```bash
sudo apt install libhdf5-mpi-dev
```

You may also want common numerical backends depending on which examples you plan to use, such as Eigen, Armadillo, or Boost.

## Get the source

Clone the repository:

```bash
git clone https://github.com/steven-varga/h5cpp.git
cd h5cpp
```

If you are working from a local checkout already, then you are ahead of the tutorial. A rare and suspiciously efficient event.

## Build the examples

The repository contains a broad collection of examples under `examples/`.

To build example code with the provided Makefiles:

```bash
cd examples/basics
make
./basics
```

For many users, this is the fastest way to verify that:

* the compiler is working
* HDF5 headers and libraries are visible
* H5CPP is usable in the current environment

## Your first program

Here is a minimal example that writes and reads a dense matrix.

```cpp
#include <armadillo>
#include <h5cpp/all>

int main() {
    arma::mat m = arma::randu<arma::mat>(4, 3);

    h5::fd_t fd = h5::create("example.h5", H5F_ACC_TRUNC);
    h5::write(fd, "/matrix", m);

    arma::mat restored = h5::read<arma::mat>(fd, "/matrix");
}
```

This program does three things:

1. creates a random matrix
2. writes it into an HDF5 file
3. reads it back into C++

That is the core H5CPP workflow in one screenful.

## Compile it

If Armadillo and HDF5 are installed system-wide, a typical compile line looks like:

```bash
g++ -std=c++17 first_example.cpp -o first_example -lhdf5 -larmadillo
```

Depending on your HDF5 installation, you may need include and library paths, for example:

```bash
g++ -std=c++17 first_example.cpp -I/usr/include/hdf5/serial -lhdf5_serial -larmadillo
```

The exact flags vary a bit by distribution. Linux packaging: eternally creative, never boring.

## Structured data example

H5CPP is not limited to matrices. It is also very useful for structured record data.

```cpp
#include <h5cpp/core>
#include <h5cpp/io>

struct quote_t {
    uint64_t ts;
    double bid;
    double ask;
    uint32_t bid_size;
    uint32_t ask_size;
};

int main() {
    auto fd = h5::create("quotes.h5", H5F_ACC_TRUNC);

    h5::pt_t pt = h5::create<quote_t>(
        fd,
        "/quotes",
        h5::max_dims{H5S_UNLIMITED},
        h5::chunk{1024}
    );

    quote_t q{1700000000, 101.25, 101.27, 400, 250};
    h5::append(pt, q);
}
```

This is the pattern you would use for appendable event streams such as:

* market trades
* quotes
* fills
* order updates
* telemetry records

That is where H5CPP becomes especially attractive for trading and research systems.

## Where to go next

Once the first example works, the next useful stops are:

* :material-play-box-outline:{.icon} **Examples** — browse practical examples in `examples/`
* :material-table-large:{.icon} **Linear Algebra** — see supported matrix backends
* :material-file-tree:{.icon} **Compound Types** — learn how structs are mapped
* :material-lan-connect:{.icon} **Parallel HDF5** — move into MPI-enabled workflows
* :material-book-open-page-variant:{.icon} **API Reference** — inspect the available interfaces

## Suggested reading path

A sensible order is:

1. `examples/basics`
2. `examples/stl`
3. `examples/compound`
4. `examples/linalg`
5. `examples/packet-table`
6. `examples/mpi`

That sequence moves from basic persistence toward structured and high-throughput workflows.

## Common problems
### HDF5 headers not found

Install the development package:

```bash
sudo apt install libhdf5-dev
```

### Linker cannot find HDF5

Some distributions package HDF5 under `serial` paths or use names like `-lhdf5_serial`. Check:

```bash
dpkg -L libhdf5-dev | grep include
dpkg -L libhdf5-dev | grep '\.so'
```

### MPI example does not build

Install the MPI-enabled HDF5 package:

```bash
sudo apt install libhdf5-mpi-dev
```

and use the appropriate MPI compiler wrapper where needed:

```bash
mpicxx --version
```

## Summary

H5CPP gives you a concise C++ interface for HDF5 that works well for both:

* :material-file-code-outline:{.icon} structured record data
* :material-grid:{.icon} homogeneous numerical datasets

It is a good fit for scientific computing, quantitative research, and trading infrastructure where the same data often needs to move cleanly between C++, Python, Julia, MATLAB, and HPC workflows.
