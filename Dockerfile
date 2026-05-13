# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# OpenMPI permits running as root inside Docker
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV OMPI_MCA_rmaps_base_oversubscribe=1

# Speed up package installation inside Docker by disabling fsync
RUN apt-get update && \
    apt-get install -y --no-install-recommends eatmydata && \
    eatmydata apt-get install -y --no-install-recommends \
        cmake \
        g++ \
        make \
        libopenmpi-dev \
        openmpi-bin \
        libhdf5-openmpi-dev \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# CMake's FindHDF5 looks for h5cc by default; the parallel package provides h5pcc
RUN ln -sf /usr/bin/h5pcc /usr/bin/h5cc

WORKDIR /h5cpp

# Copy only what's needed: h5cpp headers, MPI examples, and the standalone CMakeLists.txt
COPY h5cpp/ ./h5cpp/
COPY examples/mpi/ ./examples/mpi/
COPY docker/mpi-CMakeLists.txt ./examples/mpi/CMakeLists.txt

# Configure and build only the MPI examples
RUN cmake -B build ./examples/mpi -DCMAKE_BUILD_TYPE=Release
RUN cmake --build build --parallel

# Default: run all three MPI examples
CMD ["bash", "-c", "\
echo '=== MPI Collective ===' && \
mpiexec -n 2 ./build/collective && \
echo '' && \
echo '=== MPI Independent ===' && \
mpiexec -n 2 ./build/independent && \
echo '' && \
echo '=== MPI Throughput ===' && \
mpiexec -n 2 ./build/throughput \
"]
