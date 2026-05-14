#!/usr/bin/env bash
# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

set -euo pipefail

src_dir="${H5CPP_SOURCE_DIR:-/src}"
build_dir="${H5CPP_BUILD_DIR:-/tmp/h5cpp-build}"
test_regex="${H5CPP_TEST_REGEX:-}"
build_examples="${H5CPP_BUILD_EXAMPLES:-OFF}"
build_type="${H5CPP_BUILD_TYPE:-Release}"

rm -rf "${build_dir}"

echo "=== HDF5 configuration ==="
if command -v h5pcc >/dev/null 2>&1 && h5pcc -showconfig >/dev/null 2>&1; then
    h5pcc -showconfig | sed -n '1,45p'
else
    h5cc -showconfig | sed -n '1,45p'
fi

echo "=== Configure h5cpp ==="
cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DCMAKE_PREFIX_PATH="${HDF5_ROOT}" \
    -DHDF5_DIR="${HDF5_ROOT}/cmake" \
    -DHDF5_ROOT="${HDF5_ROOT}" \
    -DH5CPP_BUILD_TESTS=ON \
    -DH5CPP_BUILD_EXAMPLES="${build_examples}"

echo "=== Build h5cpp ==="
cmake --build "${build_dir}" --parallel

echo "=== Run h5cpp tests ==="
ctest_dir="${build_dir}"
if [[ -f "${build_dir}/test/CTestTestfile.cmake" ]]; then
    ctest_dir="${build_dir}/test"
fi
if [[ -n "${test_regex}" ]]; then
    ctest --test-dir "${ctest_dir}" --output-on-failure -R "${test_regex}"
else
    ctest --test-dir "${ctest_dir}" --output-on-failure
fi

if [[ "${H5CPP_RUN_MPI_EXAMPLES:-OFF}" == "ON" ]]; then
    echo "=== Build MPI examples ==="
    rm -rf /tmp/h5cpp-mpi-build
    mkdir -p /tmp/h5cpp-mpi-src/examples
    cp -R "${src_dir}/h5cpp" /tmp/h5cpp-mpi-src/
    cp -R "${src_dir}/examples/mpi" /tmp/h5cpp-mpi-src/examples/
    cp /opt/h5cpp/mpi-CMakeLists.txt /tmp/h5cpp-mpi-src/examples/mpi/CMakeLists.txt
    cmake -S /tmp/h5cpp-mpi-src/examples/mpi -B /tmp/h5cpp-mpi-build \
        -DCMAKE_BUILD_TYPE="${build_type}" \
        -DCMAKE_PREFIX_PATH="${HDF5_ROOT}" \
        -DHDF5_DIR="${HDF5_ROOT}/cmake" \
        -DHDF5_ROOT="${HDF5_ROOT}"
    cmake --build /tmp/h5cpp-mpi-build --parallel
    mpiexec -n "${H5CPP_MPI_NP:-2}" /tmp/h5cpp-mpi-build/collective
    mpiexec -n "${H5CPP_MPI_NP:-2}" /tmp/h5cpp-mpi-build/independent
    mpiexec -n "${H5CPP_MPI_NP:-2}" /tmp/h5cpp-mpi-build/throughput
fi
