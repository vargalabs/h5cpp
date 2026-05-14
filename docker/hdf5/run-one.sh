#!/usr/bin/env bash
# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
source "${script_dir}/matrix.env"

milestone="${1:?usage: $0 <1.10|1.12|1.14|2.x> <serial|mpi>}"
flavor="${2:?usage: $0 <1.10|1.12|1.14|2.x> <serial|mpi>}"
engine="${CONTAINER_ENGINE:-podman}"

case "${milestone}" in
    1.10) version="${H5CPP_HDF5_1_10_VERSION}" ;;
    1.12) version="${H5CPP_HDF5_1_12_VERSION}" ;;
    1.14) version="${H5CPP_HDF5_1_14_VERSION}" ;;
    2.x)  version="${H5CPP_HDF5_2_X_VERSION}" ;;
    *) echo "unknown HDF5 milestone: ${milestone}" >&2; exit 2 ;;
esac

case "${flavor}" in
    serial) mpi_examples="OFF" ;;
    mpi)    mpi_examples="${H5CPP_RUN_MPI_EXAMPLES:-ON}" ;;
    *) echo "unknown HDF5 flavor: ${flavor}" >&2; exit 2 ;;
esac

image="h5cpp-hdf5:${version}-${flavor}"
if ! "${engine}" image inspect "${image}" >/dev/null 2>&1; then
    "${script_dir}/build-one.sh" "${milestone}" "${flavor}" >/dev/null
fi

echo "=== Running ${image} with ${engine} ==="
"${engine}" run --rm \
    -e "H5CPP_BUILD_EXAMPLES=${H5CPP_BUILD_EXAMPLES:-OFF}" \
    -e "H5CPP_BUILD_TYPE=${H5CPP_BUILD_TYPE:-Release}" \
    -e "H5CPP_TEST_REGEX=${H5CPP_TEST_REGEX:-}" \
    -e "H5CPP_RUN_MPI_EXAMPLES=${mpi_examples}" \
    -e "H5CPP_MPI_NP=${H5CPP_MPI_NP:-2}" \
    -v "${repo_root}:/src:Z" \
    "${image}"
