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
    1.10) version="${H5CPP_HDF5_1_10_VERSION}"; tag="${H5CPP_HDF5_1_10_TAG}" ;;
    1.12) version="${H5CPP_HDF5_1_12_VERSION}"; tag="${H5CPP_HDF5_1_12_TAG}" ;;
    1.14) version="${H5CPP_HDF5_1_14_VERSION}"; tag="${H5CPP_HDF5_1_14_TAG}" ;;
    2.x)  version="${H5CPP_HDF5_2_X_VERSION}"; tag="${H5CPP_HDF5_2_X_TAG}" ;;
    *) echo "unknown HDF5 milestone: ${milestone}" >&2; exit 2 ;;
esac

case "${flavor}" in
    serial) parallel="OFF" ;;
    mpi)    parallel="ON" ;;
    *) echo "unknown HDF5 flavor: ${flavor}" >&2; exit 2 ;;
esac

image="h5cpp-hdf5:${version}-${flavor}"

echo "=== Building ${image} with ${engine} ==="
"${engine}" build \
    -f "${script_dir}/Containerfile" \
    -t "${image}" \
    --build-arg "HDF5_VERSION=${version}" \
    --build-arg "HDF5_TAG=${tag}" \
    --build-arg "HDF5_FLAVOR=${flavor}" \
    --build-arg "HDF5_PARALLEL=${parallel}" \
    "${repo_root}"

echo "${image}"
