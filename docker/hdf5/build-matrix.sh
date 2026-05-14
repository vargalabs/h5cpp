#!/usr/bin/env bash
# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/matrix.env"

for milestone in ${H5CPP_HDF5_MILESTONES}; do
    for flavor in ${H5CPP_HDF5_FLAVORS}; do
        "${script_dir}/build-one.sh" "${milestone}" "${flavor}"
    done
done
