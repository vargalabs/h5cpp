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
        if [[ "${milestone}" == "2.x" && "${H5CPP_ALLOW_HDF5_2_FAILURE:-ON}" == "ON" ]]; then
            "${script_dir}/run-one.sh" "${milestone}" "${flavor}" || true
        else
            "${script_dir}/run-one.sh" "${milestone}" "${flavor}"
        fi
    done
done
