#!/bin/bash
# ___________________________________________________________________________________
# SPDX-License-Identifier: MIT
# This file is part of H5CPP.
# Copyright (c) 2025-2026 Varga Labs, Toronto, ON, Canada.
# __________________________________________________________________________________

set -e

IMAGE_NAME="h5cpp-mpi"

# Build the Docker image
echo "Building ${IMAGE_NAME}..."
docker build -t "${IMAGE_NAME}" .

# Run the container (default CMD runs all three MPI examples)
echo "Running ${IMAGE_NAME}..."
docker run --rm "${IMAGE_NAME}"
