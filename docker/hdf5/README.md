# HDF5 Container Matrix

This directory provides Podman/Docker-compatible validation images for h5cpp
against serial and MPI HDF5 milestone releases.

Default container engine is `podman`. Override with:

```bash
CONTAINER_ENGINE=docker docker/hdf5/run-one.sh 1.12 serial
```

## Milestones

The milestone mapping lives in `matrix.env`.

| Milestone | Default version | Flavor |
|-----------|-----------------|--------|
| `1.10` | `1.10.10` | `serial`, `mpi` |
| `1.12` | `1.12.2` | `serial`, `mpi` |
| `1.14` | `1.14.6` | `serial`, `mpi` |
| `2.x` | `2.1.1` | `serial`, `mpi` experimental |

## Common Commands

Build one image:

```bash
docker/hdf5/build-one.sh 1.12 serial
docker/hdf5/build-one.sh 1.12 mpi
```

Run one image:

```bash
docker/hdf5/run-one.sh 1.12 serial
docker/hdf5/run-one.sh 1.12 mpi
```

Run only selected tests:

```bash
H5CPP_TEST_REGEX='test-h5rreference|test-h5dio|test-h5dio_roundtrip' \
    docker/hdf5/run-one.sh 1.12 serial
```

Build or run the full matrix:

```bash
docker/hdf5/build-matrix.sh
docker/hdf5/run-matrix.sh
```

HDF5 `2.x` lanes are forward-compatibility probes and are allowed to fail by
default in `run-matrix.sh`. Disable that behavior with:

```bash
H5CPP_ALLOW_HDF5_2_FAILURE=OFF docker/hdf5/run-matrix.sh
```

## MPI

The `mpi` flavor builds HDF5 with `HDF5_ENABLE_PARALLEL=ON` and uses OpenMPI
inside the container. `run-one.sh <milestone> mpi` builds and runs the normal
h5cpp tests, then builds and smoke-runs the MPI examples with `mpiexec -n 2`.

Override the process count with:

```bash
H5CPP_MPI_NP=4 docker/hdf5/run-one.sh 1.12 mpi
```
