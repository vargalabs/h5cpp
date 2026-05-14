# h5cpp Bench Suite

Reproducible performance benchmarks for h5cpp I/O throughput claims.

## Quick Start

```bash
cmake -B build -DH5CPP_BUILD_BENCH=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target h5cpp-bench
./build/bench/bench-sequential
```

## Design

See `tasks/h5cpp-performance-evaluation-framework-design.md` in the workspace repo for the full architecture report.

## Workloads

| Bench | What It Measures |
|-------|------------------|
| `bench-sequential` | 1D `vector<double>` write/read throughput at 1K–10M elements |
| `bench-filesystem_ceiling` | Raw `fwrite`/`fread` bandwidth ceiling on `/dev/shm` |
| `bench-containers` | h5cpp high-level vs raw HDF5 C API overhead |
| `bench-chunked_filtered` | Chunked write with gzip, zstd, lz4, shuffle+gzip |
| `bench-high_throughput` | Standard path (direct chunk comparison deferred) |
| `bench-small_object` | Latency for single values, 4 KB vectors, structs |

## Output

Each bench prints a markdown table to stdout. To suppress the table and emit JSON:

```bash
./build/bench/bench-sequential > /dev/null 2> results.json
```

*(JSON export is a v2 enhancement; currently stdout only.)*

## Methodology

- Files are created on `/dev/shm` (tmpfs) to eliminate disk variance.
- Synthetic data uses a fixed seed (`42`) for reproducible compression ratios.
- Each benchmark case runs enough iterations for nanobench to reach statistical stability.

## CI

The `bench` job runs on `ubuntu-24.04 / gcc-14` in Release mode. It is **advisory only** — never blocks merge.
