# Example → I/O Test Conversion Strategy

## Goal
Convert the `examples/` directory programs into automated I/O round-trip tests
that run under `ctest`. The examples themselves remain as documentation; the
tests exercise the same code paths and verify correctness.

## Pattern

Each converted test follows this structure:

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/all>
#include <h5cpp/all>
#include <filesystem>

TEST_CASE("[example] <name> round-trip") {
    const char* filename = "test_<name>.h5";
    std::filesystem::remove(filename);

    // 1. BUILD test data
    // 2. WRITE with h5cpp
    // 3. READ back with h5cpp
    // 4. VERIFY with doctest CHECK/CHECK_EQ/Approx
    // 5. CLEANUP file

    std::filesystem::remove(filename);
}
```

## Categorization of Examples

### Category A: Simple I/O ✅
- [x] `string/string.cpp` → `test_string_io.cpp`
- [x] `stl/vector.cpp` → `test_stl_vector_io.cpp`
- [x] `basics/basics.cpp` → `test_basics_io.cpp`
- [x] `raw_memory/raw.cpp` → `test_raw_memory_io.cpp`
- [x] `csv/csv2hdf5.cpp` → `test_csv_io.cpp`

### Category B: Compound / Custom Types ✅
- [x] `compound/struct.cpp` → `test_compound_struct_io.cpp`
- [x] `attributes/attributes.cpp` → `test_attributes_io.cpp`
- [x] `datatypes/two-bit.cpp` → `test_two_bit_io.cpp`
- [x] `datatypes/n-bit.cpp` → `test_n_bit_io.cpp`

### Category C: Linear Algebra ✅
- [x] `linalg/arma.cpp` → `test_arma_io.cpp`
- [x] `linalg/eigen3.cpp` → `test_eigen3_io.cpp` *(note: eigen3 example uses `Eigen::Matrix`; test covered via Category A/B infra)*
- [x] `linalg/blaze.cpp` → `test_blaze_io.cpp`
- [x] `linalg/dlib.cpp` → `test_dlib_io.cpp`
- [x] `linalg/ublas.cpp` → `test_ublas_io.cpp`
- [x] `linalg/blitz.cpp` → `test_blitz_io.cpp`
- [ ] `linalg/itpp.cpp` → **skipped** — upstream IT++ v4.3.1 is not C++17 compatible

### Category D: Advanced Features ✅
- [x] `packet-table/packettable.cpp` → `test_packet_table_io.cpp`
- [x] `optimized/optimized.cpp` → `test_optimized_io.cpp`
- [x] `transform/transform.cpp` → `test_transform_io.cpp`
- [x] `utf/utf.cpp` → `test_utf_io.cpp`

### Category E: MPI 🔴 TODO
- [ ] `mpi/collective.cpp` → `mpirun -n 2` collective write/read
- [ ] `mpi/independent.cpp` → independent write/read
- [ ] `mpi/throughput.cpp` → may be performance-only, hard to validate

**Blockers for Category E:**
1. `test/CMakeLists.txt` has no `find_package(MPI)` — must mirror `examples/CMakeLists.txt` logic
2. CI workflow installs serial `libhdf5-dev` only — PHDF5 (`libhdf5-mpich-dev` or `openmpi`) not available
3. `ctest` needs launcher configuration: `MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`
4. `HDF5_IS_PARALLEL` must be true for MPI tests to be relevant

### Category F: Multi-TU / Compiler-generated ✅
- [x] `multi-tu/` → `test_multi_tu.cpp` (compile-and-run verification)

## CMake Integration

The `test/examples/CMakeLists.txt` provides `add_example_test()` which mirrors
`add_test_case()` from the main `test/CMakeLists.txt`. Each example test gets:
- Its own executable
- Registered with `add_test()` for `ctest`
- Same include/link settings as core tests

## Temporary File Handling

All tests write to `<test_name>.h5` in the build directory and remove it at the
end of the test case. If the test crashes, the file may be left behind; CI
runners are ephemeral so this is acceptable.

## Future Work

1. **Category E (MPI)** — blocked on CI PHDF5 availability and CMake MPI launcher config.
2. **Property-based testing** — replace fixed test data with randomized generators
   (reuse `h5::utils::get_test_data<T>(n)` where available).
3. **Coverage** — ensure these tests are included in Codecov reporting.
