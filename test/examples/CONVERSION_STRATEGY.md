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

### Category A: Simple I/O (easy conversions)
- [x] `string/string.cpp` → `test_string_io.cpp`
- [x] `stl/vector.cpp` → `test_stl_vector_io.cpp`
- [ ] `basics/basics.cpp` → mostly API demonstration; may skip or test property creation only
- [ ] `raw_memory/raw.cpp` → raw pointer buffer round-trip
- [ ] `csv/csv2hdf5.cpp` → CSV parse + write, then read back and compare

### Category B: Compound / Custom Types
- [ ] `compound/struct.cpp` → write compound struct, read back, field-by-field compare
- [ ] `attributes/attributes.cpp` → write attrs, read back, verify values
- [ ] `datatypes/two-bit.cpp` → custom 2-bit type round-trip
- [ ] `datatypes/n-bit.cpp` → n-bit type with eigen3 dependency

### Category C: Linear Algebra (requires thirdparty libs)
These need `thirdparty/` (armadillo, eigen3, blaze, dlib, ublas, itpp, blitz).
- [ ] `linalg/arma.cpp` → arma::mat/vec round-trip
- [ ] `linalg/eigen3.cpp` → Eigen::Matrix round-trip
- [ ] `linalg/blaze.cpp` → blaze::DynamicMatrix round-trip
- [ ] `linalg/dlib.cpp` → dlib::matrix round-trip
- [ ] `linalg/ublas.cpp` → ublas::matrix round-trip
- [ ] `linalg/itpp.cpp` → itpp::mat round-trip
- [ ] `linalg/blitz.cpp` → blitz::Array round-trip

### Category D: Advanced Features
- [ ] `packet-table/packettable.cpp` → append + read back packet table
- [ ] `optimized/optimized.cpp` → verify optimized write path produces correct data
- [ ] `transform/transform.cpp` → transform-on-write/read round-trip
- [ ] `utf/utf.cpp` → UTF-8 string round-trip

### Category E: MPI (requires PHDF5 + mpirun)
- [ ] `mpi/collective.cpp` → `mpirun -n 2` collective write/read
- [ ] `mpi/independent.cpp` → independent write/read
- [ ] `mpi/throughput.cpp` → may be performance-only, hard to validate

### Category F: Multi-TU / Compiler-generated
- [ ] `multi-tu/` → compile-time test only; may not need runtime I/O test

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

1. **Linalg examples** — blocked until `thirdparty/` is merged to staging (issue #144).
2. **MPI examples** — need `ctest` launcher configuration for `mpirun`.
3. **Property-based testing** — replace fixed test data with randomized generators
   (reuse `h5::utils::get_test_data<T>(n)` where available).
4. **Coverage** — ensure these tests are included in Codecov reporting.
