# Vendored Third-Party Linear Algebra Dependencies

> **Audience:** This guide is for downstream packagers, CI maintainers, and H5CPP users choosing which linear algebra adapter to depend on.

H5CPP provides optional adapters for popular C++ linear algebra libraries, allowing seamless HDF5 I/O for matrix, vector, and array types. The adapters are located in `h5cpp/H5M*.hpp` and the corresponding library sources are vendored under `thirdparty/`.

## Design Philosophy

Our goal is to **prove that linear algebra objects work with H5CPP I/O**, not to build or ship full-featured linear algebra distributions. H5CPP adapters only need a small, stable subset of each library's API:

- `data()` / `memptr()` — raw pointer access to contiguous storage
- `size()` / `n_elem`, `n_rows`, `n_cols` — dimension queries
- `ctor(dims)` — constructors that allocate from dimensions
- `begin()` / `end()` — iterator access (where applicable)

Because our surface area is so small, many libraries that traditionally ship as compiled libraries can be used in **header-only mode**. This eliminates CMake build-time complexity, avoids version-lock issues with obsolete `cmake_minimum_required` declarations, and removes cross-platform wrapper-library headaches. The result is faster CI, fewer system-package dependencies, and reduced support burden.

---

## Library Classification: Header-Only vs Compiled

| Library | Vendored Version | H5CPP Adapter | Build Mode | Header-Only? | Rationale |
|---------|-----------------|---------------|------------|--------------|-----------|
| **Eigen3** | v3.4.0 | `H5Meigen.hpp` | `INTERFACE` | ✅ Yes | By design. Pure templates, no compiled runtime. |
| **Blaze** | v3.8 | `H5Mblaze.hpp` | `INTERFACE` | ✅ Yes | By design. Pure templates, no compiled runtime. |
| **Boost.uBLAS** | v1.74.0 | `H5Mublas.hpp` | `INTERFACE` | ✅ Yes | By design. Header-only subset of Boost. |
| **Armadillo** | **v15.2.6** | `H5Marma.hpp` | `INTERFACE` | ✅ **Yes**¹ | `ARMA_DONT_USE_WRAPPER` disables the BLAS/LAPACK wrapper library. Headers contain fallback implementations for basic operations (creation, element access, `memptr()`, `print()`). Verified by compile-and-run test without linking BLAS/LAPACK. |
| **Blitz++** | v1.0.2 | `H5Mblitz.hpp` | Compiled (`libblitz.a`) | ❌ No | Requires `globals.cpp` (global variable instantiation) and a generated `bzconfig.h` produced by its CMake configure step. |
| **dlib** | v19.24 | `H5Mdlib.hpp` | Compiled (`libdlib.a`) | ❌ No | Matrix headers pull in ODR-violation check symbols (`USER_ERROR__inconsistent_build_configuration__see_dlib_faq_*`) that require linking the compiled `libdlib`. |
| **IT++** | v4.3.1 | `H5Mitpp.hpp` | **Disabled** | ❌ No | Adapter code exists in `h5cpp/H5Mitpp.hpp` but is excluded from all CMake targets and test suites. It is not compiled or installed. Disabled because upstream v4.3.1 uses the `register` keyword extensively, which was removed in C++17, causing build failures with GCC 14. |

¹ **BLAS/LAPACK still required** for high-performance linear algebra operations (matrix multiplication, inversion, SVD, eigenvalue decomposition) if exercised by downstream code. H5CPP examples do not use these paths.

### Summary

- **4 of 7 libraries** are already or can be header-only (Eigen3, Blaze, uBLAS, Armadillo).
- **3 libraries genuinely need compiled libraries** (Blitz++, dlib, IT++).
- IT++ is **disabled** due to upstream C++17 incompatibility.

---

## Licensing & Compliance

H5CPP itself is distributed under the **MIT License**. Each vendored library retains its original license.

| Library | License | Copyleft | Vendored Path | Notes |
|---------|---------|----------|---------------|-------|
| Armadillo | Apache-2.0 | No | `thirdparty/armadillo/` | |
| Blaze | BSD-3-Clause | No | `thirdparty/blaze/` | |
| Blitz++ | Artistic-2.0 **or** BSD **or** LGPL-3.0+ | Yes (LGPL option) | `thirdparty/blitz/` | Triple-licensed; choose the license that best fits your distribution policy. |
| dlib | BSL-1.0 | No | `thirdparty/dlib/` | Boost Software License |
| Eigen3 | MPL-2.0 | Weak | `thirdparty/eigen3/` | File-level copyleft |
| IT++ | GPL-3.0+ | **Strong** | `thirdparty/itpp/` | **Disabled by default.** If you require a fully non-copyleft build, ensure IT++ sources are excluded from your distribution. |
| Boost.uBLAS | BSL-1.0 | No | `thirdparty/ublas/` | Boost Software License |

> **For Distribution Packagers:** H5CPP vendors these libraries solely so that examples compile out-of-the-box without system package manager dependencies. The core `h5cpp/` headers do not depend on any third-party linear algebra library. If your distribution policy prohibits vendored dependencies, the `thirdparty/` directory may be removed; the adapter headers in `h5cpp/H5M*.hpp` are self-contained and will compile against system packages when the corresponding headers are available on the include path.

---

## Vendoring Policy

### Why vendoring?

Vendoring guarantees that H5CPP examples compile out-of-the-box on any platform with a C++17 compiler and CMake 3.22+, without requiring users to install Eigen3, Armadillo, or Blaze through apt, Homebrew, vcpkg, or Conan. This is especially important for CI reproducibility and for users evaluating H5CPP for the first time.

### Future direction

Future releases may add `H5CPP_USE_SYSTEM_<LIB>` CMake options (e.g., `H5CPP_USE_SYSTEM_EIGEN3`) to allow builds against system packages when available. The default will remain vendored to preserve the zero-dependency experience.

### Update cadence

Vendored dependencies are reviewed quarterly. Security-critical updates are applied within 30 days of upstream release. Header-only libraries are preferred for new adapters because they eliminate compiled-library maintenance overhead.

### Library selection criteria

A new linear algebra adapter is considered when the target library meets at least two of the following:

1. **Popularity:** Top-tier usage in C++ scientific computing (GitHub stars, Conan/vcpkg downloads, citation count).
2. **Active maintenance:** Upstream releases within the last 24 months.
3. **Low build complexity:** Header-only preferred; compiled library acceptable only if build is CMake-native and cross-platform.
4. **User demand:** Explicit request via GitHub issue with demonstrated use case.

---

## Armadillo: v12.6.7 → v15.2.6 Header-Only Migration

### Motivation

Armadillo v12.6.7's `CMakeLists.txt` declared `cmake_minimum_required(VERSION 2.8.12)`, which CMake 4.x (shipped on macOS via Homebrew) rejects outright. This broke macOS CI runners. Rather than patch the obsolete minimum version in-place, we upgraded to the latest upstream release and switched to header-only usage.

### Business impact

- Eliminated a CMake 4.x configuration failure on macOS CI runners.
- Removed the need for CMake version pinning workarounds in CI.
- Reduced vendored Armadillo footprint from ~1,000+ files (including docs, tests, examples, MEX interface, and wrapper source) to ~650 headers (~6.7 MB as of v15.2.6).

### What changed

| Aspect | Before (v12.6.7) | After (v15.2.6) |
|--------|-----------------|-----------------|
| Distribution | Full source: CMake, docs, examples, `src/wrapper*.cpp`, tests, MEX interface | Headers only: `include/armadillo` + `include/armadillo_bits/` |
| Build | `add_subdirectory(v12.6.7)` → builds `libarmadillo` wrapper library | `add_library(libarmadillo INTERFACE)` with `ARMA_DONT_USE_WRAPPER` |
| CMake minimum | `2.8.12` (breaks CMake 4.x) | N/A — no CMake sub-build |
| Linking | Examples link `libarmadillo` (transitively pulls BLAS/LAPACK) | Examples link `libarmadillo` (INTERFACE target; provides includes + compile definition) |

### Verification

- ✅ `examples-arma` (CTest target) compiles and links.
- ✅ `examples-reference` (CTest target) compiles, links, and runs correctly.
- ✅ `examples-nbit` (CTest target) compiles and links.
- ✅ Standalone test (`arma::mat`, `M.ones()`, `M.print()`, `M.begin()`, `memptr()`) compiles and runs without BLAS/LAPACK.
- ✅ CMake configure succeeds (no `cmake_minimum_required` error).

### When BLAS/LAPACK are still needed

`ARMA_DONT_USE_WRAPPER` only affects whether Armadillo uses its runtime wrapper library. Operations that internally call BLAS/LAPACK (matrix multiplication, inversion, SVD, eigenvalue decomposition, etc.) still require linking `-lblas -llapack` (or `-framework Accelerate` on macOS) if used by downstream code. H5CPP's examples do not exercise these paths.

If your downstream code needs high-performance linear algebra from Armadillo, link BLAS/LAPACK explicitly:

```cmake
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(my_target BLAS::BLAS LAPACK::LAPACK)
```

---

## Systemic Risk: Other Libraries and CMake 4.x

Armadillo was the first vendored library to break on CMake 4.x, but it may not be the last. Blitz++ (v1.0.2) and dlib (v19.24) also declare legacy `cmake_minimum_required` versions in their own build systems. Because H5CPP builds these libraries via `execute_process()` during the CMake configure step, they are vulnerable to the same CMake 4.x rejection.

**Mitigation plan:**

- **Blitz++:** Adapter deprecated. Will be removed in a future release unless a maintainer steps forward. Upstream development has ceased.
- **dlib:** Candidate for upgrade to v20.x, which may have a higher `cmake_minimum_required`. Tracked as future work.

---

## Future Work

| Item | Status | Owner / Tracking |
|------|--------|-----------------|
| Blitz++ adapter removal | Planned | [#150](https://github.com/vargalabs/h5cpp/issues/150) — Adapter deprecated. Will be removed in a subsequent release unless a maintainer steps forward. |
| dlib upgrade to v20.x | Under evaluation | Evaluate v20.x for Clang 20 compatibility (`-Wmissing-template-arg-list-after-template-kw`) and CMake 4.x support. |
| IT++ adapter removal | Planned | Remove disabled `H5Mitpp.hpp` and `thirdparty/itpp/` from source tree in next major release. |
| System-package fallback | Proposed | Add `H5CPP_USE_SYSTEM_*` CMake options for header-only libraries (Eigen3, Blaze, Armadillo). |
| CMake 4.x CI gate | Proposed | Add a CI matrix entry with the latest CMake version to catch `cmake_minimum_required` breakages early. |

---

## Recommendations for Downstream Consumers

If you are building H5CPP from source and want to use the vendored linear algebra libraries:

1. **Header-only libraries** (Eigen3, Blaze, uBLAS, Armadillo) require no extra system dependencies for basic I/O.
2. **Compiled libraries** (Blitz++, dlib) are built automatically during the CMake configure step; no manual intervention is required.
3. **Armadillo consumers** who need high-performance linear algebra (matrix multiply, decompositions) should link BLAS/LAPACK explicitly in their own build.

---

*Document version: 2026-05-12. Armadillo version references are current as of v15.2.6.*
