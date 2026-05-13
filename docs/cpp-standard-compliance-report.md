# H5CPP C++ Standard Compliance Report

**Audited by:** Winston (System Architect) 🏗️  
**Date:** 2026-05-12  
**Scope:** Header-only library (`h5cpp/*.hpp`), tests (`test/`), examples (`examples/`)  
**Baseline:** C++17 (`cxx_std_17`)  
**Total header files audited:** 40+ headers (~5,800 LOC in `h5cpp/`)

---

## Executive Summary

| Metric | Finding |
|--------|---------|
| **Current baseline** | C++17 — enforced in CMake with `CMAKE_CXX_STANDARD_REQUIRED ON` and `CMAKE_CXX_EXTENSIONS OFF` |
| **Overall compliance grade** | **B+** |
| **Highest-priority findings** | 1. `std::is_pod` deprecated in C++20 (7 hits, easy fix). 2. Two latent syntax/semantic bugs in headers. 3. C++11 compatibility shims are dead weight in a C++17 baseline. |

**Verdict:** H5CPP is a **solid C++17 citizen**. It uses `if constexpr`, fold expressions, `auto` return-type deduction, `std::disjunction_v/conjunction_v`, `std::bool_constant`, and `std::string_view` correctly. The library compiles cleanly under C++17 today. The path to C++20 is **low-friction** (mostly deprecation cleanup and modernization). C++23 and C++26 introduce **opportunities, not blockers**.

---

## C++17 Compliance Audit

### Violations Found

**Severity: Medium — Latent bugs**

| File | Line | Issue | Impact |
|------|------|-------|--------|
| `H5Pall.hpp` | 209 | **Syntax error:** `impl::fapl_args<hid_t,size_t,size_t,size_t, H5Pset_fapl_direct>` is missing a closing `>`. The line reads `...size_t, H5Pset_fapl_direct>;` | Compile failure if `H5_VERSION_GE(1,14,0)` is true. |
| `H5Pdcpl.hpp` | 19 | **Undeclared identifier:** `H5Pinsert2(dapl, ...)` uses `dapl` but the enclosing function takes `dcpl`. | Compile failure or silent misuse of an unrelated variable if this code path is ever instantiated. |
| `H5config.hpp` | 61–63 | **Duplicated macro guard:** `#ifndef H5CPP_ERROR_MSG` … `#define H5CPP_ERROR_MSG(...)` appears twice (lines 57–58 and 61–63). | Harmless but indicates macro-hygiene drift. |

> **Trade-off note:** These are not C++ standard violations per se, but they are *compliance gaps* — code that does not compile under any standard when the triggering macro is defined. Fix them before expanding HDF5 version support.

### Deprecated Usage

| Construct | Count | C++17 Status | C++20 Status | Notes |
|-----------|-------|--------------|--------------|-------|
| `std::is_pod<T>::value` | 7 | Valid | **Deprecated** | Used in `H5Tmeta.hpp`, `H5Mstl.hpp`, `H5Aread.hpp`, `H5misc.hpp`. Replace with `std::is_standard_layout_v<T> && std::is_trivial_v<T>`. |
| `NULL` | ~11 | Valid (but unidiomatic) | Valid (still unidiomatic) | Used in `H5Aread.hpp`, `H5Dread.hpp`, `H5Eall.hpp`, `H5capi.hpp`, `H5Dappend.hpp`, `H5Zpipeline.hpp`. Prefer `nullptr`. |
| C-style casts `(type)expr` | ~5 | Valid | Valid | Found in `H5Dappend.hpp` (`(void*) ptr_`), `H5Zpipeline.hpp` (`(char*)aligned_alloc(...)`). Prefer `static_cast`. |

### Best-Practice Gaps

1. **C++11 compatibility shims are redundant.**  
   `h5::compat::index_sequence`, `h5::compat::make_index_sequence`, `h5::meta::compat::void_t`, and `h5::meta::compat::is_detected` are all back-ports of C++14/17 features. Since the baseline is C++17, these add maintenance surface area with no benefit. `std::index_sequence`, `std::make_index_sequence`, and `std::void_t` are available natively. `std::is_detected` is not in C++17, but the library already uses its own compatible implementation — keep it, just drop the `compat` alias dance.

2. **`h5::impl::conditional` is a hand-rolled `std::conditional` clone.**  
   Found in `H5meta.hpp` (lines 174–177). Pure dead weight in C++17.

3. **`h5cpp__constexpr` / `h5cpp__assert` macros for pre-C++17 `if constexpr`.**  
   Found in `H5config.hpp` (lines 66–72). These exist to degrade `if constexpr` to runtime `if` on C++14. The project requires C++17 — remove them.

4. **`_v` suffix modernization debt.**  
   The codebase already uses `std::disjunction_v`, `std::conjunction_v`, `std::is_pointer_v`, `std::remove_cv_t`, and `std::bool_constant` in newer code (notably `H5Tmeta.hpp`), but legacy headers still write `std::is_same<...>::value` (32 hits), `std::is_integral<...>::value` (5 hits), `std::is_arithmetic<...>::value` (3 hits), `std::is_pointer<...>::value` (4 hits), `std::extent<T,N>::value` (2 hits), and `std::rank<T>::value` (2 hits). This is not a violation, but it is **inconsistent style** that increases visual noise.

5. **Macro hygiene — missing `#undef`.**
   - `H5Pall.hpp` defines `H5CPP__capicall` but `#undef`s `H5CPP__defid` (a typo — the macro is never defined, so the `#undef` is a no-op, and `H5CPP__capicall` leaks into the global macro namespace).
   - `H5Tall.hpp` defines `H5CPP_REGISTER_TYPE_` and `H5CPP_REGISTER_STRUCT` without `#undef`. They are meant to be used by consumers, so this is acceptable, but the macro names are in the global namespace without a scoped cleanup.

6. **`const static` non-integral variables in headers.**  
   `H5Pall.hpp`, `H5Pdcpl.hpp`, `H5Pdapl.hpp`, and `H5Sall.hpp` declare `const static` objects like `h5::dcpl_t dcpl`, `h5::layout compact`, etc. Pre-C++17, `const` namespace-scope variables have internal linkage, so this is technically safe from ODR. However, C++17 `inline` variables are the explicit, modern mechanism for header-only variables. No action required for C++17 compliance, but migrating to `inline constexpr` (where possible) improves clarity.

7. **Global-namespace `operator<<` overloads are `inline`.**  
   Good — no ODR risk. Found in `H5cout.hpp`, `H5Tall.hpp`, `H5Zpipeline.hpp`, `H5Dappend.hpp`.

### Recommendations

| Priority | Action | Effort |
|----------|--------|--------|
| P0 | Fix syntax error in `H5Pall.hpp:209` and identifier bug in `H5Pdcpl.hpp:19`. | Minutes |
| P1 | Replace `std::is_pod<T>::value` with `std::is_standard_layout_v<T> && std::is_trivial_v<T>`. | ~30 min |
| P1 | Replace `NULL` with `nullptr`. | ~15 min |
| P2 | Remove `h5::compat` index-sequence shims; use `std::index_sequence` directly. | ~1 hour |
| P2 | Remove `h5cpp__constexpr` / `h5cpp__assert` macros. | ~15 min |
| P2 | Fix `H5CPP__capicall` `#undef` typo. | Minutes |
| P3 | Bulk-replace legacy `::value` / `::type` trait accessors with `_v` / `_t` aliases. | ~2 hours |
| P3 | Replace C-style casts with `static_cast`/`reinterpret_cast`. | ~15 min |

---

## C++20 Readiness

### Breaking Changes (if any)

| Feature | H5CPP Exposure | Risk |
|---------|---------------|------|
| `std::is_pod` deprecated | 7 hits | **Low** — deprecation warnings only; the trait is not removed. |
| `std::allocator::construct` removed | 0 hits | None |
| `std::result_of` removed | 0 hits | None |
| `std::is_literal_type` removed | 0 hits | None |
| `std::uncaught_exception` removed | 0 hits | None |
| `std::auto_ptr` removed | 0 hits | None |
| `std::memory_order` enum names changed | 0 hits in H5CPP headers (doctest uses them, but that's third-party) | None |
| `throw()` removed | 0 hits in H5CPP headers (found in third-party `half.hpp`) | None |
| `std::basic_string::reserve()` return value | 0 hits | None |

**Conclusion:** H5CPP has **zero breaking changes** moving to C++20. The only friction is deprecation noise from `std::is_pod`.

### Features to Adopt (with trade-offs)

| C++20 Feature | Where it helps | Trade-off | Effort |
|---------------|---------------|-----------|--------|
| **`concept` / `requires`** | Replace the heavy SFINAE / `std::enable_if` machinery in `H5Tmeta.hpp`, `H5Dread.hpp`, `H5Dwrite.hpp`, `H5Aread.hpp`, `H5Awrite.hpp`. | **Pros:** Cleaner error messages, faster compile times, self-documenting constraints. **Cons:** Breaks C++17 compatibility; must be behind feature-test macros or a major-version bump. | Medium — refactoring the trait vocabulary into concepts is design work, not just syntax. |
| **`std::span<T>`** | The `write(const T* ptr, ...)` and `read(..., T* ptr)` overloads currently take raw pointers + `h5::count_t`. `std::span` would unify contiguous-buffer APIs. | **Pros:** Safer, self-describing bounds. **Cons:** Adds a new public-API surface; existing `T*` overloads must be preserved for backward compatibility. | Low — add `std::span` overloads as thin wrappers. |
| **`explicit(bool)`** | The `hid_t` conversion constructors in `H5Iall.hpp` use `H5CPP__EXPLICIT` macro to toggle `explicit` on/off. `explicit(bool)` would eliminate the macro entirely. | **Pros:** Removes a macro. **Cons:** C++20-only; needs macro fallback for C++17. | Low |
| **Designated initializers** | `h5::offset_t`, `h5::count_t` etc. are initialized with braced lists. Designated init would not add much here because these are homogeneous arrays. | Minimal benefit for this codebase. | N/A |
| **`using enum`** | The `filter_direction_t` enum in `H5Zpipeline.hpp` and HDF5 C enums like `H5D_layout_t` could be imported locally. | Minor readability win in switch statements. | Low |
| **`std::format`** | Error message construction in `H5config.hpp` uses `std::string` concatenation. `std::format` would be cleaner. | **Pros:** Type-safe formatting. **Cons:** `<format>` header is heavy; the project targets header-only minimalism. | Low priority. |
| **`consteval`** | `h5::meta::get_extent()` and `h5::meta::extent_to_string` are already `constexpr`. `consteval` would guarantee compile-time only. | **Pros:** Stronger contract. **Cons:** C++20-only; marginal value. | Low |
| **`std::bit_cast`** | Could replace some `reinterpret_cast` patterns when mapping POD types to HDF5 types. | **Pros:** Safer, well-defined. **Cons:** Only useful if the library does raw byte-level type punning (it mostly delegates to HDF5 CAPI). | Low priority. |

### Migration Effort Estimate

- **Minimal C++20 compile (fix warnings only):** ~1 hour (replace `std::is_pod`, silence any new warnings).
- **Adopt concepts for core traits:** ~1–2 days (design + refactor `H5Tmeta.hpp` vocabulary).
- **Add `std::span` APIs:** ~half a day.

---

## C++23 Readiness

### Breaking Changes (if any)

| Feature | H5CPP Exposure | Risk |
|---------|---------------|------|
| `std::aligned_storage` deprecated | 0 hits | None |
| `std::aligned_union` deprecated | 0 hits | None |
| `std::is_pod` still deprecated | Already noted | None |
| `throw()` still removed | Already noted | None |

**Conclusion:** No breaking changes.

### Features to Adopt (with trade-offs)

| C++23 Feature | Where it helps | Trade-off | Effort |
|---------------|---------------|-----------|--------|
| **`std::expected<T,E>`** | Many H5CPP functions return `h5::ds_t` and throw on error. `std::expected` would offer a non-throwing alternative for performance-critical or embedded consumers. | **Pros:** Zero-cost error paths, no exceptions. **Cons:** Major API design decision; would need parallel overloads or a policy-based design. | High — strategic, not tactical. |
| **`if consteval`** | The `constexpr` metaprogramming in `H5meta.hpp` and `H5Tmeta.hpp` is already `if constexpr`. `if consteval` is not needed here because the code does not branch on "am I in a constant-evaluated context?" | Not applicable. | N/A |
| **`static_assert(false)` relaxed** | `H5Dgather.hpp` uses `detail::dependent_false_t<T,E>` to force a compile error in an unconstrained `gather` fallback. C++23 allows `static_assert(false, "...")` directly in templates. | **Pros:** Removes boilerplate `dependent_false_t`. **Cons:** C++23-only; trivial macro-guard possible. | Minutes |
| **`deducing this` / explicit object parameter** | The CRTP `pipeline_t<Derived>` in `H5Zpipeline.hpp` and `filter_t<Derived>` in `H5Zall.hpp` could be simplified with deducing `this`, avoiding the `static_cast<Derived*>(this)` pattern. | **Pros:** Eliminates CRTP boilerplate and potential misuse. **Cons:** C++23-only; the current CRTP is well-tested and stable. | Medium — design change. |
| **`std::print` / `std::println`** | Could replace `std::cout << ...` in debug/diagnostic code (`H5cout.hpp`). | **Pros:** Cleaner. **Cons:** `<print>` requires formatting support; H5CPP is header-only and targets embedded/HPC where iostream may already be avoided. | Low priority. |
| **`std::mdspan`** | Multidimensional array descriptors (`h5::count_t`, `h5::offset_t`) are custom `array`-like types. `std::mdspan` is a standardized multidimensional span. | **Pros:** Interoperability with the broader C++23 ecosystem. **Cons:** `mdspan` is complex; H5CPP's current descriptors are simple, constexpr-friendly, and tightly coupled to HDF5 CAPI. | Low priority — watch, don't adopt yet. |
| **`std::generator`** | Not applicable to a synchronous I/O library. | N/A | N/A |

### Migration Effort Estimate

- **Tactical C++23 cleanups (`static_assert(false)`, `deducing this` experiments):** ~half a day.
- **Strategic `std::expected` error-handling path:** ~1 week (design + implementation + tests).

---

## C++26 Preview

### Anticipated Breaking Changes

As of the 2026-05-12 draft status, **no breaking changes** are anticipated that would affect H5CPP. The library does not use:
- `= delete(void)` (not applicable)
- Trigraphs (already extinct)
- `register` keyword
- Dynamic exception specifications

### Features to Watch

| C++26 Feature | Relevance to H5CPP |
|---------------|-------------------|
| **Pack indexing** | `h5::meta::detail::cat` (tuple concatenation) and the `static_for` fold-expression patterns could be simplified with pack indexing `Ts...[0]`. This would reduce template-metaprogramming complexity. |
| **`std::is_within_lifetime`** | Not applicable — the library does not do low-level lifetime introspection beyond RAII handles. |
| **Reflection (std::meta)** | If reflection lands in C++26, the manual `H5CPP_REGISTER_STRUCT` macro and the `h5::name<T>` specialization dance could be replaced with compile-time reflection. This would be a **game-changer** for HDF5 type registration. |

> **Trade-off note:** C++26 reflection is the single most exciting feature for H5CPP's long-term architecture. It could eliminate the need for the external `h5cpp` compiler tool for POD struct registration. However, wait for standardization before designing around it.

---

## Appendix: File-by-File Notes

### Core Metaprogramming (`H5meta.hpp`, `H5Tmeta.hpp`, `compat.hpp`)
- **`H5meta.hpp`:** Uses C++17 fold expressions correctly (`(callback(...), ...)`). Custom `conditional` and `cat` tuple-concatenation are present but redundant. `static_for` is a solid C++17 implementation.
- **`H5Tmeta.hpp`:** The newest code. Uses `std::bool_constant`, `std::remove_cv_t`, `std::remove_reference_t`, `std::disjunction_v`, `std::conjunction_v`, `std::is_pointer_v` — all C++17. This is the **model** for how the rest of the library should look.
- **`compat.hpp`:** Back-ports `index_sequence`, `make_index_sequence`, `void_t`, `is_detected`. All redundant in C++17 except `is_detected` (which is not in the standard). Recommend keeping `is_detected` but dropping the `index_sequence` / `void_t` shims.

### Type System (`H5Tall.hpp`)
- `H5CPP_REGISTER_TYPE_` macro is clean but leaks into the global macro namespace.
- `if constexpr` used correctly inside the macro to branch on pointer types.
- `std::is_pointer<C_TYPE>::value` should be `std::is_pointer_v<C_TYPE>` for consistency.

### Handle System (`H5Iall.hpp`)
- CRTP-like `hid_t` template is well-structured.
- `H5CPP__EXPLICIT` macro could become `explicit(bool)` in C++20.
- Copy/move assignment operators correctly guard against self-assignment.

### I/O Dispatch (`H5Dwrite.hpp`, `H5Dread.hpp`, `H5Dcreate.hpp`, `H5Dappend.hpp`)
- Heavy use of `if constexpr` for compile-time argument parsing (offset/stride/block/count presence). This is **exactly** what C++17 was designed for.
- `arg::tpos` and `arg::get` provide a clean "named optional arguments" pattern.
- `try` / `catch` blocks around large template bodies are present; this is valid C++17 but can bloat code size. Consider `noexcept` contracts in the future.

### Property Lists (`H5Pall.hpp`, `H5Pdcpl.hpp`, `H5Pdapl.hpp`)
- **Bug:** `H5Pall.hpp:209` syntax error (`impl::fapl_args<..., H5Pset_fapl_direct>` missing `>`).
- **Bug:** `H5Pdcpl.hpp:19` uses `dapl` instead of `dcpl`.
- `const static` factory objects (`h5::gzip`, `h5::chunk`, etc.) are idiomatic for this library but could be `inline constexpr` in C++17.

### Pipeline (`H5Zpipeline.hpp`, `H5Zall.hpp`)
- CRTP `pipeline_t<Derived>` and `filter_t<Derived>` are stable.
- `aligned_alloc` + C-style cast `(char*)` is present. Prefer `static_cast<char*>(aligned_alloc(...))`.
- `std::unique_ptr<char>` with custom `free` deleter is used; this is correct but `std::unique_ptr<char, h5::impl::free>` would be more explicit.

### Error Handling (`H5Eall.hpp`)
- `thread_local` function pointer and `void*` for HDF5 error stack callback — correct and thread-safe.
- Exception hierarchy is deep but well-organized (`h5::error::io::dataset::write`, etc.).

### Output (`H5cout.hpp`)
- `operator<<` overloads are `inline` and in the global namespace — safe from ODR.
- `std::ostream& operator<<(std::ostream&, const std::vector<T>&)` is defined globally. This **could** clash with user-defined overloads in other libraries. Consider moving it to `namespace h5` and relying on ADL, or scoping it.

### Linear Algebra Adapters (`H5Marma.hpp`, `H5Meigen.hpp`, etc.)
- Clean specialization of `h5::impl::data`, `h5::impl::size`, `h5::impl::rank`, `h5::impl::get`.
- `const_cast` in Eigen adapter (`H5Meigen.hpp:39,43`) is noted with a `// TODO: remove const_cast`. This is a code-quality issue, not a standard-compliance issue.

### Tests (`test/H5compat.cpp`, `test/H5meta.cpp`)
- Tests already use `std::is_same_v`, confirming the project's *intent* is to use modern idioms — the legacy headers just haven't been fully modernized yet.
- `test/H5compat.cpp` tests the back-port shims. If the shims are removed, these tests become unnecessary.

---

*End of report.*
