[![CI](https://github.com/vargalabs/h5cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/vargalabs/h5cpp/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/vargalabs/h5cpp/branch/main/graph/badge.svg)](https://app.codecov.io/gh/vargalabs/h5cpp/tree/main)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20123216.svg)](https://doi.org/10.5281/zenodo.20123216)
[![GitHub release](https://img.shields.io/github/v/release/vargalabs/h5cpp.svg)](https://github.com/vargalabs/h5cpp/releases)
[![Documentation](https://img.shields.io/badge/docs-stable-blue)](https://vargalabs.github.io/h5cpp)

Easy to use  [HDF5][hdf5] C++ templates for Serial and Paralell HDF5  
----------------------------------------------------------------------

## Build Matrix

| OS / Compiler | GCC 13        | GCC 14        | GCC 15      | Clang 17   | Clang 18     | Clang 19     | Clang 20     | Apple Clang | MSVC         |
|---------------|---------------|---------------|-----------|--------------|--------------|--------------|--------------|-------------|--------------|
| Ubuntu 22.04  | ![gcc13][200] | ![NA][NA]     | ![NA][NA] | ![cl17][250] | ![cl18][251] | ![cl19][252] | ![cl20][253] | ![NA][NA]   | ![NA][NA]    |
| Ubuntu 24.04  | ![gcc13][300] | ![gcc14][301] | ![NA][NA] | ![NA][NA]    | ![cl18][351] | ![cl19][352] | ![cl20][353] | ![NA][NA]   | ![NA][NA]    |
| macOS 15      | ![NA][NA]     | ![NA][NA]     | ![NA][NA] | ![NA][NA]    | ![NA][NA]    | ![NA][NA]    | ![NA][NA]    | ![ac][400]  | ![NA][NA]    |
| Windows       | ![NA][NA]     | ![NA][NA]     | ![NA][NA] | ![NA][NA]    | ![NA][NA]    | ![NA][NA]    | ![NA][NA]    | ![NA][NA]   | ![msvc][500] |


**H5CPP is a modern C++ template library for serial and parallel HDF5 I/O**, designed to make high-performance persistence natural for C++ applications. It provides type-safe RAII wrappers around HDF5 resources, high-level create, read, write, and append operations, and seamless interoperability with the native HDF5 C API. H5CPP supports primitive numeric types, POD/C/C++ structs, std::vector, std::string, and major linear algebra containers including Armadillo, Eigen, Blaze, Blitz++, Boost uBLAS, IT++, and dlib.

**The library includes compiler-assisted object persistence** for complex user-defined types, **enabling efficient storage of structured C++ data in portable HDF5 containers**. It supports chunked and compressed datasets, extendable packet-table style streams, hyperslab selection, **custom HDF5 datatypes, and MPI-enabled parallel HDF5 workflows using independent or collective I/O**.

H5CPP **is intended for** scientific computing, **high-performance computing**, machine learning, numerical analysis, and data-intensive C++ applications **that require portable**, inspectable, **and efficient binary storage** compatible with the broader HDF5 ecosystem, **including Python, R, MATLAB, Fortran, Julia, and other HDF5-capable environments**.



## Migration and Contribution Workflow

H5CPP is being migrated from its current stable state toward a modernized implementation with broader compiler coverage, improved platform support, and
new library features. The migration is expected to be uneventful for existing users, but regressions, portability issues, and missing feature coverage should
be tracked through GitHub issues. When you encounter a problem please file an issue in the pattern described under `issue naming`; if you have a fix, follow the steps further down.

Fixes and contributions are welcome, including external feature branches and pull requests. All contributions are accepted under the MIT license.

### Issue naming

Issues should follow the form:

```text
<group/category>, <description>
```

Examples:

```text
ci, add macOS AppleClang runner
fix, resolve HDF5 discovery on Homebrew
feature, add compiler metadata contract
refactor, simplify datatype synthesis
docs, update migration notes
```

The category should describe the intent of the work, for example:

```text
fix
feature
refactor
ci
docs
test
legal
cleanup
```

### Branch naming

Feature branches should be created from the current `staging` branch and follow
the form:

```text
<#issue>-<category>-<description>
```

Examples:

```text
138-ci-add-macos-appleclang-runner
139-fix-hdf5-homebrew-discovery
140-refactor-datatype-synthesis
```

This is the same style expected when creating an issue-linked development branch
with `gh issue develop <#issue>`. If the generated branch name differs, rename or
create the branch using the project convention.

New work should always start from the current `staging` branch:

```bash
git fetch origin
git checkout staging
git pull --ff-only origin staging
gh issue develop <#issue> --checkout
```

If needed, ensure the resulting branch follows the project naming convention.

### Commit format

Each commit should follow this pattern:

```text
[#issue]:author:category, description
```

Examples:

```text
[#138]:svarga:ci, add macOS AppleClang runner
[#139]:svarga:fix, resolve HDF5 discovery on Homebrew
[#140]:svarga:refactor, simplify datatype synthesis
```

Keep commits focused. Prefer several precise commits over one heroic commit
that tries to solve the universe and accidentally invents another one.

### Rebase-based development

H5CPP uses a rebase-based feature-branch workflow. Feature branches are carried
forward by rebasing onto the current `staging` branch rather than merging
`staging` into the feature branch.

Before opening or updating a pull request:

```bash
git fetch origin
git checkout <feature-branch>
git rebase origin/staging
```

Resolve conflicts locally, rerun the relevant build and test matrix where
possible, then force-push safely:

```bash
git push --force-with-lease
```

Use `--force-with-lease`, not plain `--force`; we are civilized barbarians.

### Pull requests

Once the feature branch is complete, open a pull request targeting `staging`.

A pull request should:

* reference the issue it resolves;
* follow the branch and commit naming convention;
* keep the change scoped to the issue;
* pass the relevant CI jobs;
* document user-visible behavior changes.

The preferred merge target for active development is `staging`. Release branches are promoted from `staging` after validation.


[NA]: https://vargalabs.github.io/h5cpp/docs/assets/na.svg

<!-- Ubuntu 22.04 -->
[200]: https://vargalabs.github.io/h5cpp/badges/ubuntu-22.04-gcc-13.svg
[250]: https://vargalabs.github.io/h5cpp/badges/ubuntu-22.04-clang-17.svg
[251]: https://vargalabs.github.io/h5cpp/badges/ubuntu-22.04-clang-18.svg
[252]: https://vargalabs.github.io/h5cpp/badges/ubuntu-22.04-clang-19.svg
[253]: https://vargalabs.github.io/h5cpp/badges/ubuntu-22.04-clang-20.svg

<!-- Ubuntu 24.04 -->
[300]: https://vargalabs.github.io/h5cpp/badges/ubuntu-24.04-gcc-13.svg
[301]: https://vargalabs.github.io/h5cpp/badges/ubuntu-24.04-gcc-14.svg
[351]: https://vargalabs.github.io/h5cpp/badges/ubuntu-24.04-clang-18.svg
[352]: https://vargalabs.github.io/h5cpp/badges/ubuntu-24.04-clang-19.svg
[353]: https://vargalabs.github.io/h5cpp/badges/ubuntu-24.04-clang-20.svg

<!-- macOS 15 -->
[400]: https://vargalabs.github.io/h5cpp/badges/macos-15-apple-clang.svg

<!-- Windows -->
[500]: https://vargalabs.github.io/h5cpp/badges/windows-latest-msvc.svg
