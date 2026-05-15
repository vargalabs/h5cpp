---
skill: includes
tags: [includes, core, io, reflection]
hdf5_min: 1.10
tested: true
source: test/H5Dio_roundtrip.cpp, test/H5Dappend.cpp
---

# Include order and entry points

h5cpp has three entry points. The choice depends on whether compiler-assisted
struct reflection is needed.

## `<h5cpp/all>` — no reflection, simplest

Use when working with scalars, STL containers, and linalg types only.

```cpp
#include <h5cpp/all>
```

Equivalent to `core` + `io` with no injection point. Any compiler-generated
shim (`h5cpp-compiler`) cannot be sandwiched here.

## `core` → shim → `io` — with reflection (compiler plugin)

The include order is load-bearing. `core` builds the type/meta system.
The compiler-generated shim specializes into it. `io` then sees the
completed type registry.

```cpp
#include <h5cpp/core>
#include "mystructs.h"   // compiler-generated shim from h5cpp-compiler
#include <h5cpp/io>
```

## Rules

- Never reverse `core` / `io` — `io` depends on `core` being complete.
- `H5Tmeta.hpp` is the central trait file — only modify within the issue chain
  that owns the type system.
- `<h5cpp/all>` and the `core` + `io` sandwich are mutually exclusive per
  translation unit.

## Minimal test include (used in all ctests)

```cpp
#include <h5cpp/core>
#include <h5cpp/io>
```
