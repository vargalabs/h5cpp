---
hide:
  - toc
---

# Design Considerations

## Guaranteed No Vendor Lock-In

H5CPP was designed around a simple principle:

> **make HDF5 easier to use in modern C++ without trapping the user inside a wrapper-specific world**

That means two things from the outset:

- :material-license:{.icon} **MIT licensed**
- :material-link-variant:{.icon} **fully interoperable with the native HDF5 C API**

Together, these two properties eliminate vendor lock-in in both the legal and technical sense. You are free to use H5CPP in open-source, research, commercial, and production systems, and you are never forced to abandon the native HDF5 programming model in order to benefit from the library.

## No lock-in at the license level

H5CPP is released under the **MIT license**.

That matters because infrastructure libraries tend to live for years and often become deeply embedded in research platforms, market-data pipelines, simulation engines, and HPC applications. A restrictive license at the storage layer is the sort of thing that quietly becomes tomorrow’s migraine.

With MIT licensing:

- :material-check-circle-outline:{.icon} H5CPP can be used in commercial and proprietary codebases
- :material-check-circle-outline:{.icon} internal forks and extensions are straightforward
- :material-check-circle-outline:{.icon} adoption does not force a future relicensing event
- :material-check-circle-outline:{.icon} teams retain full freedom over how the library is integrated and maintained

In plain English: no legal handcuffs, no “gotcha” later.

## No lock-in at the API level

Many wrapper libraries reduce boilerplate at the cost of introducing a private object model that users must fully adopt. Once code enters that abstraction, getting back to the underlying system becomes awkward or impossible. H5CPP was designed specifically to avoid that trap.  With one deliberate exception, **all H5CPP handle types are fully binary compatible with native HDF5 C API handles**. The library does not invent a disconnected runtime identity for files, groups, datasets, dataspaces, datatypes, references, or property lists. Instead, it remains aligned with the underlying HDF5 `hid_t` handle model. That means H5CPP can be used as a modern C++ layer **without giving up direct access to the HDF5 C API**.

Conceptually:

```cpp
h5::fd_t fd = h5::open("data.h5", H5F_ACC_RDONLY);

/* direct native HDF5 call */
H5Fflush(static_cast<hid_t>(fd), H5F_SCOPE_GLOBAL);
```

This is not a convenience feature added afterward. It is a foundational design decision.

## The deliberate exception: `h5::pt_t`

The one intentional exception is:

* :material-table-arrow-down:{.icon} `h5::pt_t` — the packet table abstraction

Unlike the rest of the handle layer, `h5::pt_t` is a higher-level convenience abstraction for append-oriented workflows. It is not a one-to-one mirror of a native HDF5 handle in the same strict sense as files, groups, datasets, dataspaces, datatypes, references, and property lists. That exception is deliberate and contained. Everywhere else, the design remains tightly aligned with the native HDF5 object model.

## Full interoperability with the HDF5 C API

Because H5CPP preserves handle compatibility with the native API, users retain **full interaction with the HDF5 C layer**. That matters in real systems, because no wrapper can or should try to predict every advanced use case. Users often need to access:

* :material-tune:{.icon} low-level property-list tuning
* :material-memory:{.icon} explicit storage and layout control
* :material-database-settings-outline:{.icon} advanced dataset and datatype features
* :material-lan-connect:{.icon} MPI and parallel HDF5 integration details
* :material-wrench-outline:{.icon} niche or newly added HDF5 functionality not yet wrapped at the C++ level

H5CPP does not block that path. It improves the common case while preserving the escape hatch. That is exactly what you want from a systems library: help where possible, and zero obstruction where precision matters.

## A thin C++ layer, not a competing runtime

H5CPP was never meant to create a second universe around HDF5. Its role is much narrower and more useful:

* :material-lightning-bolt:{.icon} reduce boilerplate
* :material-book-open-variant:{.icon} improve readability
* :material-language-cpp:{.icon} integrate naturally with templates, containers, and RAII
* :material-link:{.icon} remain directly connected to the native HDF5 system

That distinction matters. H5CPP should feel like a **typed, modern C++ layer over HDF5**, not a proprietary storage runtime with its own rules, identities, and blind spots. The storage format remains HDF5.
The underlying runtime remains HDF5. The low-level escape path remains HDF5. So users gain convenience without surrendering control.

## Preserve the HDF5 mental model

A compatibility-first design also means preserving the HDF5 mental model itself. H5CPP does not try to erase the underlying concepts:

* files
* groups
* datasets
* dataspaces
* datatypes
* attributes
* property lists
* references

Those ideas remain visible because they are the real structure of the system. That has obvious practical value:

* :material-school-outline:{.icon} developers familiar with HDF5 remain productive immediately
* :material-file-replace-outline:{.icon} existing HDF5 documentation stays relevant
* :material-swap-horizontal:{.icon} moving between H5CPP and the C API is straightforward
* :material-head-cog-outline:{.icon} debugging is easier because the wrapper has not disguised the underlying model

This is especially important in scientific computing, HPC, and quantitative research, where many users already know HDF5 from C, Fortran, Python, MATLAB, or Julia.

## RAII without sacrificing control

H5CPP brings modern C++ ergonomics to HDF5, but not at the price of opacity.

A central design goal was to combine:

* the lifetime safety of RAII
* the explicitness of the HDF5 C API
* full interoperability with native HDF5 handles

This gives users a better day-to-day interface:

* fewer manual cleanup paths
* less repetitive boilerplate
* safer resource handling
* more readable code

while keeping native HDF5 semantics visible and accessible.

That balance matters. In systems programming, abstraction that hides too much eventually turns into debugging tax, and nobody enjoys paying that invoice.

## Incremental adoption in real codebases

Another consequence of the no-lock-in design is that H5CPP supports **incremental adoption**. A team does not need to rewrite an entire HDF5 codebase just to start benefiting from modern C++ improvements. Instead, they can:

* keep existing HDF5 C code where it already works
* use H5CPP in newly written modules
* gradually replace repetitive low-level code
* continue calling native HDF5 functions wherever needed

That makes H5CPP practical for long-lived systems and production environments, not just toy examples and fresh repositories full of optimism.

## Why this matters in quantitative and HPC systems

This no-lock-in approach is particularly valuable in:

* :material-finance:{.icon} quantitative research platforms
* :material-pulse:{.icon} market-data storage systems
* :material-chart-scatter-plot:{.icon} numerical analysis pipelines
* :material-server-network:{.icon} HPC batch and MPI workflows

These environments tend to care deeply about:

* explicit control
* predictable performance
* interoperability
* stable storage formats
* incremental migration paths
* long-term maintainability

A wrapper that forces a proprietary abstraction boundary works against those goals. H5CPP instead gives users a modern C++ interface while keeping them anchored to the stable and widely deployed HDF5 ecosystem.

## Design summary

The core design choices are simple:

* :material-check-circle-outline:{.icon} **MIT licensed** for maximum freedom of use
* :material-check-circle-outline:{.icon} preserve the native HDF5 object model
* :material-check-circle-outline:{.icon} keep H5CPP handles binary compatible with HDF5 C API handles
* :material-check-circle-outline:{.icon} guarantee direct interoperability with the native C API
* :material-check-circle-outline:{.icon} use RAII to reduce cleanup boilerplate and improve safety
* :material-check-circle-outline:{.icon} keep `h5::pt_t` as the intentional higher-level exception
* :material-check-circle-outline:{.icon} support incremental adoption in existing HDF5 codebases

## Bottom line

H5CPP was designed to make HDF5 easier to use, **not harder to leave**. That is the real no-lock-in guarantee. The library is MIT licensed, stays fully grounded in the HDF5 C API, and preserves direct access to the native runtime. With the exception of the packet-table abstraction `h5::pt_t`, H5CPP handles remain binary compatible with HDF5 C API handles, so mixed C and C++ usage is not a workaround — it is part of the design.

