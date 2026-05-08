---
hide:
  - toc
---

H5CPP is designed as a **thin modern C++ layer over the native HDF5 handle model**, not as a replacement object system. Its goal is to reduce boilerplate, improve readability, and make HDF5 feel natural in modern C++ — **without cutting users off from the underlying HDF5 C API**.

Many wrapper libraries take the opposite approach: they reduce boilerplate by imposing a private object model, and once you adopt it, getting back to the underlying system becomes awkward or impossible. H5CPP was built to avoid that trap. With one deliberate exception, **all H5CPP handle types are fully binary compatible with native HDF5 C API handles**, remaining aligned with the underlying `hid_t` model rather than inventing a disconnected runtime identity for files, groups, datasets, dataspaces, datatypes, references, or property lists. The deliberate exception is `h5::pt_t`, a higher-level packet-table abstraction introduced where a dedicated wrapper genuinely improves usability.

This design keeps mixed C and C++ usage straightforward, preserves the relevance of existing HDF5 documentation and expertise, and supports incremental adoption in existing codebases. It also means that H5CPP reduces abstraction lock-in in the technical sense, while the **MIT license** removes it in the legal one. For infrastructure libraries that tend to live for years inside research platforms, market-data pipelines, simulation engines, and HPC systems, that is a design choice with long-term consequences.

---

<div class="flex gap-4 w-full" markdown="1">
  <div class="w-1/2 border border-solid border-gray-300 rounded-lg shadow-lg overflow-hidden flex flex-col" markdown="1">
  <div class="bg-gray-400 text-white px-4 h-9 rounded-t flex items-center" markdown="1">:fontawesome-solid-compass-drafting:{.icon} What H5CPP Adds</div>
  <div class="p-4 flex-1" markdown="1">

* :material-lightning-bolt:{.icon} reduces descriptor and cleanup boilerplate
* :material-book-open-variant:{.icon} makes HDF5 code easier to read and maintain
* :material-language-cpp:{.icon} fits naturally into modern C++ templates, containers, and RAII
* :material-shield-check:{.icon} improves safety through scoped resource management
* :material-school-outline:{.icon} keeps existing HDF5 knowledge immediately useful
* :material-table-arrow-down:{.icon} introduces higher-level abstractions such as `h5::pt_t` only where they provide clear value

  </div> </div>

  <div class="w-1/2 border border-solid border-gray-300 rounded-lg shadow-lg overflow-hidden flex flex-col" markdown="1">
  <div class="bg-gray-400 text-white px-4 h-9 rounded-t flex items-center" markdown="1">:material-link-variant:{.icon} What H5CPP Preserves</div>
  <div class="p-4 flex-1" markdown="1">

* :material-link:{.icon} direct interoperability with the native HDF5 C API
* :material-check-circle-outline:{.icon} alignment with the native `hid_t` handle model
* :material-tune:{.icon} access to low-level property-list tuning and layout control
* :material-database-settings-outline:{.icon} use of advanced dataset and datatype features without wrapper lock-in
* :material-swap-horizontal:{.icon} straightforward incremental adoption in existing HDF5 codebases
* :material-license:{.icon} MIT-licensed use in open-source, research, commercial, and internal systems

  </div></div>
</div>


## Bottom line

H5CPP was **designed** to make HDF5 easier to use **without replacing the native model** beneath it. That is the core design decision: modern C++ ergonomics on top, native HDF5 semantics underneath.
