---
hide:
  - toc
---

**H5CPP** is a modern C++ interface to HDF5 for persisting vectors, matrices, tensors, strings, and structs with concise template-based code, while preserving full interoperability with the native HDF5 C API. In practice, it removes much of the needless plumbing that direct HDF5 use in C and C++ so often entails. From research and simulation to market data and production, storage lies on the critical path of many quantitative systems. With LLVM-based compiler-assisted reflection, even complex structs can often be handled with a one-line expression. Direct HDF5 use in C and C++ often means:

**H5CPP began from a practical requirement:** efficient storage for large numerical datasets with both indexed block access and sequential streaming. Existing serialization systems handled streams reasonably well, but not the kind of partial, multidimensional, file-backed access needed in numerical computing and financial engineering. HDF5 already had most of the right storage primitives: partial I/O, extendable datasets, compression, and broad interoperability across operating systems and scientific environments. What was missing was a modern C++ interface.

The earliest implementation began as a collection of templates. Later, contact with The HDF Group helped shape the first H5CPP11 project, followed in time by the C++17 version. I am especially thankful to Gerd Heber for his sustained guidance and generosity over the years; to Elena Pourmal and David Pareah for their encouragement, support, and influence on the project’s direction; and, from Fermilab, to Mark Paterno and Chris for their thoughtful input. Many of the stronger ideas in H5CPP were sharpened through those discussions; any mistakes, omissions, or rough edges are entirely my own.

<div class="flex gap-4 w-full" markdown="1">
  <div class="w-1/2 border border-solid border-gray-300 rounded-lg shadow-lg overflow-hidden flex flex-col" markdown="1">
  <div class="bg-gray-400 text-white px-4 h-9 rounded-t flex items-center" markdown="1">:material-file-document-edit-outline: What You End Up Writing by Hand</div>
  <div class="p-4 flex-1" markdown="1">

  * :material-vector-square:{.icon} manual datatype and dataspace construction
  * :material-ruler-square:{.icon} repetitive shape bookkeeping
  * :material-code-braces:{.icon} verbose read/write setup
  * :material-lifebuoy:{.icon} resource lifetime management
  * :material-call-split:{.icon} impedance mismatch with modern C++ container types
  * :material-lan-connect:{.icon} added complexity when moving to MPI-enabled parallel I/O
  </div> </div>

  <div class="w-1/2 border border-solid border-gray-300 rounded-lg shadow-lg overflow-hidden flex flex-col" markdown="1">
  <div class="bg-gray-400 text-white px-4 h-9 rounded-t flex items-center" markdown="1">:material-alert-circle-outline: Where That Complexity Shows Up</div>
  <div class="p-4 flex-1" markdown="1">
  * :material-file-code-outline:{.icon} structured event data such as trades, quotes, fills, and order updates
  * :material-grid:{.icon} dense numerical datasets such as returns, features, factors, and signal matrices
  * :material-chart-line:{.icon} simulation, backtesting, and numerical research code
  * :material-archive-outline:{.icon} shared datasets reused across research, replay, and production workflows
  </div></div>
</div>

<div class="mt-4 w-full border border-solid border-gray-300 rounded-lg shadow-lg overflow-hidden flex flex-col" markdown="1">
<div class="bg-green-700 text-white px-4 h-9 rounded-t flex items-center" markdown="1">:material-server: What You Get with H5CPP</div>
<div class="flex w-full" markdown="1">
<div class="w-1/2 p-0" markdown="1">

* :material-lightning-bolt:{.icon} **fast to integrate**
* :material-language-cpp:{.icon} **natural in modern C++**
* :material-shield-check:{.icon} **safer through RAII**
</div>
<div class="w-1/2 p-0" markdown="1">

* :material-file-tree:{.icon} **well suited to structured records**
* :material-table-large:{.icon} **well suited to dense numerical arrays**
* :material-lan-connect:{.icon} **ready to scale from local workflows to parallel I/O**
</div></div></div>




H5CPP [has been presented](https://steven-varga.ca/site/talks/) in HDF5 and C++ community venues over multiple years, including HUG sessions, HDF Group events, C++ community talks, and ISC-related material. Topics include compiler-assisted reflection, POD introspection, MPI/parallel I/O, throughput/latency trade-offs, and practical HDF5 workflows.

