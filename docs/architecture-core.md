---
hide:
  - toc
---

Within an HDF5 container, datasets may be stored using compact, chunked, or contiguous layouts. Objects are identified by slash-separated `(/)` paths. Non-leaf nodes are groups `(h5::gr_t)`, while terminal leaf nodes are datasets `(h5::ds_t)` and named types `(h5::dt_t)`. Groups, datasets, and named types may also carry attached attributes `(h5::att_t)`. At first glance, an HDF5 container resembles a regular file system, but with a richer and more specialized API.

#### :material-view-module:{.icon} Chunked Layout and Partial I/O
<img src="../assets/hdf5-dataset-layout-chunked.svg" alt="Chunked Layout"
   class="float-right w-2/5 border-2 border-solid border-red rounded-lg shadow-lg ml-4 transition-transform duration-500 hover:scale-110 bg-transparent"/>

Chunked layout is the standard way to work efficiently with large datasets when full-array I/O is unnecessary. Instead of storing data as one contiguous region, the dataset is divided into fixed-size chunks, which enables partial reads and writes and allows the use of filters such as compression. In H5CPP, chunked storage is requested by adding `h5::chunk{...}` to the dataset creation property list. This implicitly selects `h5::layout_chunked`.

```cpp
h5::ds_t ds = h5::create<double>(fd, "dataset", ...,
	h5::chunk{4, 8} | h5::fill_value<double>{3} | h5::gzip{9});
```

Here, `fd` is an open HDF5 file descriptor of type `h5::fd_t`, and the omitted arguments define the dataset extent. H5CPP supports practical partial I/O through hyperslab-style selections using `h5::offset{...}`, `h5::stride{...}`, and `h5::block{...}`. For example, let: `arma::mat M(20, 16);` To write `M` matrix into a larger dataset at a given offset, optionally with a stride, use: `h5::write(ds, M, h5::offset{4, 8}, h5::stride{2, 2});` H5CPP infers the memory address, datatype, and dimensions of supported objects and forwards them to the underlying HDF5 calls. When working with raw pointers, or with types not yet recognized by H5CPP, the logical shape must be specified explicitly with `h5::count{...}` as in `h5::write(ds, M.memptr(), h5::count{5, 10});` For common cases, dataset creation and full-object write can be expressed in a single call. The following creates a chunked dataset, applies a fill value and gzip compression, and writes the full contents of `M` matrix:

```cpp
h5::write(fd, "dataset", M, h5::chunk{4, 8} | h5::fill_value<double>{3} | h5::gzip{9});
```
See the [examples](examples.md) section for complete working cases.

#### :material-arrow-expand-horizontal:{.icon} Contiguous Layout and IO Access
<img src="../assets/hdf5-dataset-layout-contiguous.svg" alt="Contiguous Layout"
   class="float-left w-2/5 border-2 border-solid border-red rounded-lg shadow-lg mr-4 transition-transform duration-500 hover:scale-110 bg-transparent"/>

Contiguous layout is the simplest storage mode in HDF5. The dataset is stored as one continuous block on disk, which makes it a good fit when the entire object is typically written or read in a single operation. This layout works well for **small datasets** and for workloads that do not require compression, chunking, or partial I/O. It also keeps metadata overhead low, which can be beneficial when working with many small objects. When no filtering is requested and the object is written in one shot, H5CPP will generally choose this layout automatically.

**Example:** the following call opens `arma.h5`, creates a dataset with the appropriate dimensions and datatype, and writes the full content of the vector in a single operation.

```cpp
arma::vec V({1.,2.,3.,4.,5.,6.,7.,8.});
h5::write("arma.h5", "one shot create write", V);
```

To request contiguous layout explicitly, pass the `h5::contiguous` flag in the dataset creation property list.

The resulting dataset layout is conceptually equivalent to:

```text
DATASET "one shot create write" {
   DATATYPE  H5T_IEEE_F64LE
   DATASPACE  SIMPLE { ( 8 ) / ( 8 ) }
   STORAGE_LAYOUT {
      CONTIGUOUS
      SIZE 64
      OFFSET 5888
   }
   FILTERS {
      NONE
   }
   FILLVALUE {
      FILL_TIME H5D_FILL_TIME_IFSET
      VALUE  H5D_FILL_VALUE_DEFAULT
   }
   ALLOCATION_TIME {
      H5D_ALLOC_TIME_LATE
   }
}
```

The trade-off is straightforward: contiguous layout is simple and efficient for full-object I/O, but it does not support chunk-based filtering or efficient partial access. For small, dense datasets, that is often exactly the right choice.


#### :material-package-variant-closed:{.icon} Compact Layout
<img src="../assets/hdf5-dataset-layout-compact.svg" alt="Compact Layout"
   class="float-right w-2/5 border-2 border-solid border-red rounded-lg shadow-lg ml-4 transition-transform duration-500 hover:scale-110 bg-transparent"/>

Compact layout stores the entire dataset payload directly in the object header rather than in a separate data block. This minimizes indirection and can be very efficient for **small datasets** that are written and read as a whole. Because both metadata and data live together, access is simple and overhead is low. The trade-off is capacity and flexibility. Compact datasets are size-limited by the available object header space, so they are not suitable for larger arrays, append-style workflows, compression, or partial I/O. In practice, compact layout is most useful for very small fixed-size objects where minimizing storage overhead matters more than scalability. For anything expected to grow, be filtered, or accessed in pieces, chunked or contiguous layout is the better choice.


#### :material-ruler-square:{.icon} Dataspaces and Dimensions

A dataspace describes how data is shaped and mapped between memory and file storage. In practical terms, it tells HDF5 how an in-memory region corresponds to a dataset on disk, or how a region on disk should be read back into memory. For example, a contiguous block of memory may be interpreted as a vector, matrix, or cube-shaped dataset depending on the associated dataspace. A dataspace may have fixed extent, bounded extensible extent, or unlimited extent along one or more dimensions. When working with supported objects, H5CPP derives the in-memory dataspace automatically. When passing raw pointers to I/O operations, the relevant shape must be provided explicitly, and the file selection determines how much memory is transferred. The following descriptors define dataset dimensions:

* `h5::current_dims{i,j,k,...}` — current extent of the dataset
* `h5::max_dims{i,j,k,...}` — maximum extent; use `H5S_UNLIMITED` for an unbounded dimension
* `h5::chunk{i,j,k,...}` — chunk shape for chunked datasets; a suitable chunk layout can significantly improve performance

The following descriptors define read or write selections within a dataset:

* `h5::offset{i,j,k,...}` — starting coordinates of the selection
* `h5::stride{i,j,k,...}` — step between selected elements
* `h5::block{i,j,k,...}` — size of each selected block
* `h5::count{i,j,k,...}` — number of blocks to transfer

**Note:** `h5::stride`, `h5::block`, and scatter/gather-style selections are not available when `h5::high_throughput` is enabled, as that mode prioritizes simplified high-bandwidth transfer paths.


#  :material-tag-multiple-outline:{.icon} Attributes

Attributes are the standard HDF5 mechanism for attaching side-band metadata to groups (`h5::gr_t`), datasets (`h5::ds_t`), and named datatypes (`h5::dt_t`). In H5CPP, the same core I/O model applies here as well, so attributes support the same family of storage-capable object types described in the [Supported Types](#supported-objects) section. In the definitions below, `P ::= h5::ds_t | h5::gr_t | h5::ob_t | h5::dt_t` denotes any valid HDF5 object handle, while `std::tuple<T...>` represents a sequence of pairwise attribute operations, where consecutive tuple elements are interpreted as attribute-name / attribute-value pairs.

```cpp
h5::at_t acreate( const P& parent, const std::string& name, args_t&&... args );
h5::at_t aopen(const  P& parent, const std::string& name, const h5::acpl_t& acpl);
T aread( const P& parent, const std::string& name, const h5::acpl_t& acpl)
h5::att_t awrite( const P& parent, const std::string& name, const T& ref, const h5::acpl_t& acpl)
void awrite( const P& parent, const std::tuple<Field...>& fields, const h5::acpl_t& acpl)
```

**Example:** the following snippet creates three attributes and attaches them to an `h5::gr_t` group. It also illustrates mixed use of the HDF5 C API and H5CPP, including how a raw `hid_t` handle returned by the C API may be wrapped in the RAII-enabled `h5::gr_t` handle.
```cpp
auto attributes = std::make_tuple(
        "author", "steven varga",  "company","vargaconsulting.ca",  "year", 2019);
// freely mixed CAPI call returned `hid_t id` is wrapped in RAII capable H5CPP template class: 
h5::gr_t group{H5Gopen(fd,"/my-data-set", H5P_DEFAULT)};
// templates will synthesize code from std::tuple, resulting in 3 attributes created and written into `group` 
h5::awrite(group, attributes);
```