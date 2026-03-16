---
hide:
  - toc
---

H5CPP uses C++ template metaprogramming to build a Pythonic I/O interface: calls are composed from strongly typed descriptors rather than rigid positional argument lists, while remaining fully resolved at compile time and incurring no runtime overhead. The following sections describe the core H5CPP I/O operators: `h5::create`, `h5::read`, `h5::write`, `h5::append`, `h5::open`, `h5::close` and their syntax is presented in EBNF-style notation. We begin with a few common terms. 

- **File or container** An HDF5 file may be viewed as a container with its own internal hierarchy of groups, datasets, named types, and attributes. The file itself is an ordinary file in the host file system and can be copied, moved, or archived using standard tools, while its contents are accessed through HDF5-specific I/O calls. To refer to such a container, pass either an open file handle `h5::fd_t`, a group handle `h5::gr_t`, or the path to the HDF5 file:

```cpp
file ::= const h5::gr_t | const h5::fd_t& fd | const std::string& file_path;
```

- **Groups and datasets** HDF5 groups (`h5::gr_t`) are similar to file-system directories: they act as containers for datasets, named types, and subgroups. For this reason, the terms *group* and *directory* are used interchangeably in this documentation. Like datasets, groups may also carry attributes. A dataset is a concrete object within the container. To identify one uniquely, pass either an open dataset handle `h5::ds_t`, or specify both the container and the dataset path. In the latter case, H5CPP generates the necessary shim code at compile time to obtain the corresponding `h5::fd_t` handle and resolve the dataset.

```cpp
dataset ::= (const h5::fd_t& fd, const std::string& dataset_path)
          | (const std::string& file_path, const std::string& dataset_path)
          | const h5::ds_t& ds;
```

- **Dataspace** HDF5 datasets may have different shapes in memory and on disk. A **dataspace** describes that shape: it specifies the current extent of the object and, where applicable, the maximum extent to which it may be extended.

```cpp
dataspace ::= const h5::sp_t& dataspace
            | const h5::current_dims& current_dims
            | const h5::current_dims& current_dims,
              const h5::max_dims& max_dims
            | const h5::max_dims& max_dims;
```

- **Type parameter `T`** denotes the template type parameter associated with the object being read or written. In H5CPP, the element type is deduced at compile time, which allows the library to support a flexible range of storage-capable C++ types. Broadly speaking, supported objects fall into two categories:

    * types backed by contiguous memory, such as vectors, matrices, and C-style POD or standard-layout records
    * more complex C++ object types with non-trivial internal structure

The latter are not generally supported as first-class storage objects. See the [Supported Types](#types) section for details.

#### :material-folder-open-outline:{.icon}  Open
The [previous section](#pythonic-syntax) introduced the common EBNF terms such as `file` and `dataspace`. Opening HDF5 objects follows the same general pattern: each object type has a corresponding open function, optionally parameterized with a property list.

```cpp
h5::fd_t h5::open( const std::string& path,  H5F_ACC_RDWR | H5F_ACC_RDONLY [, const h5::fapl_t& fapl] );
h5::gr_t h5::gopen( const h5::fd_t | h5::gr_t& location, const std::string& path [, const h5::gapl_t& gapl] );
h5::ds_t h5::open( const h5::fd_t | h5::gr_t& location, const std::string& path [, const h5::dapl_t& dapl] );
h5::at_t h5::aopen(const h5:ds_t | h5::gr_t& node, const std::string& name [, const & acpl] );
```

The optional property lists include [`h5::fapl_t`][602] for files and [`h5::dapl_t`][605] for datasets. These allow access behavior to be refined when needed. The supported file access flags are:

* `H5F_ACC_RDWR` — open with read/write access
* `H5F_ACC_RDONLY` — open with read-only access

In the standard HDF5 model, only one process may open a file for writing at a time, while multiple readers may access it concurrently. This is the usual single-writer / multiple-reader constraint. If exclusive write access is too restrictive, use **parallel HDF5**. In that model, multiple MPI processes may write to the same file using collective or independent access patterns such as `h5::collective` or `h5::independent`. On MPI-based clusters and supercomputers, parallel HDF5 is the preferred approach and can scale to very high aggregate throughput.

Here is a cleaner revision:

#### :material-plus-box-outline:{.icon} Create

Creation follows the same pattern as opening: each HDF5 object type has a corresponding create function, optionally parameterized with property lists.

```cpp
[file]
h5::fd_t h5::create( const std::string& path, H5F_ACC_TRUNC | H5F_ACC_EXCL, 
            [, const h5::fcpl_t& fcpl] [, const h5::fapl_t& fapl]);
[group]
h5::fd_t h5::gcreate( const h5::fd_t | const h5::gr_t, const std::string& name
            [, const h5::lcpl_t& lcpl] [, const h5::gcpl_t& gcpl] [, const h5::gapl_t& gapl]);
[dataset]
template <typename T> h5::ds_t h5::create<T>( 
    const h5::fd_t | const h5::gr_t& location, const std::string& dataset_path, dataspace, 
    [, const h5::lcpl_t& lcpl] [, const h5::dcpl_t& dcpl] [, const h5::dapl_t& dapl] );
[attribute]
template <typename T>
h5::at_t acreate<T>( const h5::ds_t | const h5::gr_t& | const h5::dt_t& node, const std::string& name
     [, const h5::current_dims{...} ] [, const h5::acpl_t& acpl]);
```
The most commonly used property lists are: [`h5::fcpl_t`][601], [`h5::fapl_t`][602], [`h5::lcpl_t`][603], [`h5::dcpl_t`][604], and [`h5::dapl_t`][605]. The file creation flags are `H5F_ACC_TRUNC` — create the file, truncating it if it already exists, `H5F_ACC_EXCL` — create the file only if it does not already exist. The following creates an HDF5 file and then defines a chunked, extensible dataset.
```cpp
#include <h5cpp/all>
...
arma::mat M(2,3);
h5::fd_t fd = h5::create("arma.h5", H5F_ACC_TRUNC);
h5::ds_t ds = h5::create<short>(fd, "dataset/path/object.name",
    h5::current_dims{10, 20}, h5::max_dims{10, H5S_UNLIMITED}, h5::chunk{2, 3} | h5::fill_value<short>{3} | h5::gzip{9});
ds["attribute-name"] = std::vector<int>(10);
...
```

In this example:

* the file is created with truncation enabled
* the dataset is created with current dimensions `{10, 20}`
* the second dimension is declared extensible via `H5S_UNLIMITED`
* chunking, fill value, and gzip compression are specified through the dataset creation property list

A couple of technical notes.

* I changed `h5::fd_t` to `h5::gr_t` for `gcreate`, because returning a file descriptor there would be mildly alarming.
* `h5::at_t` may need to be `h5::att_t` if that is your actual attribute handle type elsewhere.
* `const h5::ds_t | const h5::gr_t | const h5::dt_t& node` is a bit uneven as written in the source; in real docs I would probably normalize all of those references visually.

The section is now much easier to scan: syntax first, property lists next, example after. Less soup, more structure.

Here is a cleaner revision of the section:

#### :material-database-arrow-down-outline:{.icon} Read

H5CPP provides two forms of read operators: one returns a **newly created object** for convenient one-shot access, while the other reads into **an existing object** for repeated use where avoiding allocation is important. When a value is returned, H5CPP creates the object with the appropriate shape and relies on return value optimization (RVO) so that no extra copy is introduced. When reading into an existing object, the caller is responsible for providing a correctly sized destination. H5CPP then uses the supplied memory directly whenever possible, avoiding unnecessary temporaries. Internally, HDF5 may still use a chunk-sized transfer buffer, for example for filtering or datatype conversion. In practice this overhead is usually small, and chunk sizes should generally be chosen so they remain cache-friendly.

```cpp
template <typename T> T h5::read( const h5::ds_t& ds
    [, const h5::offset_t& offset]  [, const h5::stride_t& stride] [, const h5::count_t& count]
    [, const h5::dxpl_t& dxpl ] ) const;
template <typename T> h5::err_t h5::read( const h5::ds_t& ds, T& ref 
    [, const [h5::offset_t& offset]  [, const h5::stride_t& stride] [, const h5::count_t& count]
    [, const h5::dxpl_t& dxpl ] ) const;

template <typename T> T aread( const h5::ds_t& | const h5::gr_t& | const h5::dt_t& node, 
    const std::string& name [, const h5::acpl_t& acpl]) const;
template <typename T> T aread( const h5::at_t& attr [, const h5::acpl_t& acpl]) const;
```

Dimensions of extent `1` are dropped automatically, so the result is returned as a 2D object. For repeated access inside a tight loop, read into a preallocated buffer:

```cpp
arma::mat buffer(3, 4);
for (int i = 0; i < 10'000'000; i++)
    h5::read(ds, buffer, h5::offset{i + chunk});
```

In this form, `h5::count{...}` is inferred from the size of `buffer`. Passing an explicit count that does not match the target object can produce a useful compile-time diagnostic. When working with raw or managed pointers, the target memory extent must be specified explicitly:

```cpp
std::unique_ptr<double[]> ptr(new double[buffer_size]);

for (int i = 0; i < 10'000'000; i++)
    h5::read(ds, ptr.get(), h5::offset{i + buffer_size}, h5::count{buffer_size});
```

#### :material-database-arrow-up-outline:{.icon} WRITE
H5CPP provides two forms of write operators. **Reference-based writes** are intended for object types the library understands directly, such as supported containers, matrices, and records, and **are generally the preferred** interface because they combine convenience with the same performance characteristics as pointer-based access. **Pointer-based writes** operate on arbitrary contiguous memory regions and **are** mainly **useful in cases where no higher-level type mapping is available**; in that case, you provide the memory address and transfer properties, and H5CPP handles the rest. As with reads, the underlying HDF5 library may still allocate a chunk-sized transfer buffer, typically for filtering or datatype conversion, but this overhead is usually small when chunk sizes are chosen sensibly—ideally no larger than the processor’s last-level cache.

```cpp
[dataset]
template <typename T> void h5::write( dataset,  const T& ref
    [,const h5::offset_t& offset] [,const h5::stride_t& stride]  [,const& h5::dxcpl_t& dxpl] );
template <typename T> void h5::write( dataset, const T* ptr
    [,const hsize_t* offset] [,const hsize_t* stride] ,const hsize_t* count [, const h5::dxpl_t dxpl ]);

[attribute]
template <typename T> void awrite( const h5::ds_t& | const h5::gr_t& | const h5::dt_t& node, 
    const std::string &name, const T& obj  [, const h5::acpl_t& acpl]);
template <typename T> void awrite( const h5::at_t& attr, const T* ptr [, const h5::acpl_t& acpl]);

```
Property lists are:  [`dxpl_t`][606]

```cpp
#include <Eigen/Dense>
#include <h5cpp/all>

h5::fd_t fd = h5::create("some_file.h5",H5F_ACC_TRUNC);
h5::write(fd,"/result",M);
```

####  :material-table-row-plus-after:{.icon} Append
For streaming or append-only workloads, **packet tables are** often the most **practical** choice. Although the append operator uses its own `h5::pt_t` handle, the underlying storage object is still an ordinary HDF5 dataset, and conversion between `h5::pt_t` and `h5::ds_t` is seamless. The main difference is internal: unlike the other H5CPP handles, `h5::pt_t` maintains its own buffer and uses a specialized transfer pipeline optimized for incremental writes. The same pipeline may also be enabled for regular dataset I/O through the `h5::experimental` data-transfer property, as described in the [experimental pipeline documentation][308].

``` c++
#include <h5cpp/core>
    #include "your_data_definition.h"
#include <h5cpp/io>
template <typename T> void h5::append(h5::pt_t& ds, const T& ref);
```

**Example:**
``` c++
#include <h5cpp/core>
    #include "your_data_definition.h"
#include <h5cpp/io>
auto fd = h5::create("NYSE high freq dataset.h5");
h5::pt_t pt = h5::create<ns::nyse_stock_quote>( fd, 
        "price_quotes/2018-01-05.qte",h5::max_dims{H5S_UNLIMITED}, h5::chunk{1024} | h5::gzip{9} );
quote_update_t qu;

bool having_a_good_day{true};
while( having_a_good_day ){
    try{
        recieve_data_from_udp_stream( qu )
        h5::append(pt, qu);
    } catch ( ... ){
      if( cant_fix_connection() )
            having_a_good_day = false; 
    }
}
```

