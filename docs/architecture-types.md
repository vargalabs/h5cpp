---
hide:
  - toc
---

# Type System
At the heart of H5CPP lies the type mapping mechanism to HDF5 NATIVE types. All type requests are redirected to this segment in one way or another. That includes supported vectors, matrices, cubes, C like structs, etc. While HDF5 internally supports type translations among various binary representation, H5CPP deals only with native types. This is not a violation of the HDF5 use-anywhere policy, only type conversion is delegated to hosts with different binary representations. Since the most common processors are Intel and AMD, with this approach conversion is unnescessary most of the time. In summary, H5CPP uses NATIVE types exclusively.

```yacc
integral 		:= [ unsigned | signed ] [int_8 | int_16 | int_32 | int_64 | float | double ] 
vectors  		:=  *integral
rugged_arrays 	:= **integral
string 			:= **char
linalg 			:= armadillo | eigen | ... 
scalar 			:= integral | pod_struct | string

# not handled yet: long double, complex, specialty types
```

Here is the relevant part responsible for type mapping:
```cpp
#define H5CPP_REGISTER_TYPE_( C_TYPE, H5_TYPE )                                           \
namespace h5 { namespace impl { namespace detail { 	                                      \
	template <> struct hid_t<C_TYPE,H5Tclose,true,true,hdf5::type> : public dt_p<C_TYPE> {\
		using parent = dt_p<C_TYPE>;                                                      \
		using parent::hid_t;                                                              \
		using hidtype = C_TYPE;                                                           \
		hid_t() : parent( H5Tcopy( H5_TYPE ) ) { 										  \
			hid_t id = static_cast<hid_t>( *this );                                       \
			if constexpr ( std::is_pointer<C_TYPE>::value )                               \
					H5Tset_size (id,H5T_VARIABLE), H5Tset_cset(id, H5T_CSET_UTF8);        \
		}                                                                                 \
	};                                                                                    \
}}}                                                                                       \
namespace h5 {                                                                            \
	template <> struct name<C_TYPE> {                                                     \
		static constexpr char const * value = #C_TYPE;                                    \
	};                                                                                    \
}                                                                                         \
```
Arithmetic types are associated with their NATIVE HDF5 equivalent:
```cpp
H5CPP_REGISTER_TYPE_(bool,H5T_NATIVE_HBOOL)

H5CPP_REGISTER_TYPE_(unsigned char, H5T_NATIVE_UCHAR) 			H5CPP_REGISTER_TYPE_(char, H5T_NATIVE_CHAR)
H5CPP_REGISTER_TYPE_(unsigned short, H5T_NATIVE_USHORT) 		H5CPP_REGISTER_TYPE_(short, H5T_NATIVE_SHORT)
H5CPP_REGISTER_TYPE_(unsigned int, H5T_NATIVE_UINT) 			H5CPP_REGISTER_TYPE_(int, H5T_NATIVE_INT)
H5CPP_REGISTER_TYPE_(unsigned long int, H5T_NATIVE_ULONG) 		H5CPP_REGISTER_TYPE_(long int, H5T_NATIVE_LONG)
H5CPP_REGISTER_TYPE_(unsigned long long int, H5T_NATIVE_ULLONG) H5CPP_REGISTER_TYPE_(long long int, H5T_NATIVE_LLONG)
H5CPP_REGISTER_TYPE_(float, H5T_NATIVE_FLOAT) 					H5CPP_REGISTER_TYPE_(double, H5T_NATIVE_DOUBLE)
H5CPP_REGISTER_TYPE_(long double,H5T_NATIVE_LDOUBLE)

H5CPP_REGISTER_TYPE_(char*, H5T_C_S1)
```
Record/POD struct types are registered through this macro:
```cpp
#define H5CPP_REGISTER_STRUCT( POD_STRUCT ) \
	H5CPP_REGISTER_TYPE_( POD_STRUCT, h5::register_struct<POD_STRUCT>() )
```
**FYI:** there are no other public/unregistered macros other than `H5CPP_REGISTER_STRUCT`

### Resource Handles and CAPI Interop
By default the `hid_t` type is automatically converted to / from H5CPP `h5::hid_t<T>` templated identifiers. All HDF5 CAPI identifiers are wrapped via the `h5::impl::hid_t<T>` internal template, maintaining binary compatibility, with the exception of `h5::pt_t` packet table handle.


Here are some properties of descriptors,`h5::ds_t` and `h5::open` are used only for the demonstration, replace them with arbitrary descriptors.

* `h5::ds_t ds` default CTOR, is initialized with `H5I_UNINIT`
* `h5::ds_t res = h5::open(...)` RVO  with RAII maneged resource `H5Iinc_ref` not called, best way to initialize
* `h5::ds_t res = std::move(h5::ds_t(...))` move assignment with RAII enabled for each handles, `H5Iinc_ref` called once
* `h5::ds_t res = h5::ds_t(...)` copy assignment with RAII enabled for each handles, `H5Iinc_ref` called twice, considered inexpensive
* `h5::ds_t res(hid_t)` creates resource with RAII, increments reference counts of `hid_t` with `H5Iinc_ref` 
* `h5::ds_t res{hid_t}` create resource without RAII, for managed CAPI handles, try not using it, a full copy of a handle is cheap

The Packet Table interface has an additional conversion CTOR defined from `h5::ds_t` to `h5::pt_t` and is the preffered way to obtain a handle:

* `h5::pt_t pt = h5::open(fd, 'path/to/dataset', ...)` obtains a `h5::ds_t` resource, converts it to `h5::pt_t` type then assigns it to names variable, this process will result in 2 increments to underlying CAPI handle or 2 calls to `H5Iinc_ref`. 
*  `h5::pt_t pt = h5::ds_t(...)` this pattern is used when creating packet table from existing and opened `h5::ds_t`, resources are copied (not stolen), in other words `h5::ds_t` remains valid after the operation, all state changes will be reflected.

During the lifespan of H5CPP handles all of the underlying HDF5 desciptors are **properly initialized, reference count maintained and closed.**

Resources may be grouped into `handles` and `property_lists`-s. 
```yacc
T := [ handles | property_list ]
handles   := [ fd_t | ds_t | att_t | err_t | grp_t | id_t | obj_t ]
property_lists := [ file | dataset | attrib | group | link | string | type | object ]

#            create       access       transfer     copy 
file    := [ h5::fcpl_t | h5::fapl_t                            ] 
dataset := [ h5::dcpl_t | h5::dapl_t | h5::dxpl_t               ]
attrib  := [ h5::acpl_t                                         ] 
group   := [ h5::gcpl_t | h5::gapl_t                            ]
link    := [ h5::lcpl_t | h5::lapl_t                            ]
string  := [              h5::scpl_t                            ] 
type    := [              h5::tapl_t                            ]
object  := [ h5::ocpl_t                           | h5::ocpyl_t ]
```
# Property Lists

The functions, macros, and subroutines listed here are used to manipulate property list objects in various ways, including to reset property values. With the use of property lists, HDF5 functions have been implemented and can be used in applications with fewer parameters than would be required without property lists, this mechanism is similar to [POSIX fcntl][700]. Properties are grouped into classes, and each class member may be daisy chained to obtain a property list.

Here is an example of how to obtain a data creation property list with chunk, fill value, shuffling, nbit, fletcher32 filters and gzip compression set:
```cpp
h5::dcpl_t dcpl = h5::chunk{2,3} 
  | h5::fill_value<short>{42} | h5::fletcher32 | h5::shuffle | h5::nbit | h5::gzip{9};
auto ds = h5::create("container.h5","/my/dataset.dat", h5::create_path | h5::utf8, dcpl, h5::default_dapl);
```
Properties may be passed in arbitrary order, by reference, or directly by daisy chaining them. The following property list descriptors are available:
```yacc
#            create       access       transfer     copy 
file    := [ h5::fcpl_t | h5::fapl_t                            ] 
dataset := [ h5::dcpl_t | h5::dapl_t | h5::dxpl_t               ]
attrib  := [ h5::acpl_t                                         ] 
group   := [ h5::gcpl_t | h5::gapl_t                            ]
link    := [ h5::lcpl_t | h5::lapl_t                            ]
string  := [              h5::scpl_t                            ] 
type    := [              h5::tapl_t                            ]
object  := [ h5::ocpl_t                           | h5::ocpyl_t ]
```

#### Default Properties:
Set to a non-default value (different from HDF5 CAPI):

* `h5::default_lcpl =  h5::utf8 | h5::create_intermediate_group;`

Reset to default (same as HDF5 CAPI):

* `h5::default_acpl`, `h5::default_dcpl`, `h5::default_dxpl`, `h5::default_fapl`,  `h5::default_fcpl`


## File Operations
#### [File Creation Property List][1001]
```cpp
// you may pass CAPI property list descriptors whose properties are daisy chained with the '|' operator 
auto fd = h5::create("002.h5", H5F_ACC_TRUNC, 
    h5::file_space_page_size{4096} | h5::userblock{512},  // file creation properties
    h5::fclose_degree_weak | h5::fapl_core{2048,1} );     // file access properties
```

* [`h5::userblock{hsize_t}`][1001] Sets the user block size of a file creation property list
* [`h5::sizes{size_t,size_t}`][1002] Sets the byte size of the offsets and lengths used to address objects in an HDF5 file.
* [`h5::sym_k{unsigned,unsigned}`][1003] Sets the size of parameters used to control the symbol table nodes.
* [`h5::istore_k{unsigned}`][1004] Sets the size of the parameter used to control the B-trees for indexing chunked datasets.
* [`h5::file_space_page_size{hsize_t}`][1005] Sets the file space page size for a file creation property list.
* [`h5::file_space_page_strategy{H5F_fspace_strategy_t strategy, hbool_t persist, hsize_t threshold}`][1010] Sets the file space handling strategy and persisting free-space values for a file creation property list. <br/>
* [`h5::shared_mesg_nindexes{unsigned}`][1007] Sets number of shared object header message indexes.
* [`h5::shared_mesg_index{unsigned,unsigned,unsigned}`][1008] Configures the specified shared object header message index.
* [`h5::shared_mesg_phase_change{unsigned,unsigned}`][1009] Sets shared object header message storage phase change thresholds.

#### [File Access Property List][1020]
**Example:**
```cpp
h5::fapl_t fapl = h5::fclose_degree_weak | h5::fapl_core{2048,1} | h5::core_write_tracking{false,1} 
      | h5::fapl_family{H5F_FAMILY_DEFAULT,0};
      
```
* [`h5::driver{hid_t new_driver_id, const void *new_driver_info}`][1055] Sets a file driver.<br/>
* [`h5::fclose_degree{H5F_close_degree_t}`][1022] Sets the file close degree.<br/>
  **Flags:** `h5::fclose_degree_weak`, `h5::fclose_degree_semi`, `h5::fclose_degree_strong`, `h5::fclose_degree_default`
* [`h5::fapl_core{size_t increment, hbool_t backing_store}`][1023] Modifies the file access property list to use the H5FD_CORE driver.
* [`h5::core_write_tracking{hbool_t is_enabled, size_t page_size}`][1024] Sets write tracking information for core driver, H5FD_CORE. 
* [`h5::fapl_direct{size_t alignment, size_t block_size, size_t cbuf_size}`][1025] Sets up use of the direct I/O driver.
* [`h5::fapl_family{hsize_t memb_size, hid_t memb_fapl_id}`][1026] Sets the file access property list to use the family driver.
* [`h5::family_offset{hsize_t offset}`][1027] Sets offset property for low-level access to a file in a family of files.
* [`h5::fapl_log{const char *logfile, unsigned long long flags, size_t buf_size}`][1028] Sets up the logging virtual file driver (H5FD_LOG) for use.
* [`h5::fapl_mpiio{MPI_Comm comm, MPI_Info info}`][1029] Stores MPI IO communicator information to the file access property list.
* [`h5::multi{const H5FD_mem_t *memb_map, const hid_t *memb_fapl, const char * const *memb_name, const haddr_t *memb_addr, hbool_t relax}`][1030] Sets up use of the multi-file driver.
* [`h5::multi_type{H5FD_mem_t}`][1031] Specifies type of data to be accessed via the MULTI driver, enabling more direct access. <br>
**Flags:** `h5::multi_type_super`, `h5::multi_type_btree`, `h5::multi_type_draw`, `h5::multi_type_gheap`, `h5::multi_type_lheap`, `h5::multi_type_ohdr`
* [`h5::fapl_split{const char *meta_ext, hid_t meta_plist_id, const char *raw_ext, hid_t raw_plist_id}`][1032] Emulates the old split file driver.
* [`h5::sec2`][1033] (flag) Sets the sec2 driver.
* [`h5::stdio`][1034] (flag) Sets the standard I/O driver.
* [`h5::windows`][1035] (flag) Sets the Windows I/O driver.
* [`h5::file_image{void*,size_t}`][1036] Sets an initial file image in a memory buffer.
* [`h5::file_image_callback{H5_file_image_callbacks_t *callbacks_ptr}`][1037] Sets the callbacks for working with file images.
* [`h5::meta_block_size{hsize_t}`][1038] Sets the minimum metadata block size.
* [`h5::page_buffer_size{size_t,unsigned,unsigned}`][1039] Sets the maximum size for the page buffer and the minimum percentage for metadata and raw data pages.
* [`h5::sieve_buf_size{size_t}`][1040] Sets the maximum size of the data sieve buffer.
* [`h5::alignment{hsize_t, hsize_t}`][1041] Sets alignment properties of a file access property list.
* [`h5::cache{int,size_t,size_t,double}`][1042] Sets the raw data chunk cache parameters.
* [`h5::elink_file_cache_size{unsigned}`][1043] Sets the number of files that can be held open in an external link open file cache.
* [`h5::evict_on_close{hbool_t}`][1044] Controls the library's behavior of evicting metadata associated with a closed object
* [`h5::metadata_read_attempts{unsigned}`][1045] Sets the number of read attempts in a file access property list.
* [`h5::mdc_config{H5AC_cache_config_t*}`][1046] Set the initial metadata cache configuration in the indicated File Access Property List to the supplied value.
* [`h5::mdc_image_config{H5AC_cache_image_config_t * config_ptr}`][1047] Sets the metadata cache image option for a file access property list.
* [`h5::mdc_log_options{const char*,hbool_t}`][1048] Sets metadata cache logging options.
* [`h5::all_coll_metadata_ops{hbool_t}`][1049] Sets metadata I/O mode for read operations to collective or independent (default).
* [`h5::coll_metadata_write{hbool_t}`][1050] Sets metadata write mode to collective or independent (default).
* [`h5::gc_references{unsigned}`][1050] H5Pset_gc_references sets the flag for garbage collecting references for the file.
* [`h5::small_data_block_size{hsize_t}`][1052] Sets the size of a contiguous block reserved for small data.
* [`h5::libver_bounds {H5F_libver_t,H5F_libver_t}`][1053] Sets bounds on library versions, and indirectly format versions, to be used when creating objects.
* [`h5::object_flush_cb{H5F_flush_cb_t,void*}`][1054] Sets a callback function to invoke when an object flush occurs in the file.
* **`h5::fapl_rest_vol`** or to request KITA/RestVOL services both flags are interchangeable you only need to specify one of them follow [instructions:][701] to setup RestVOL, once the required modules are included
* **`h5::kita`** same as above

## Group Operations
#### [Group Creation Property List][1100]

* [`h5::local_heap_size_hint{size_t}`][1101] Specifies the anticipated maximum size of a local heap.
* [`h5::link_creation_order{unsigned}`][1102] Sets creation order tracking and indexing for links in a group.
* [`h5::est_link_info{unsigned, unsigned}`][1103] Sets estimated number of links and length of link names in a group.
* [`h5::link_phase_change{unsigned, unsigned}`][1104] Sets the parameters for conversion between compact and dense groups. 

#### [Group Access Property List][1200]
* [`local_heap_size_hint{hbool_t is_collective}`][1201] Sets metadata I/O mode for read operations to collective or independent (default).


## Link Operations
#### [Link Creation Property List][1300]
* [`h5::char_encoding{H5T_cset_t}`][1301] Sets the character encoding used to encode link and attribute names. <br/>
**Flags:** `h5::utf8`, `h5::ascii`
* [`h5::create_intermediate_group{unsigned}`][1302] Specifies in property list whether to create missing intermediate groups. <br/>
**Flags** `h5::create_path`, `h5::dont_create_path`

#### [Link Access Property List][1400]
* [`h5::nlinks{size_t}`][1401] Sets maximum number of soft or user-defined link traversals.
* [`h5::elink_cb{H5L_elink_traverse_t, void*}`][1402] Sets metadata I/O mode for read operations to collective or independent (default).
* [`h5::elink_fapl{hid_t}`][1403] Sets a file access property list for use in accessing a file pointed to by an external link.
* [`h5::elink_acc_flags{unsigned}`][1403] Sets the external link traversal file access flag in a link access property list.<br/>
  **Flags:** 	`h5::acc_rdwr`, `h5::acc_rdonly`, `h5::acc_default`

## Dataset Operations
#### [Dataset Creation Property List][1500]
Filtering properties require `h5::chunk{..}` set to sensible values. Not having set `h5::chunk{..}` set is equal with requesting `h5::layout_contigous` a dataset layout without chunks, filtering and sub-setting capabilities. This layout is useful for single shot read/write operations, and is the prefered method to save small linear algebra objects.
**Example:**
```cpp
h5::dcpl_t dcpl = h5::chunk{1,4,5} | h5::deflate{4} | h5::compact | h5::dont_filter_partial_chunks
    | h5::fill_value<my_struct>{STR} | h5::fill_time_never | h5::alloc_time_early 
    | h5::fletcher32 | h5::shuffle | h5::nbit;
```
* [`h5::layout{H5D_layout_t layout}`][1501] Sets the type of storage used to store the raw data for a dataset. <br/>
**Flags:** `h5::compact`, `h5::contigous`, `h5::chunked`, `h5::virtual`
* [`h5::chunk{...}`][1502] control chunk size, takes in initializer list with rank matching the dataset dimensions, sets `h5::chunked` layout
* [`h5::chunk_opts{unsigned}`][1503] Sets the edge chunk option in a dataset creation property list.<br/>
**Flags:** `h5::dont_filter_partial_chunks`
* [`h5::deflate{0-9}`][1504] | `h5::gzip{0-9}` set deflate compression ratio
* [`h5::fill_value<T>{T* ptr}`][1505] sets fill value
* [`h5::fill_time{H5D_fill_time_t fill_time}`][1506] Sets the time when fill values are written to a dataset.<br/>
**Flags:** `h5::fill_time_ifset`, `h5::fill_time_alloc`, `h5::fill_time_never`
* [`h5::alloc_time{H5D_alloc_time_t alloc_time}`][1507] Sets the timing for storage space allocation.<br/>
**Flags:** `h5::alloc_time_default`, `h5::alloc_time_early`, `h5::alloc_time_incr`, `h5::alloc_time_late`
* [`h5::fletcher32`][1509] Sets up use of the Fletcher32 checksum filter.
* [`h5::nbit`][1510] Sets up the use of the N-Bit filter.
* [`h5::shuffle`][1512] Sets up use of the shuffle filter.

#### [Dataset Access Property List][1600]
In addition to CAPI properties the follwoing properties are added to provide finer control dataset layouts, and filter chains.

* Mutually exclusive Sparse Matrix Properties to alternate between formats

  * **`h5::csc`** compressed sparse column format
  * **`h5::csr`** compressed sparse row format
  * **`h5::coo`** coordinates list

* **`h5::multi{["field 0", .., "field n"]}`** objects are persisted in multiple datasets with optional field name definition. When specified the system will break up compound classes into basic components and writes each of them into a directory/folder specified with `dataset_name` argument.

* **`h5::high_throughput`** Sets high throughput H5CPP custom filter chain. HDF5 library comes with a complex, feature rich environment to index data-sets by strides, blocks, or even by individual coordinates within chunk boundaries - less fortunate the performance impact on throughput.  Setting this flag will replace the built in filter chain with a simpler one (without complex indexing features), then delegates IO calls to the recently introduced [HDF5 Optimized API][400] calls.<br/>

  The implementation is based on BLAS level 3 blocking algorithm, supports data access only at chunk boundaries, edges are handled as expected. For maximum throughput place edges at chunk boundaries.
  **Note:** This feature and indexing within a chunk boundary such as `h5::stride` is mutually exclusive.

* [`h5::chunk_cache{size_t rdcc_nslots, size_t rdcc_nbytes, double rdcc_w0}`][1601] Sets the raw data chunk cache parameters.
* [`all_coll_metadata_ops{hbool_t is_collective}`][1602] Sets metadata I/O mode for read operations to collective or independent (default).
* [`h5::efile_prefix{const char*}`][1603] Sets the external dataset storage file prefix in the dataset access property list.
* [*`h5::append_flush`*][not-implemented] N/A, H5CPP comes with custom high performance packet table implementation
* [`h5::virtual_view{H5D_vds_view_t}`][1604] Sets the view of the virtual dataset (VDS) to include or exclude missing mapped elements.
* [`h5::virtual_printf_gap{hsize_t}`][1605] Sets the maximum number of missing source files and/or datasets with the printf-style names when getting the extent of an unlimited virtual dataset.

#### [Dataset Transfer Property List][1700]

* [`h5::buffer{size_t, void*, void*}`][1701] Sets type conversion and background buffers.
* [`h5::edc_check{H5Z_EDC_t}`][1702] Sets whether to enable error-detection when reading a dataset.
* [`h5::filter_callback{H5Z_filter_func_t func, void *op_data}`][1703] Sets user-defined filter callback function.
* [`h5::data_transform{const char *expression}`][1704] Sets a data transform expression.
* [`h5::type_conv_cb{H5T_conv_except_func_t func, void *op_data}`][1705] Sets user-defined datatype conversion callback function.
* [`h5::hyper_vector_size{size_t vector_size}`][1706] Sets number of I/O vectors to be read/written in hyperslab I/O.
* [`h5::btree_ratios{double left, double middle, double right}`][1707] Sets B-tree split ratios for a dataset transfer property list.
* [`h5::vlen_mem_manager{H5MM_allocate_t alloc, void *alloc_info, H5MM_free_t free, void *free_info}`][1708]
* [`h5::dxpl_mpiio{H5FD_mpio_xfer_t xfer_mode}`][1709] Sets data transfer mode <br/>
**Flags:** `h5::collective`, `h5::independent`
* [`h5::dxpl_mpiio_chunk_opt{H5FD_mpio_chunk_opt_t opt_mode}`][1710] Sets a flag specifying linked-chunk I/O or multi-chunk I/O.<br/>
**Flags:** `h5::chunk_one_io`, `h5::chunk_multi_io`
* [`h5::dxpl_mpiio_chunk_opt_num{unsigned num_chunk_per_proc}`][1711] Sets a numeric threshold for linked-chunk I/O.
* [`h5::dxpl_mpiio_chunk_opt_ratio{unsigned percent_proc_per_chunk}`][1712] Sets a ratio threshold for collective I/O.
* [`h5::dxpl_mpiio_collective_opt{H5FD_mpio_collective_opt_t opt_mode}`][1713] Sets a flag governing the use of independent versus collective I/O.<br/>
**Flags:** `h5::collective_io`, `h5::individual_io`

## Misc Operations
#### [Object Creation Property List][1800]
* [`h5::ocreate_intermediate_group{unsigned}`][1801] Sets tracking and indexing of attribute creation order.
* [`h5::obj_track_times{hbool_t}`][1802] Sets the recording of times associated with an object.
* [`h5::attr_phase_change{unsigned, unsigned}`][1803] Sets attribute storage phase change thresholds.
* [`h5::attr_creation_order{unsigned}`][1804] Sets tracking and indexing of attribute creation order. <br/>
**Flags:** `h5::crt_order_tracked`, `h5::crt_order_indexed`

#### [Object Copy Property List][1900]
* [`h5::copy_object{unsigned}`][1901] Sets properties to be used when an object is copied. <br/>
**Flags:** `h5::shallow_hierarchy`, `h5::expand_soft_link`, `h5::expand_ext_link`, `h5::expand_reference`, `h5::copy_without_attr`,
`h5::merge_commited_dtype`

## MPI / Parallel Operations
* [`h5::fapl_mpiio{MPI_Comm comm, MPI_Info info}`][1029] Stores MPI IO communicator information to the file access property list.
* [`h5::all_coll_metadata_ops{hbool_t}`][1049] Sets metadata I/O mode for read operations to collective or independent (default).
* [`h5::coll_metadata_write{hbool_t}`][1050] Sets metadata write mode to collective or independent (default).
* [`h5::gc_references{unsigned}`][1051] H5Pset_gc_references sets the flag for garbage collecting references for the file.
* [`h5::small_data_block_size{hsize_t}`][1052] Sets the size of a contiguous block reserved for small data.
* [`h5::object_flush_cb{H5F_flush_cb_t,void*}`][1054] Sets a callback function to invoke when an object flush occurs in the file.
* [`h5::fapl_coll_metadata_ops{hbool_t}`] 
* [`h5::gapl_coll_metadata_ops{hbool_t}`] 
* [`h5::dapl_coll_metadata_ops{hbool_t}`] 
* [`h5::tapl_coll_metadata_ops{hbool_t}`] 
* [`h5::lapl_coll_metadata_ops{hbool_t}`] 
* [`h5::aapl_coll_metadata_ops{hbool_t}`] 
* [`h5::dxpl_mpiio{H5FD_mpio_xfer_t xfer_mode}`][1709] Sets data transfer mode <br/>
**Flags:** `h5::collective`, `h5::independent`
* [`h5::dxpl_mpiio_chunk_opt{H5FD_mpio_chunk_opt_t opt_mode}`][1710] Sets a flag specifying linked-chunk I/O or multi-chunk I/O.<br/>
**Flags:** `h5::chunk_one_io`, `h5::chunk_multi_io`
* [`h5::dxpl_mpiio_chunk_opt_num{unsigned num_chunk_per_proc}`][1711] Sets a numeric threshold for linked-chunk I/O.
* [`h5::dxpl_mpiio_chunk_opt_ratio{unsigned percent_proc_per_chunk}`][1712] Sets a ratio threshold for collective I/O.
* [`h5::dxpl_mpiio_collective_opt{H5FD_mpio_collective_opt_t opt_mode}`][1713] Sets a flag governing the use of independent versus collective I/O.<br/>
**Flags:** `h5::collective_io`, `h5::individual_io`


