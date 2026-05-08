---
hide:
  - toc
---

# Supported Types

H5CPP maps C++ types to HDF5 datasets at compile time. The type dispatch is entirely template-based — no runtime type registry, no virtual calls. Before looking at what is supported, it helps to understand the two core traits that drive dispatch.

## Type Traits

### `h5::is_contiguous<T>`

Evaluates to `std::true_type` when `T` exposes a contiguous memory block that H5CPP can hand directly to the HDF5 C API — that is, when a pointer to the first element plus a byte count is sufficient to describe the object's data.

| Result | Examples |
|--------|---------|
| `true`  | `std::vector<T>`, `std::array<T,N>`, all supported linalg types |
| `false` | `std::list<T>`, `std::map<K,V>`, `std::vector<bool>` |

`std::vector<bool>` is `false` because the standard mandates bit-packing for that specialization — there is no `bool*` to pass to HDF5.

### `h5::storage_representation_t`

A compile-time enum that classifies every recognised type into one of four storage strategies:

| Value | Meaning |
|-------|---------|
| `scalar` | a single arithmetic value or POD struct |
| `linear_value_dataset` | a 1-D sequence of scalars (contiguous) |
| `jagged_value_dataset` | a sequence accessed through iterators |
| `unsupported` | type cannot be stored (e.g. `std::vector<bool>`, raw function pointers) |

## Scalar and Arithmetic Types

All C++ arithmetic types map to their HDF5 NATIVE equivalents:

```cpp
bool                                     → H5T_NATIVE_HBOOL
char / unsigned char                     → H5T_NATIVE_CHAR / H5T_NATIVE_UCHAR
short / unsigned short                   → H5T_NATIVE_SHORT / H5T_NATIVE_USHORT
int / unsigned int                       → H5T_NATIVE_INT / H5T_NATIVE_UINT
long / unsigned long                     → H5T_NATIVE_LONG / H5T_NATIVE_ULONG
long long / unsigned long long           → H5T_NATIVE_LLONG / H5T_NATIVE_ULLONG
float                                    → H5T_NATIVE_FLOAT
double                                   → H5T_NATIVE_DOUBLE
long double                              → H5T_NATIVE_LDOUBLE
char* / const char*                      → H5T_C_S1 (variable-length UTF-8)
```

## STL Containers

H5CPP identifies containers via the [detection idiom (WG21 N4421)](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4436.pdf): any class that provides the expected begin/end/size interface is treated as a container, regardless of whether it is literally from `<vector>` or a user-defined type.

### Contiguous storage — direct I/O

These containers expose `.data()` and allow zero-copy transfers:

| Container | Notes |
|-----------|-------|
| `std::vector<T>` | `T` must be arithmetic or a registered POD struct |
| `std::array<T,N>` | fixed-size; shape inferred at compile time |

```cpp
std::vector<double> v(1024);
h5::write(fd, "signal", v);

std::array<float, 32> coeffs;
h5::read(fd, "filter/coeffs", coeffs);
```

**`std::vector<bool>` is not supported.** The C++ standard mandates bit-packing for this specialization; there is no contiguous `bool*`, so H5CPP classifies it as `unsupported`. Use `std::vector<uint8_t>` or `std::vector<int>` instead.

### Iterator-based storage — chunked I/O

These containers have no `.data()` accessor. H5CPP moves data through iterators in chunk-sized blocks. Memory overhead is bounded by `element_size × (container_size + chunk_size)`.

```cpp
std::list<T>               std::forward_list<T>
std::deque<T>              std::set<T>
std::multiset<T>           std::unordered_set<T>
std::map<K,V>              std::multimap<K,V>
std::unordered_map<K,V>    std::unordered_multimap<K,V>
```

### Container adaptors

`std::stack<T>`, `std::queue<T>`, and `std::priority_queue<T>` are adaptors. Their storage strategy is determined by the underlying container.

### Nested containers

`std::vector<std::vector<T>>` and similar nesting maps to HDF5 ragged/fractal arrays.

## Strings

Both `std::string` and C strings (`char*`, `const char*`) are supported. HDF5 variable-length strings are used, with UTF-8 encoding. Fixed-length strings and `char[N]` buffers require a custom `h5::dt_t<T>`.

```cpp
std::vector<std::string> labels = {"open", "high", "low", "close"};
h5::write(fd, "columns", labels);
```

**Not supported:** `wchar_t`, `char16_t`, `char32_t`.

## Linear Algebra Types

See [Linalg & STL](architecture-integration.md) for the full list of supported Armadillo, Eigen, Blaze, Blitz++, IT++, dlib, and uBLAS types, together with the memory layout conventions used for zero-copy I/O.

## POD Structs and Compound Types

Arbitrarily nested Plain Old Data structs are supported either through the [h5cpp compiler](compiler.md) (which generates the HDF5 datatype descriptor automatically) or by manually calling `H5CPP_REGISTER_STRUCT(Type)`.

```cpp
H5CPP_REGISTER_STRUCT(my_trade_record);

std::vector<my_trade_record> trades = load_trades();
h5::write(fd, "trades/2024-01-05", trades);
```

## Custom Datatypes

Pass an explicit `h5::dt_t<T>` to any I/O call to override the default type mapping. This is useful for fixed-length strings, opaque blobs, and non-standard binary layouts:

```cpp
using fixed_str = char[42];
h5::dt_t<fixed_str> dt{H5Tcreate(H5T_STRING, sizeof(fixed_str))};
H5Tset_cset(dt, H5T_CSET_UTF8);

h5::ds_t ds = h5::create<fixed_str>(fd, "labels",
    h5::current_dims{1000}, h5::chunk{256}, dt);
```

## Raw Pointers

When a type is not directly supported, you can fall back to a raw pointer plus an explicit `h5::count`:

```cpp
my_object obj(100);
h5::read("file.h5", "dataset", obj.ptr(), h5::count{10, 10}, h5::offset{100, 0});
```

`T*` must point to a contiguous block of an arithmetic type or a registered POD struct.
