---
hide:
  - toc
---

H5CPP is a modern C++ interface for HDF5 that lets developers persist vectors, matrices, tensors, strings, and structs with concise template-based code, while preserving full interoperability with the native HDF5 C API. It works naturally with the C++ type system and, through LLVM-based compiler-assisted reflection, reduces persistence of even complex structs to a single-line expression.

## Problem

In quantitative finance, data storage is not a side concern — it sits directly on the critical path of research, simulation, market-data engineering, and production trading systems. A large fraction of quantitative developers come from physics, mathematics, engineering, and other scientific disciplines, where **HDF5 is already familiar ground**. It underpins data workflows across **SciPy**, **MATLAB**, **Julia**, and many other numerical environments, and fits naturally beside optimization and modeling tools such as **Convex.jl**, **JuMP.jl**, and **CVXPy**. That makes HDF5 a particularly attractive foundation for quant workflows: talented people already know the model, trust the format, and can move comfortably between research and production environments.

The problem is that, in low-level C and C++, working with HDF5 directly often turns routine persistence into unnecessary plumbing:

- :material-vector-square:{.icon} manual datatype and dataspace construction
- :material-ruler-square:{.icon} repetitive shape bookkeeping
- :material-code-braces:{.icon} verbose read/write setup
- :material-lifebuoy:{.icon} resource lifetime management
- :material-call-split:{.icon} impedance mismatch between HDF5 and modern C++ container types
- :material-lan-connect:{.icon} friction when moving from local workflows to MPI-enabled parallel I/O

For day-to-day development, the issue is not that HDF5 is weak — it is that too much developer effort goes into plumbing rather than data.

## Impact

That friction shows up everywhere in trading and quantitative research:

- :material-file-code-outline:{.icon} structured event data such as trades, quotes, fills, and order updates becomes harder to persist cleanly
- :material-grid:{.icon} dense homogeneous datasets such as price, volume, returns, feature, and signal matrices require too much ceremony
- :material-chart-line:{.icon} simulation, backtesting, and numerical research code becomes harder to maintain
- :material-content-save-edit:{.icon} checkpoints, experiment outputs, and intermediate research artifacts accumulate boilerplate
- :material-archive-outline:{.icon} historical datasets become harder to archive once and reuse across research, replay, and production workflows
- :material-source-branch-sync:{.icon} codebases that need both serial and parallel HDF5 accumulate complexity quickly

The result is slower iteration, less readable code, and more opportunities for bugs in exactly the places where persistence should be boring, reliable, and fast.

## Approach

H5CPP addresses this by giving HDF5 a **modern C++ interface** centered on types, containers, and RAII. It is particularly well suited to the two dominant data shapes in trading systems and quant workflows:

- :material-file-tree:{.icon} **structured datasets** for trades, quotes, fills, order events, and other record-oriented data
- :material-table-large:{.icon} **homogeneous dense datasets** for returns, features, factors, covariances, signals, and other matrix-like numerical data

H5CPP gives trading-system developers and quantitative analysts a persistence layer that is:

- :material-lightning-bolt:{.icon} **fast to integrate into low-latency and research code**
- :material-file-code-outline:{.icon} **natural for structured event data**
- :material-grid:{.icon} **well suited for dense homogeneous numerical datasets**
- :material-language-cpp:{.icon} **aligned with modern C++ container and type systems**
- :material-shield-check:{.icon} **safer through RAII and reduced descriptor boilerplate**
- :material-lan-connect:{.icon} **capable of scaling from local research workflows to MPI-enabled batch processing**

## Outcome

In practice, that means:

- :material-database-arrow-down:{.icon} market-data handlers can store **structured records** directly as appendable HDF5 datasets
- :material-chart-box:{.icon} quantitative workflows can persist and reload **homogeneous matrices** for returns, features, factors, and signals with minimal ceremony
- :material-clock-fast:{.icon} less engineering time is spent wiring up storage, and more time goes into execution logic, modeling, and analysis
- :material-source-merge:{.icon} research artifacts move more naturally between **C++**, **Python**, **Julia**, **MATLAB**, and **R**
- :material-archive-outline:{.icon} historical datasets can be stored once in a compact, structured format and reused across backtests, simulations, and live strategy development
- :material-rocket-launch-outline:{.icon} teams get a smoother path from prototype to production, especially when the same datasets feed both strategy research and trading infrastructure

- :material-vector-polyline:{.icon} dlib matrix/vector types


```cpp
#include <armadillo>
#include <h5cpp/all>
...
auto fd = h5::open("some_file.h5", H5F_ACC_RDWR);
/* the RVO arma::Mat<double> object will have the size 10x5 filled*/
try {
	/* will drop extents of unit dimension returns a 2D object */
	auto M = h5::read<arma::mat>(fd,"path/to/object", 
        h5::offset{3,4,1}, h5::count{10,1,5}, h5::stride{3,1,1} ,h5::block{2,1,1} );
} catch (const std::runtime_error& ex ){
	...
} // fd closes underlying resource (raii idiom) auto result = h5::read<Eigen::MatrixXd>(fd, "/result");
```

### Structured append workflows

H5CPP supports appendable chunked datasets for structured records:

```cpp
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

## Background

H5CPP has been presented in HDF5 and C++ community venues over multiple years, including HUG sessions, HDF Group events, C++ community talks, and ISC-related material. Topics include compiler-assisted reflection, POD introspection, MPI/parallel I/O, throughput/latency trade-offs, and practical HDF5 workflows.

[:material-presentation-play:{.icon} **Presentations Portfolio**](https://steven-varga.ca/site/talks/)

