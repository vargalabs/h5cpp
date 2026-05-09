### Armadillo: C++ Library for Linear Algebra & Scientific Computing  
https://arma.sourceforge.net

Copyright 2008-2023 Conrad Sanderson (https://conradsanderson.id.au)  
Copyright 2008-2016 National ICT Australia (NICTA)  
Copyright 2017-2023 Data61 / CSIRO  

---

### Quick Links  

- [download latest stable release](https://arma.sourceforge.net/download.html)
- [documentation for functions and classes](https://arma.sourceforge.net/docs.html)
- [bug reports & questions](https://arma.sourceforge.net/faq.html)  

---
 
### Contents

1.  [Introduction](#1-introduction)
2.  [Citation Details](#2-citation-details)
3.  [Distribution License](#3-distribution-license)

4.  [Prerequisites and Dependencies](#4-prerequisites-and-dependencies)

5.  [Linux and macOS: Installation](#5-linux-and-macos-installation)
6.  [Linux and macOS: Compiling and Linking](#6-linux-and-macos-compiling-and-linking)

7.  [Windows: Installation](#7-windows-installation)
8.  [Windows: Compiling and Linking](#8-windows-compiling-and-linking)

9.  [Support for OpenBLAS and Intel MKL](#9-support-for-openblas-and-intel-mkl)
10. [Caveat on use of C++11 auto Keyword](#10-caveat-on-use-of-c11-auto-keyword)
11. [Support for OpenMP](#11-support-for-openmp)

12. [Documentation of Functions and Classes](#12-documentation-of-functions-and-classes)
13. [API Stability and Version Policy](#13-api-stability-and-version-policy)
14. [Bug Reports and Frequently Asked Questions](#14-bug-reports-and-frequently-asked-questions)

15. [MEX Interface to Octave/Matlab](#15-mex-interface-to-octavematlab)
16. [Related Software Using Armadillo](#16-related-software-using-armadillo)

---

### 1. Introduction

Armadillo is a high quality C++ library for linear algebra and scientific computing,
aiming towards a good balance between speed and ease of use.

It's useful for algorithm development directly in C++,
and/or quick conversion of research code into production environments.
It has high-level syntax and functionality which is deliberately similar to Matlab.

The library provides efficient classes for vectors, matrices and cubes,
as well as 200+ associated functions covering essential and advanced functionality
for data processing and manipulation of matrices.

Various matrix decompositions (eigen, SVD, QR, etc) are provided through
integration with LAPACK, or one of its high performance drop-in replacements
(eg. OpenBLAS, Intel MKL, Apple Accelerate framework, etc).

A sophisticated expression evaluator (via C++ template meta-programming)
automatically combines several operations (at compile time) to increase speed
and efficiency.

The library can be used for machine learning, pattern recognition, computer vision,
signal processing, bioinformatics, statistics, finance, etc.

Authors:
  * Conrad Sanderson - https://conradsanderson.id.au
  * Ryan Curtin      - https://ratml.org

---

### 2: Citation Details

Please cite the following papers if you use Armadillo in your research and/or software.  
Citations are useful for the continued development and maintenance of the library.

  * Conrad Sanderson and Ryan Curtin.  
    Armadillo: a template-based C++ library for linear algebra.  
    Journal of Open Source Software, Vol. 1, No. 2, pp. 26, 2016.  
  
  * Conrad Sanderson and Ryan Curtin.  
    A User-Friendly Hybrid Sparse Matrix Class in C++.  
    Lecture Notes in Computer Science (LNCS), Vol. 10931, pp. 422-430, 2018.

---

### 3: Distribution License

Armadillo can be used in both open-source and proprietary (closed-source) software.

Armadillo is licensed under the Apache License, Version 2.0 (the "License").
A copy of the License is included in the "LICENSE.txt" file.

Any software that incorporates or distributes Armadillo in source or binary form
must include, in the documentation and/or other materials provided with the software,
a readable copy of the attribution notices present in the "NOTICE.txt" file.
See the License for details. The contents of the "NOTICE.txt" file are for
informational purposes only and do not modify the License.

---

### 4: Prerequisites and Dependencies

The functionality of Armadillo is partly dependent on other libraries:
- OpenBLAS (or standard BLAS)
- LAPACK
- ARPACK
- SuperLU

Use of OpenBLAS (instead of standard BLAS) is strongly recommended on all systems.
On macOS, the Accelerate framework can be used for BLAS and LAPACK functions.

If sparse matrices are not needed, ARPACK and SuperLU are not required.

Armadillo requires a C++ compiler that supports at least the C++11 standard.

On Linux-based systems, install the GCC C++ compiler, which is available as a pre-built package.
The package name might be `g++` or `gcc-c++` depending on your system.

On macOS systems, a C++ compiler can be obtained by first installing Xcode (at least version 8)
and then running the following command in a terminal window:  

    xcode-select --install

On Windows systems, the MinGW toolset or Visual Studio C++ 2019 (MSVC) can be used.

Caveats on the use of SuperLU:
- SuperLU must be available as a shared library
- Only the following SuperLU versions are supported: 5.2.x, 5.3.x, 6.0.x
- SuperLU 6.0.x must be compiled with default integer size (32 bits)

---

### 5: Linux and macOS: Installation

Armadillo can be installed in several ways: either manually or via cmake, with or without root access.
The cmake based installation is preferred.

The cmake tool can be downloaded from https://www.cmake.org
or (preferably) installed using the package manager on your system;
on macOS systems, cmake can be installed through MacPorts or Homebrew.

Before installing Armadillo, first install OpenBLAS and LAPACK, and optionally ARPACK and SuperLU.
It is also necessary to install the corresponding development files for each library.
For example, when installing the `libopenblas` package, also install the `libopenblas-dev` package.


#### 5a: Installation via CMake

The cmake based installer detects which relevant libraries
are installed on your system (eg. OpenBLAS, LAPACK, SuperLU, ARPACK, etc)
and correspondingly modifies Armadillo's configuration.
The installer also generates the Armadillo runtime library,
which is a wrapper for all the detected libraries.

Change into the directory that was created by unpacking the armadillo archive
(eg. `cd armadillo-10.6.1`) and then run cmake using:

    cmake .

**NOTE:** the full stop (.) separated from `cmake` by a space is important.

On macOS, to enable the detection of OpenBLAS, 
use the additional `ALLOW_OPENBLAS_MACOS` option when running cmake:

    cmake -DALLOW_OPENBLAS_MACOS=ON .

Depending on your installation, OpenBLAS may masquerade as standard BLAS.
To detect standard BLAS and LAPACK, use the `ALLOW_BLAS_LAPACK_MACOS` option:

    cmake -DALLOW_BLAS_LAPACK_MACOS=ON .

By default, cmake assumes that the Armadillo runtime library and the corresponding header files 
will be installed in the default system directory (eg. in the `/usr` hierarchy in Linux-based systems).
To install the library and headers in an alternative directory,
use the additional option `CMAKE_INSTALL_PREFIX` in this form:

    cmake . -DCMAKE_INSTALL_PREFIX:PATH=alternative_directory

If cmake needs to be re-run, it's a good idea to first delete the `CMakeCache.txt` file
(not `CMakeLists.txt`).

**Caveat:** if Armadillo is installed in a non-system directory,
make sure that the C++ compiler is configured to use the `lib` and `include`
sub-directories present within this directory.
Note that the `lib` directory might be named differently on your system.
On recent 64 bit Debian & Ubuntu systems it is `lib/x86_64-linux-gnu`.
On recent 64 bit Fedora & RHEL systems it is `lib64`.

If you have sudo access (ie. root/administrator/superuser privileges)
and didn't use the `CMAKE_INSTALL_PREFIX` option, run the following command:

    sudo make install

If you don't have sudo access, make sure to use the `CMAKE_INSTALL_PREFIX` option
and run the following command:

    make install


#### 5b: Manual Installation

Manual installation involves simply copying the `include/armadillo` header
**and** the associated `include/armadillo_bits` directory to a location
such as `/usr/include/` which is searched by your C++ compiler.
If you don't have sudo access or don't have write access to `/usr/include/`,
use a directory within your own home directory (eg. `/home/user/include/`).

If required, modify `include/armadillo_bits/config.hpp`
to indicate which libraries are currently available on your system.
Comment or uncomment the following lines:

    #define ARMA_USE_LAPACK  
    #define ARMA_USE_BLAS  
    #define ARMA_USE_ARPACK  
    #define ARMA_USE_SUPERLU  

If support for sparse matrices is not needed, ARPACK and SuperLU are not necessary.

Note that the manual installation will not generate the Armadillo runtime library,
and hence you will need to link your programs directly with OpenBLAS, LAPACK, etc.

---

### 6: Linux and macOS: Compiling and Linking

If you have installed Armadillo via the cmake installer,
use the following command to compile your programs:

    g++ prog.cpp -o prog -O2 -std=c++11 -larmadillo

If you have installed Armadillo manually, link with OpenBLAS and LAPACK
instead of the Armadillo runtime library:

    g++ prog.cpp -o prog -O2 -std=c++11 -lopenblas -llapack

If you have manually installed Armadillo in a non-standard location,
such as `/home/user/include/`, you will need to make sure 
that your C++ compiler searches `/home/user/include/` 
by explicitly specifying the directory as an argument/option. 
For example, using the `-I` switch in GCC and Clang:

    g++ prog.cpp -o prog -O2 -std=c++11 -I /home/user/include/ -lopenblas -llapack

If you're getting linking issues (unresolved symbols),
enable the `ARMA_DONT_USE_WRAPPER` option:

    g++ prog.cpp -o prog -O2 -std=c++11 -I /home/user/include/ -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

If you don't have OpenBLAS, on Linux change `-lopenblas` to `-lblas`;
on macOS change `-lopenblas -llapack` to `-framework Accelerate`

The `examples` directory contains a short example program that uses Armadillo.

We recommend that compilation is done with optimisation enabled,
in order to make best use of the extensive template meta-programming
techniques employed in Armadillo.
For GCC and Clang compilers use `-O2` or `-O3` to enable optimisation.

For more information on compiling and linking, see the Questions page: 
https://arma.sourceforge.net/faq.html

---

### 7: Windows: Installation

The installation is comprised of 3 steps:

* Step 1:
  Copy the entire `include` folder to a convenient location
  and tell your compiler to use that location for header files
  (in addition to the locations it uses already).
  Alternatively, the `include` folder can be used directly.

* Step 2:
  If required, modify `include/armadillo_bits/config.hpp`
  to indicate which libraries are currently available on your system:

    #define ARMA_USE_LAPACK  
    #define ARMA_USE_BLAS  
    #define ARMA_USE_ARPACK  
    #define ARMA_USE_SUPERLU  

  If support for sparse matrices is not needed, ARPACK or SuperLU are not necessary.

* Step 3:
  Configure your compiler to link with LAPACK and BLAS
  (and optionally ARPACK and SuperLU).
  Note that OpenBLAS can be used as a high-performance substitute
  for both LAPACK and BLAS.

---

### 8: Windows: Compiling and Linking

Within the `examples` folder, the MSVC project named `example1_win64`
can be used to compile `example1.cpp`.
The project needs to be compiled as a 64 bit program:
the active solution platform must be set to x64, instead of win32.

The MSVC project was tested on Windows 10 (64 bit) with Visual Studio C++ 2019.
Adaptations may be required for 32 bit systems, later versions of Windows and/or the compiler.
For example, options such as `ARMA_BLAS_LONG` and `ARMA_BLAS_UNDERSCORE`,
defined in `include/armadillo_bits/config.hpp`, may need to be either enabled or disabled.

The folder `examples/lib_win64` contains a copy of lib and dll files
obtained from a pre-compiled release of OpenBLAS:
https://github.com/xianyi/OpenBLAS/releases/  
The compilation was done by a third party.  USE AT YOUR OWN RISK.

**Caveat:** 
for any high performance scientific/engineering workloads,
we strongly recommend using a Linux-based operating system, such as:
  * Fedora  https://fedoraproject.org/
  * Ubuntu  https://www.ubuntu.com/
  * CentOS  https://centos.org/

---

### 9: Support for OpenBLAS and Intel MKL

Armadillo can use OpenBLAS or Intel Math Kernel Library (MKL) as high-speed
replacements for BLAS and LAPACK. In essence this involves linking with the
replacement libraries instead of BLAS and LAPACK.

Minor modifications to `include/armadillo_bits/config.hpp` may be required
to ensure Armadillo uses the same integer sizes and style of function names
as used by the replacement libraries. Specifically, the following defines
may need to be enabled or disabled:

    ARMA_USE_WRAPPER  
    ARMA_BLAS_CAPITALS  
    ARMA_BLAS_UNDERSCORE  
    ARMA_BLAS_LONG  
    ARMA_BLAS_LONG_LONG  

See the documentation for more information on the above defines.

On Linux-based systems, MKL might be installed in a non-standard location such as `/opt`
which can cause problems during linking.
Before installing Armadillo, the system should know where the MKL libraries are located.
For example, `/opt/intel/mkl/lib/intel64/`.
This can be achieved by setting the `LD_LIBRARY_PATH` environment variable,
or for a more permanent solution, adding the directory locations to `/etc/ld.so.conf`.
It may also be possible to store a text file with the locations
in the `/etc/ld.so.conf.d` directory. For example, `/etc/ld.so.conf.d/mkl.conf`.
If `/etc/ld.so.conf` is modified or `/etc/ld.so.conf.d/mkl.conf` is created,
`/sbin/ldconfig` must be run afterwards.

Below is an example of `/etc/ld.so.conf.d/mkl.conf`
where Intel MKL is installed in `/opt/intel`

    /opt/intel/lib/intel64  
    /opt/intel/mkl/lib/intel64  

If MKL is installed and it is persistently giving problems during linking,
Support for MKL can be disabled by editing the CMakeLists.txt file,
deleting CMakeCache.txt and re-running the cmake based installation.
Comment out the line containing:

    INCLUDE(ARMA_FindMKL)

---

### 10: Caveat on use of C++11 auto Keyword

Use of the C++11 `auto` keyword is not recommended with Armadillo objects and expressions.

Armadillo has a template meta-programming framework which creates lots of short lived temporaries
that are not properly handled by `auto`.

---

### 11: Support for OpenMP

Armadillo can use OpenMP to automatically speed up computationally
expensive element-wise functions such as exp(), log(), cos(), etc.
This requires a C++ compiler with OpenMP 3.1+ support.

For GCC and Clang compilers, use the following option to enable OpenMP:
`-fopenmp`

---

### 12: Documentation of Functions and Classes

The documentation of Armadillo functions and classes is available at:  
https://arma.sourceforge.net/docs.html

The documentation is also in the `docs.html` file distributed with Armadillo.
Use a web browser to view it.

---

### 13: API Stability and Version Policy

Each release of Armadillo has its public API (functions, classes, constants)
described in the accompanying API documentation (docs.html) specific
to that release.

Each release of Armadillo has its full version specified as A.B.C,
where A is a major version number, B is a minor version number, and C is a patch level.
The version specification has explicit meaning
(similar to [Semantic Versioning](https://semver.org/)), as follows:

* Within a major version (eg. 10), each minor version has a public API that
  strongly strives to be backwards compatible (at the source level) with the
  public API of preceding minor versions. For example, user code written for
  version 10.0 should work with version 10.1, 10.2, etc.
  However, subsequent minor versions may have more features (API additions and extensions)
  than preceding minor versions. As such, user code _specifically_
  written for version 10.2 may not work with 10.1.

* An increase in the patch level, while the major and minor versions are retained,
  indicates modifications to the code and/or documentation which aim to fix bugs
  without altering the public API.

* We don't like changes to existing public API and strongly prefer not to break
  any user software. However, to allow evolution, the public API in future major versions
  while remaining backwards compatible in as many cases as possible
  (eg. major version 11 may have slightly different public API than major version 10).

**CAVEAT:**
the above policy applies only to the public API described in the documentation.
Any functionality within Armadillo which is _not explicitly_ described
in the public API documentation is considered as internal implementation detail,
and may be changed or removed without notice.

---

### 14: Bug Reports and Frequently Asked Questions

Armadillo has gone through extensive testing and has been successfully
used in production environments. However, as with almost all software,
it's impossible to guarantee 100% correct functionality.

If you find a bug in the library or the documentation, we are interested
in hearing about it. Please make a _small_ and _self-contained_ program
which exposes the bug, and then send the program source and the bug description
to the developers. The small program must have a main() function and use only
functions/classes from Armadillo and the standard C++ library (no other libraries).

The contact details are at:  
https://arma.sourceforge.net/contact.html

Further information about Armadillo is on the frequently asked questions page:  
https://arma.sourceforge.net/faq.html

---

### 15: MEX Interface to Octave/Matlab

The `mex_interface` folder contains examples of how to interface
Octave/Matlab with C++ code that uses Armadillo matrices.

---

### 16: Related Software Using Armadillo

* ensmallen: C++ library for non-linear numerical optimisation (L-BFGS, SGD, CMA-ES, etc)  
  https://ensmallen.org/

* Bandicoot: C++ library for accelerated linear algebra on GPUs  
  https://coot.sourceforge.io

* MLPACK: extensive library of machine learning algorithms  
  https://mlpack.org

* RcppArmadillo: integration of Armadillo with R  
  https://dirk.eddelbuettel.com/code/rcpp.armadillo.html

* CARMA: interface between Armadillo and Python / NumPy  
  https://github.com/RUrlus/carma

