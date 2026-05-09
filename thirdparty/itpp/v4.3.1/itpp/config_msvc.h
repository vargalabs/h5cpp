/* itpp/config.h.  Generated from config.h.in by configure.  */
/* itpp/config.h.in.  Generated from configure.ac by autoheader.  */


#ifndef CONFIG_H
#define CONFIG_H


/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the `acosh' function. */
/* #undef HAVE_ACOSH */

/* Define to 1 if you have the `asinh' function. */
/* #undef HAVE_ASINH */

/* Define to 1 if you have the `atanh' function. */
/* #undef HAVE_ATANH */

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have a BLAS library. */
#  define HAVE_BLAS 1
#endif

#ifdef HAVE_ACML
/* Define if you have an ACML BLAS library. */
#  define HAVE_BLAS_ACML 1
#endif

/* Define if you have an ATLAS BLAS library. */
/* #undef HAVE_BLAS_ATLAS */

#ifdef HAVE_MKL
/* Define if you have an MKL BLAS library. */
#  define HAVE_BLAS_MKL 1
#endif

/* Define to 1 if you have the `cbrt' function. */
/* #undef HAVE_CBRT */

/* Define to 1 if you have the <cmath> header file. */
#define HAVE_CMATH 1

/* Define to 1 if you have the <complex> header file. */
#define HAVE_COMPLEX 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#define HAVE_DECL_ISFINITE 0

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
   */
#define HAVE_DECL_ISINF 0

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
   */
#define HAVE_DECL_ISNAN 0

/* Define to 1 if you have the declaration of `signgam', and to 0 if you
   don't. */
#define HAVE_DECL_SIGNGAM 0

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `erf' function. */
/* #undef HAVE_ERF */

/* Define to 1 if you have the `erfc' function. */
/* #undef HAVE_ERFC */

/* Define to 1 if you have the `expm1' function. */
/* #undef HAVE_EXPM1 */

/* Define if the compiler supports extern template */
/* #undef HAVE_EXTERN_TEMPLATE */

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have FFT library. */
#  define HAVE_FFT 1
#endif

/* Define if you have FFTW3 library. */
/* #undef HAVE_FFTW3 */

#ifdef HAVE_ACML
/* Define if you have ACML FFT library. */
#  define HAVE_FFT_ACML 1
#endif

#ifdef HAVE_MKL
/* Define if you have MKL FFT library. */
#  define HAVE_FFT_MKL 1
#endif

/* Define to 1 if you have the `finite' function. */
#define HAVE_FINITE 1

/* Define to 1 if you have the `fpclass' function. */
#define HAVE_FPCLASS 1

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* Define to 1 if you have the `isfinite' function. */
#define HAVE_ISFINITE 1

/* Define to 1 if you have the `isinf' function. */
/* #undef HAVE_ISINF */

/* Define to 1 if you have the `isnan' function. */
#define HAVE_ISNAN 1

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have LAPACK library. */
#  define HAVE_LAPACK 1
#endif

/* Define to 1 if you have the `lgamma' function. */
/* #undef HAVE_LGAMMA */

/* Define to 1 if you have the `log1p' function. */
/* #undef HAVE_LOG1P */

/* Define to 1 if you have the `log2' function. */
/* #undef HAVE_LOG2 */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `rint' function. */
/* #undef HAVE_RINT */

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `std::isfinite' function. */
/* #undef HAVE_STD_ISFINITE */

/* Define to 1 if you have the `std::isinf' function. */
/* #undef HAVE_STD_ISINF */

/* Define to 1 if you have the `std::isnan' function. */
/* #undef HAVE_STD_ISNAN */

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `tgamma' function. */
/* #undef HAVE_TGAMMA */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if you use zdotusub_ Fortran wrapper. */
/* #undef HAVE_ZDOTUSUB */

/* Define if "void zdotu_()" should be used. */
#define HAVE_ZDOTU_VOID 1

/* Define if you want exceptions handling */
/* #undef ITPP_EXCEPTIONS */

/* Name of package */
#define PACKAGE "itpp"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ediap@users.sourceforge.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "IT++"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "IT++ 4.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "itpp"

/* Define to the version of this package. */
#define PACKAGE_VERSION "4.0.0"

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 4

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 4

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* The size of `unsigned short', as computed by sizeof. */
#define SIZEOF_UNSIGNED_SHORT 2

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
/* #undef TIME_WITH_SYS_TIME */

/* Version number of package */
#define VERSION "4.0.0"


#if defined(HAVE_CMATH)
#  include <cmath>
#endif

#endif /* #ifndef CONFIG_H */

