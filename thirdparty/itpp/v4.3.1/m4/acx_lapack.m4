dnl @synopsis ACX_LAPACK
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2007-02-15
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/). On
dnl success, it sets the LAPACK_LIBS output variable to hold the
dnl requisite library linkages. Besides, it defines HAVE_LAPACK
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS
dnl
dnl in that order. BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>. In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as was
dnl used to compile the LAPACK and BLAS libraries.

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])

test "x$cheev" = x && cheev=cheev_

# Initialise local variables
# We cannot use LAPACK if BLAS is not found
if test "$acx_blas_ok" != yes; then
  acx_lapack_ok=noblas
else
  acx_lapack_ok=no
fi

# Parse "--with-lapack=<lib>" option
AC_ARG_WITH(lapack,
  [AS_HELP_STRING([--with-lapack@<:@=LIB@:>@], [use LAPACK library, optionally specified by LIB])])
case $with_lapack in
  yes | "") ;;
  no) acx_lapack_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
  *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
  save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
  AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
  AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
  AC_MSG_RESULT($acx_lapack_ok)
  LIBS="$save_LIBS"
fi

# LAPACK linked to by default?  (it is sometimes included in BLAS)
if test "$acx_lapack_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
  LIBS="$save_LIBS"
fi

# LAPACK in MKL library?
# (http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm)
if test "$acx_lapack_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  AC_CHECK_LIB(mkl_lapack32, $cheev,
    [acx_lapack_ok=yes; LAPACK_LIBS="-lmkl_lapack32 -lmkl_lapack64"],
    [AC_CHECK_LIB(mkl_lapack, $cheev,
      [acx_lapack_ok=yes; LAPACK_LIBS="-lmkl_lapack"])],
    [-lmkl_lapack64])
  LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
  if test "$acx_lapack_ok" = no; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
    AC_CHECK_LIB($lapack, $cheev, [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"])
    LIBS="$save_LIBS"
  fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, define HAVE_LAPACK
if test "$acx_lapack_ok" = yes; then
  AC_DEFINE(HAVE_LAPACK, 1, [Define if you have LAPACK library.])
else
  if test "$acx_lapack_ok" != disabled; then
    AC_MSG_ERROR([cannot find any LAPACK library.
You can override this error by using "--without-lapack" option, but the
functionality of the IT++ library will be limited. You have been warned!])
  fi
fi

])dnl ACX_LAPACK
