dnl @synopsis ACX_FFT
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2007-02-15
dnl @license GPLWithACException
dnl
dnl This macro looks for some FFT implementation, e.g. FFTW3, MKL, ACML,
dnl etc. On success, it sets the FFT_LIBS output variable to hold the
dnl requisite library linkages. Besides it defines one of the following:
dnl HAVE_FFT, HAVE_FFTW, HAVE_MKL_FFT, HAVE_ACML_FFT
dnl
dnl The user may also use --with-fft=<lib> in order to use some
dnl specific FFT library <lib>.

AC_DEFUN([ACX_FFT], [

# Initialise local variables
acx_fft_ok=no
fft_mkl_ok=no
fft_acml_ok=no
fftw3_ok=no


# Parse "--with-fft=<lib>" option
AC_ARG_WITH(fft, [AS_HELP_STRING([--with-fft@<:@=LIB@:>@], [use FFT library, optionally specified by LIB])])
case $with_fft in
  yes | "") ;;
  no) acx_fft_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) FFT_LIBS="$with_fft" ;;
  *) FFT_LIBS="-l$with_fft" ;;
esac

# Parse "--with-fft-include=<path>" option
AC_ARG_WITH(fft_include, [AS_HELP_STRING([--with-fft-include=DIR],
    [path to FFT header files])],
  [CPPFLAGS="$CPPFLAGS -I$with_fft_include"])

# First, check FFT_LIBS environment variable
if test "x$FFT_LIBS" != x; then
  save_LIBS="$LIBS"; LIBS="$FFT_LIBS $LIBS"
  AC_MSG_CHECKING([for DftiComputeForward in $FFT_LIBS])
  AC_TRY_LINK_FUNC(DftiComputeForward, [acx_fft_ok=yes])
  AC_MSG_RESULT($acx_fft_ok)
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([mkl_dfti.h], [fft_mkl_ok=yes], [acx_fft_ok=no])
  fi
  if test "$acx_fft_ok" = no; then
    AC_MSG_CHECKING([for zfft1dx in $FFT_LIBS])
    AC_TRY_LINK_FUNC(zfft1dx, [acx_fft_ok=yes])
    AC_MSG_RESULT($acx_fft_ok)
    if test "$acx_fft_ok" = yes; then
      AC_CHECK_HEADER([acml.h], [fft_acml_ok=yes], [acx_fft_ok=no])
    fi
  fi
  if test "$acx_fft_ok" = no; then
    AC_MSG_CHECKING([for fftw_plan_dft_1d in $FFT_LIBS])
    AC_TRY_LINK_FUNC(fftw_plan_dft_1d, [acx_fft_ok=yes], [FFT_LIBS=""])
    AC_MSG_RESULT($acx_fft_ok)
    if test "$acx_fft_ok" = yes; then
      AC_CHECK_HEADER([fftw3.h], [fftw3_ok=yes], [acx_fft_ok=no; FFT_LIBS=""])
    fi
  fi
  LIBS="$save_LIBS"
fi

# FFT in BLAS (MKL) library?
if test "$acx_fft_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC(DftiComputeForward, [acx_fft_ok=yes])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([mkl_dfti.h], [fft_mkl_ok=yes], [acx_fft_ok=no])
  fi
  LIBS="$save_LIBS"
fi

# FFT in BLAS (ACML) library?
if test "$acx_fft_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC(zfft1dx, [acx_fft_ok=yes])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([acml.h], [fft_acml_ok=yes], [acx_fft_ok=no])
  fi
  LIBS="$save_LIBS"
fi

# FFT in FFTW3 library?
if test "$acx_fft_ok" = no; then
  AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, [acx_fft_ok=yes])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([fftw3.h], [fftw3_ok=yes; FFT_LIBS="-lfftw3"],
      [acx_fft_ok=no])
  fi
fi

# FFT in FFTW3 library (extra -lm)?
if test "$acx_fft_ok" = no; then
  AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, [acx_fft_ok=yes], [], [-lm])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([fftw3.h], [fftw3_ok=yes; FFT_LIBS="-lfftw3 -lm"],
      [acx_fft_ok=no])
  fi
fi

# FFT in MKL library?
if test "$acx_fft_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS"
  AC_CHECK_LIB(mkl, DftiComputeForward, [acx_fft_ok=yes], [],
    [-lguide -lpthread])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([mkl_dfti.h],
      [fft_mkl_ok=yes; FFT_LIBS="-lmkl -lguide -lpthread"], [acx_fft_ok=no])
  fi
  LIBS="$save_LIBS"
fi

# FFT in ACML library?
if test "$acx_fft_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS$MY_FLIBS"
  AC_CHECK_LIB(acml, zfft1dx, [acx_fft_ok=yes])
  if test "$acx_fft_ok" = yes; then
    AC_CHECK_HEADER([acml.h],
      [fft_acml_ok=yes; FFT_LIBS="-lacml$MY_FLIBS"], [acx_fft_ok=no])
  fi
  LIBS="$save_LIBS"
fi

AC_SUBST(FFT_LIBS)

# Finally, define HAVE_*
if test "$acx_fft_ok" = yes; then
  AC_DEFINE(HAVE_FFT, 1, [Define if you have FFT library.])
  if test "$fft_mkl_ok" = yes; then
    AC_DEFINE(HAVE_FFT_MKL, 1, [Define if you have MKL FFT library.])
  fi
  if test "$fft_acml_ok" = yes; then
    AC_DEFINE(HAVE_FFT_ACML, 1, [Define if you have ACML FFT library.])
  fi
  if test "$fftw3_ok" = yes; then
    AC_DEFINE(HAVE_FFTW3, 1, [Define if you have FFTW3 library.])
  fi
else
  if test "$acx_fft_ok" != disabled; then
    AC_MSG_ERROR([cannot find any FFT library.
You can override this error by using "--without-fft" option, but the
functionality of the IT++ library will be limited. You have been warned!])
  fi
fi

])dnl ACX_FFT
