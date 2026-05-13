include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libmath.la 
if ENABLE_DEBUG
noinst_LTLIBRARIES += libmath_debug.la
endif

libmath_la_SOURCES = $(h_base_math_sources) $(cpp_base_math_sources)
libmath_la_CXXFLAGS = $(CXXFLAGS_OPT)

libmath_debug_la_SOURCES = $(h_base_math_sources) $(cpp_base_math_sources)
libmath_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/base/math
pkginclude_HEADERS = $(h_base_math_sources)
