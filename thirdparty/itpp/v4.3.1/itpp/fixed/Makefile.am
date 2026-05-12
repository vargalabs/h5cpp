include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libfixed.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libfixed_debug.la
endif

libfixed_la_SOURCES = $(h_fixed_sources) $(cpp_fixed_sources) 
libfixed_la_CXXFLAGS = $(CXXFLAGS_OPT)

libfixed_debug_la_SOURCES = $(h_fixed_sources) $(cpp_fixed_sources) 
libfixed_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/fixed
pkginclude_HEADERS = $(h_fixed_sources)
