include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libsrccode.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libsrccode_debug.la
endif

libsrccode_la_SOURCES = $(h_srccode_sources) $(cpp_srccode_sources) 
libsrccode_la_CXXFLAGS = $(CXXFLAGS_OPT)

libsrccode_debug_la_SOURCES = $(h_srccode_sources) $(cpp_srccode_sources) 
libsrccode_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/srccode
pkginclude_HEADERS = $(h_srccode_sources)
