include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libcomm.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libcomm_debug.la
endif

libcomm_la_SOURCES = $(h_comm_sources) $(cpp_comm_sources) 
libcomm_la_CXXFLAGS = $(CXXFLAGS_OPT)

libcomm_debug_la_SOURCES = $(h_comm_sources) $(cpp_comm_sources) 
libcomm_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/comm
pkginclude_HEADERS = $(h_comm_sources)
