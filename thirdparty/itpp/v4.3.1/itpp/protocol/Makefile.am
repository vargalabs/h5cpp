include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libprotocol.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libprotocol_debug.la
endif

libprotocol_la_SOURCES = $(h_protocol_sources) $(cpp_protocol_sources) 
libprotocol_la_CXXFLAGS = $(CXXFLAGS_OPT)

libprotocol_debug_la_SOURCES = $(h_protocol_sources) $(cpp_protocol_sources) 
libprotocol_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/protocol
pkginclude_HEADERS = $(h_protocol_sources)
