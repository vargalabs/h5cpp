include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = liboptim.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += liboptim_debug.la
endif

liboptim_la_SOURCES = $(h_optim_sources) $(cpp_optim_sources) 
liboptim_la_CXXFLAGS = $(CXXFLAGS_OPT)

liboptim_debug_la_SOURCES = $(h_optim_sources) $(cpp_optim_sources) 
liboptim_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/optim
pkginclude_HEADERS = $(h_optim_sources)
