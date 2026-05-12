include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libstat.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libstat_debug.la
endif

libstat_la_SOURCES = $(h_stat_sources) $(cpp_stat_sources) 
libstat_la_CXXFLAGS = $(CXXFLAGS_OPT)

libstat_debug_la_SOURCES = $(h_stat_sources) $(cpp_stat_sources) 
libstat_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/stat
pkginclude_HEADERS = $(h_stat_sources)
