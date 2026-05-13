include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libsignal.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libsignal_debug.la
endif

libsignal_la_SOURCES = $(h_signal_sources) $(cpp_signal_sources) 
libsignal_la_CXXFLAGS = $(CXXFLAGS_OPT)

libsignal_debug_la_SOURCES = $(h_signal_sources) $(cpp_signal_sources) 
libsignal_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/signal
pkginclude_HEADERS = $(h_signal_sources)
