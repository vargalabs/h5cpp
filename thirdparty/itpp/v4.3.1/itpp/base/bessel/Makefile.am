include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libbessel.la 
if ENABLE_DEBUG
noinst_LTLIBRARIES += libbessel_debug.la
endif

libbessel_la_SOURCES = $(noinst_h_base_bessel_sources) \
	$(cpp_base_bessel_sources)
libbessel_la_CXXFLAGS = $(CXXFLAGS_OPT)

libbessel_debug_la_SOURCES = $(noinst_h_base_bessel_sources) \
	$(cpp_base_bessel_sources)
libbessel_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

