include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

noinst_LTLIBRARIES = libalgebra.la 
if ENABLE_DEBUG
noinst_LTLIBRARIES += libalgebra_debug.la
endif

libalgebra_la_SOURCES = $(noinst_h_base_algebra_sources) \
	$(h_base_algebra_sources) \
	$(cpp_base_algebra_sources)
libalgebra_la_CXXFLAGS = $(CXXFLAGS_OPT)

libalgebra_debug_la_SOURCES = $(noinst_h_base_algebra_sources) \
	$(h_base_algebra_sources) \
	$(cpp_base_algebra_sources)
libalgebra_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)

pkgincludedir = $(includedir)/@PACKAGE@/base/algebra
pkginclude_HEADERS = $(h_base_algebra_sources)
