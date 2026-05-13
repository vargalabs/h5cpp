include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

if ENABLE_COMM
  COMM_DIR = comm
  COMM_LIB = $(top_builddir)/@PACKAGE@/comm/libcomm.la
  COMM_D_LIB = $(top_builddir)/@PACKAGE@/comm/libcomm_debug.la
endif

if ENABLE_FIXED
  FIXED_DIR = fixed
  FIXED_LIB = $(top_builddir)/@PACKAGE@/fixed/libfixed.la
  FIXED_D_LIB = $(top_builddir)/@PACKAGE@/fixed/libfixed_debug.la
endif

if ENABLE_OPTIM
  OPTIM_DIR = optim
  OPTIM_LIB = $(top_builddir)/@PACKAGE@/optim/liboptim.la
  OPTIM_D_LIB = $(top_builddir)/@PACKAGE@/optim/liboptim_debug.la
endif

if ENABLE_PROTOCOL
  PROTOCOL_DIR = protocol
  PROTOCOL_LIB = $(top_builddir)/@PACKAGE@/protocol/libprotocol.la
  PROTOCOL_D_LIB = $(top_builddir)/@PACKAGE@/protocol/libprotocol_debug.la
endif

if ENABLE_SIGNAL
  SIGNAL_DIR = signal
  SIGNAL_LIB = $(top_builddir)/@PACKAGE@/signal/libsignal.la
  SIGNAL_D_LIB = $(top_builddir)/@PACKAGE@/signal/libsignal_debug.la
endif

if ENABLE_SRCCODE
  SRCCODE_DIR = srccode
  SRCCODE_LIB = $(top_builddir)/@PACKAGE@/srccode/libsrccode.la
  SRCCODE_D_LIB = $(top_builddir)/@PACKAGE@/srccode/libsrccode_debug.la
endif

SUBDIRS = base stat $(COMM_DIR) $(FIXED_DIR) $(OPTIM_DIR) $(PROTOCOL_DIR) \
	$(SIGNAL_DIR) $(SRCCODE_DIR)

EXTRA_DIST = config_msvc.h

lib_LTLIBRARIES = lib@PACKAGE@.la
if ENABLE_DEBUG
  lib_LTLIBRARIES += lib@PACKAGE@_debug.la
endif

# Dummy C++ source to cause C++ linking
nodist_EXTRA_lib@PACKAGE@_la_SOURCES = dummy.cxx
lib@PACKAGE@_la_SOURCES =
lib@PACKAGE@_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION) \
	$(NO_UNDEFINED)
lib@PACKAGE@_la_LIBADD = $(top_builddir)/@PACKAGE@/base/libbase.la \
	$(top_builddir)/@PACKAGE@/stat/libstat.la \
	$(COMM_LIB) $(FIXED_LIB) $(OPTIM_LIB) $(PROTOCOL_LIB) $(SIGNAL_LIB) \
	$(SRCCODE_LIB)

# Dummy C++ source to cause C++ linking
nodist_EXTRA_lib@PACKAGE@_debug_la_SOURCES = dummy.cxx
lib@PACKAGE@_debug_la_SOURCES =
lib@PACKAGE@_debug_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION) \
	$(NO_UNDEFINED)
lib@PACKAGE@_debug_la_LIBADD = $(top_builddir)/@PACKAGE@/base/libbase_debug.la \
	$(top_builddir)/@PACKAGE@/stat/libstat_debug.la \
	$(COMM_D_LIB) $(FIXED_D_LIB) $(OPTIM_D_LIB) $(PROTOCOL_D_LIB) \
	$(SIGNAL_D_LIB) $(SRCCODE_D_LIB)

nobase_pkginclude_HEADERS = $(h_sources)
