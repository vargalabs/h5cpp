include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS = dist-bzip2

if HTML_DOC
  HTML_DOC_DIR = doc
endif

if RUN_TESTS
  TEST_DIR = tests
endif

SUBDIRS = itpp extras win32 $(HTML_DOC_DIR) $(TEST_DIR)

EXTRA_DIST = \
	autogen.sh \
	ChangeLog-2005 \
	ChangeLog-2006 \
	ChangeLog-2007 \
	ChangeLog-2008 \
	configure.ac.in \
	itpp-config.1.in \
	itpp.spec \
	itpp.spec.in \
	Makefile.common \
	NEWS-3.10 \
	NEWS-3.9 \
	NEWS-3.99 \
	NEWS-4.0 \
	TODO \
	VERSION

dist_man_MANS = itpp-config.1

bin_SCRIPTS = itpp-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = itpp.pc
if ENABLE_DEBUG
  pkgconfig_DATA += itpp_debug.pc
endif

clean-local:
	rm -f config_summary.log

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/itpp

dist-html:
	$(MAKE) -C doc dist-html

