#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

PREFIX = /usr/local

VERSION = 1.10.4.6
DEB = -1~exp1
PROGRAM_NAME = libh5cpp-dev
BIN_DIR = $(PREFIX)/bin
INCLUDE_DIR = $(PREFIX)/include
EXAMPLE_DIR = $(PREFIX)/share

MAN_BASE_DIR = $(PREFIX)/man
MAN_DIR = $(MAN_BASE_DIR)/man1
MAN_EXT = 1

INSTALL = install
RM = rm -f
MKDIR = mkdir -p

CC ?= gcc
COMPILER_OPTIONS = -Wall -O -g

INSTALL_PROGRAM = $(INSTALL) -c -m 0755
INSTALL_DATA = $(INSTALL) -c -m 0644
INSTALL_EXAMPLE = $(INSTALL) -c -m 0755
INSTALL_INCLUDE = $(INSTALL) -m 0755

OBJECT_FILES =

# ------------------------------------------------------------------
# documentation / mkdocs
# ------------------------------------------------------------------
PYTHON ?= python3
PIP ?= $(PYTHON) -m pip
MKDOCS ?= mkdocs
MKDOCS_CONFIG ?= mkdocs.yml
DOCS_DIR ?= docs
SITE_DIR ?= site

# Optional:
# If you later use Material or plugins, put them here or in requirements-docs.txt
DOCS_REQUIREMENTS ?= requirements-docs.txt

PYTHON ?= python3
VENV_DIR ?= .venv
VENV_PYTHON = $(VENV_DIR)/bin/python
VENV_PIP = $(VENV_PYTHON) -m pip
MKDOCS = $(VENV_PYTHON) -m mkdocs
MKDOCS_CONFIG ?= mkdocs.yml
DOCS_DIR ?= docs
SITE_DIR ?= site
DOCS_REQUIREMENTS ?= requirements-docs.txt

DOXY          ?= doxygen
DOXY_DIR       = doxy
DOXY_OUT       = docs/doxygen/html

.PHONY: docs docs-build docs-serve docs-strict docs-deploy \
        docs-clean docs-install docs-venv docs-bootstrap venv-clean \
        docs-doxygen docs-doxygen-view docs-doxygen-clean

$(VENV_PYTHON):
	$(PYTHON) -m venv $(VENV_DIR)
	$(VENV_PIP) install --upgrade pip

docs-venv: $(VENV_PYTHON)

docs-bootstrap: docs-venv docs-install

docs: docs-build

docs-install: $(VENV_PYTHON)
ifneq ("$(wildcard $(DOCS_REQUIREMENTS))","")
	$(VENV_PIP) install -r $(DOCS_REQUIREMENTS)
else
	$(VENV_PIP) install mkdocs mkdocs-material
endif

docs-build: $(VENV_PYTHON)
	$(MKDOCS) build --config-file $(MKDOCS_CONFIG) --site-dir $(SITE_DIR) 

docs-strict: $(VENV_PYTHON)
	$(MKDOCS) build --strict --config-file $(MKDOCS_CONFIG) --site-dir $(SITE_DIR)

docs-serve: $(VENV_PYTHON)
	$(MKDOCS) serve --livereload --config-file $(MKDOCS_CONFIG)

docs-deploy: $(VENV_PYTHON)
	$(MKDOCS) gh-deploy --clean --config-file $(MKDOCS_CONFIG)

docs-clean:
	rm -rf $(SITE_DIR)

docs-doxygen:
	cd $(DOXY_DIR) && $(DOXY) Doxyfile 2>&1 | tee doxygen.log

docs-doxygen-view: docs-doxygen
	xdg-open $(DOXY_OUT)/index.html

docs-doxygen-clean:
	$(RM) -rf docs/doxygen $(DOXY_DIR)/doxygen.log

venv-clean:
	rm -rf $(VENV_DIR)

tar-gz:
	tar --exclude='.[^/]*' --exclude-vcs-ignores -czvf ../libh5cpp_${VERSION}.orig.tar.gz ./
	gpg --detach-sign --armor ../libh5cpp_${VERSION}.orig.tar.gz
	scp ../libh5cpp_${VERSION}.orig.tar.* osaka:h5cpp.org/download/

H5CPP_AMALGAMATED = h5cpp.hpp

amalgamate:
	$(PYTHON) scripts/amalgamate.py h5cpp $(H5CPP_AMALGAMATED)

clean: docs_clean
	@$(RM) h5cpp-* $(H5CPP_AMALGAMATED)

dist-debian-src: tar-gz
	debuild -i -us -uc -S

dist-debian-bin:
	debuild -i -us -uc -b

dist-debian-src-upload: dist-debian-src
	debsign -k 1B04044AF80190D78CFBE9A3B971AC62453B78AE ../libh5cpp_${VERSION}${DEB}_source.changes
	#dput mentors ../libh5cpp_${VERSION}${DEB}_source.changes

dist-rpm: dist-debian
	sudo alien -r ../libh5cpp-dev_${VERSION}_all.deb
