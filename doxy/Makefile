# H5CPP Doxygen — run from this directory
# Outputs to ../docs/doxygen/html/

DOXY    ?= doxygen
OUT     = ../docs/doxygen/html

.PHONY: all view clean

all:
	$(DOXY) Doxyfile 2>&1 | tee doxygen.log

view: all
	xdg-open $(OUT)/index.html

clean:
	$(RM) -rf ../docs/doxygen doxygen.log
