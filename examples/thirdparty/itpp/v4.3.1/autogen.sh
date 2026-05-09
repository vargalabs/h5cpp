#!/bin/sh
# Run this to generate all the initial makefiles, etc.

DIE=no

check_tool() {
    ($1 --version) </dev/null >/dev/null 2>&1 || {
	echo "*** Error: You must have \"$1\" installed to compile IT++ SVN sources"
        DIE=yes
    }
}

check_tool "autoconf"
check_tool "automake"
check_tool "libtoolize"
check_tool "sed"

test "$DIE" = yes && exit 1

srcdir=$(dirname $0)
test -z "${srcdir}" && srcdir=.
ORIGDIR=$(pwd)
cd "${srcdir}" || exit $?


PV=$(cat VERSION | cut -d' ' -f1)
LV=$(cat VERSION | cut -d' ' -f2)
if test "x$(cat VERSION | cut -d' ' -f3)" = "xsvn"; then
    if test -d ".git/svn"; then
        REV=$(LC_ALL=C git svn find-rev HEAD)
    elif test -d ".svn"; then
        REV=$(LC_ALL=C svn info $0 | sed -n 's/^Revision: //p')
    fi
    if test "x${REV}" != x; then
        PV="${PV}.r${REV}"
    else
        PV="${PV}.d`LC_ALL=C date +%Y%m%d`"
    fi
fi
PD=$(LC_ALL=C date +"%B %Y")

echo "Bootstapping IT++ version ${PV}"

sed -e "s/@PACKAGE_VERSION@/${PV}/" -e "s/@LIBRARY_VERSION@/${LV}/" \
    < configure.ac.in > configure.ac || exit $?
sed -e "s/@PACKAGE_VERSION@/${PV}/" \
    < itpp.spec.in > itpp.spec || exit $?
sed -e "s/@PACKAGE_VERSION@/${PV}/" -e "s/@PACKAGE_DATE@/${PD}/" \
    < itpp-config.1.in > itpp-config.1 || exit $?

test ! -d build-aux && (mkdir build-aux || exit $?)

aclocal -I m4 || exit $?
libtoolize --copy --force --automake || exit $?
aclocal -I m4 || exit $?
autoconf || exit $?
autoheader || exit $?
automake --add-missing --copy || exit $?

cd "${ORIGDIR}" || exit $?
