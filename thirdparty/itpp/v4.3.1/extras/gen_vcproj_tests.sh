#!/bin/sh
TESTS_DIR="../tests"
MKL_VCPROJ_DIR="../win32/itpp_mkl_tests"
ACML_VCPROJ_DIR="../win32/itpp_acml_tests"
for i in "${TESTS_DIR}"/*.cpp; do
	TEST_FILE=`basename "${i}" .cpp`
	sed -e "s/TEST_BASENAME/${TEST_FILE}/g" acml_tests_vcproj.template > \
		"${ACML_VCPROJ_DIR}/${TEST_FILE}.vcproj"
	sed -e "s/TEST_BASENAME/${TEST_FILE}/g" mkl_tests_vcproj.template > \
		"${MKL_VCPROJ_DIR}/${TEST_FILE}.vcproj"
done
