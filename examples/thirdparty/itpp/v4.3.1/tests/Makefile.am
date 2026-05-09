include $(top_srcdir)/Makefile.common
LDADD = $(top_builddir)/itpp/libitpp.la
AM_CXXFLAGS = $(CXXFLAGS_OPT)

BASE_TESTS = array_test bessel_test blas_test circular_buffer_test \
  converters_test fastmath_test gf2mat_test integration_test itfile_test \
  mat_test parser_test rand_test sort_test sparse_test specmat_test \
  timer_test vec_test linspace_test operators_test
BASE_LAP_TESTS = cholesky_test det_test eigen_test inv_test ls_solve_test \
  lu_test matfunc_test qr_test schur_test svd_test

COMM_TESTS = bch_test commfunc_test convcode_test error_count_test \
  galois_test interleaver_test ldpc_test llr_test modulator_test \
  pulse_shape_test rec_syst_conv_code_test reedsolomon_test turbo_test siso_test \
  exit_test stc_test demapper_test
COMM_LAP_TESTS = modulator_nd_test
COMM_FFT_TESTS = channel_test

FIXED_TESTS = fix_test

OPTIM_TESTS = newton_search_test

SIGNAL_TESTS = filter_test source_test
SIGNAL_LAP_TESTS = fastica_test poly_test
SIGNAL_FFT_TESTS = freq_filt_test sigfun_test transforms_test window_test
SIGNAL_LAP_FFT_TESTS = filter_design_test

STAT_TESTS = histogram_test
STAT_LAP_TESTS = stat_test

ALL_TESTS = $(BASE_TESTS) $(BASE_LAP_TESTS) $(COMM_TESTS) $(COMM_LAP_TESTS) \
  $(COMM_FFT_TESTS) $(FIXED_TESTS) $(OPTIM_TESTS) $(SIGNAL_TESTS) \
  $(SIGNAL_LAP_TESTS) $(SIGNAL_FFT_TESTS) $(SIGNAL_LAP_FFT_TESTS) \
  $(STAT_TESTS) $(STAT_LAP_TESTS)

USED_TESTS = $(BASE_TESTS) $(STAT_TESTS)
if RUN_LAPACK_TESTS
  USED_TESTS += $(BASE_LAP_TESTS) $(STAT_LAP_TESTS)
endif # RUN_LAPACK_TESTS
if ENABLE_COMM
  USED_TESTS += $(COMM_TESTS)
if RUN_LAPACK_TESTS
  USED_TESTS += $(COMM_LAP_TESTS)
endif # RUN_LAPACK_TESTS
if RUN_FFT_TESTS
  USED_TESTS += $(COMM_FFT_TESTS)
endif # RUN_FFT_TESTS
endif # ENABLE_COMM
if ENABLE_FIXED
  USED_TESTS += $(FIXED_TESTS)
endif # ENABLE_FIXED
if ENABLE_OPTIM
  USED_TESTS += $(OPTIM_TESTS)
endif # ENABLE_OPTIM
if ENABLE_SIGNAL
  USED_TESTS += $(SIGNAL_TESTS)
if RUN_LAPACK_TESTS
  USED_TESTS += $(SIGNAL_LAP_TESTS)
endif # RUN_LAPACK_TESTS
if RUN_FFT_TESTS
  USED_TESTS += $(SIGNAL_FFT_TESTS)
if RUN_LAPACK_TESTS
  USED_TESTS += $(SIGNAL_LAP_FFT_TESTS)
endif # RUN_LAPACK_TESTS
endif # RUN_FFT_TESTS
endif # ENABLE_SIGNAL

array_test_SOURCES = array_test.cpp
bch_test_SOURCES = bch_test.cpp
bessel_test_SOURCES = bessel_test.cpp
blas_test_SOURCES = blas_test.cpp
channel_test_SOURCES = channel_test.cpp
converters_test_SOURCES = converters_test.cpp
itfile_test_SOURCES = itfile_test.cpp
itfile_test_CPPFLAGS = -DITFILE_TEST_FILE=\"$(srcdir)/itfile_test_data.it\" $(AM_CPPFLAGS)
cholesky_test_SOURCES = cholesky_test.cpp
circular_buffer_test_SOURCES = circular_buffer_test.cpp
commfunc_test_SOURCES = commfunc_test.cpp
convcode_test_SOURCES = convcode_test.cpp
det_test_SOURCES = det_test.cpp
eigen_test_SOURCES = eigen_test.cpp
error_count_test_SOURCES = error_count_test.cpp
fastica_test_SOURCES = fastica_test.cpp
fastica_test_CPPFLAGS = -DFASTICA_TEST_FILE=\"$(srcdir)/fastica_test_data.txt\" $(AM_CPPFLAGS)
fastmath_test_SOURCES = fastmath_test.cpp
filter_test_SOURCES = filter_test.cpp
filter_design_test_SOURCES = filter_design_test.cpp
fix_test_SOURCES = fix_test.cpp
freq_filt_test_SOURCES = freq_filt_test.cpp
galois_test_SOURCES = galois_test.cpp
gf2mat_test_SOURCES = gf2mat_test.cpp
histogram_test_SOURCES = histogram_test.cpp
integration_test_SOURCES = integration_test.cpp
interleaver_test_SOURCES = interleaver_test.cpp
inv_test_SOURCES = inv_test.cpp
ldpc_test_SOURCES = ldpc_test.cpp
llr_test_SOURCES = llr_test.cpp
ls_solve_test_SOURCES = ls_solve_test.cpp
lu_test_SOURCES = lu_test.cpp
matfunc_test_SOURCES = matfunc_test.cpp
mat_test_SOURCES = mat_test.cpp
modulator_test_SOURCES = modulator_test.cpp
modulator_nd_test_SOURCES = modulator_nd_test.cpp
newton_search_test_SOURCES = newton_search_test.cpp
parser_test_SOURCES = parser_test.cpp
parser_test_CPPFLAGS = -DPARSER_TEST_FILE=\"$(srcdir)/parser_test_data.txt\" $(AM_CPPFLAGS)
poly_test_SOURCES = poly_test.cpp
pulse_shape_test_SOURCES = pulse_shape_test.cpp
qr_test_SOURCES = qr_test.cpp
rand_test_SOURCES = rand_test.cpp
rec_syst_conv_code_test_SOURCES = rec_syst_conv_code_test.cpp
reedsolomon_test_SOURCES = reedsolomon_test.cpp
schur_test_SOURCES = schur_test.cpp
sigfun_test_SOURCES = sigfun_test.cpp
sort_test_SOURCES = sort_test.cpp
source_test_SOURCES = source_test.cpp
sparse_test_SOURCES = sparse_test.cpp
specmat_test_SOURCES = specmat_test.cpp
stat_test_SOURCES = stat_test.cpp
svd_test_SOURCES = svd_test.cpp
timer_test_SOURCES = timer_test.cpp
transforms_test_SOURCES = transforms_test.cpp
turbo_test_SOURCES = turbo_test.cpp
vec_test_SOURCES = vec_test.cpp
linspace_test_SOURCES = linspace_test.cpp
operators_test_SOURCES = operators_test.cpp
window_test_SOURCES = window_test.cpp
siso_test_SOURCES = siso_test.cpp
exit_test_SOURCES = exit_test.cpp
stc_test_SOURCES = stc_test.cpp
demapper_test_SOURCES = demapper_test.cpp

check_PROGRAMS = $(USED_TESTS)
tmp_files = $(USED_TESTS:=.tmp)
ref_files = $(ALL_TESTS:=.ref)

EXTRA_DIST = itfile_test_data.it fastica_test_data.txt parser_test_data.txt \
	$(ref_files)

check-local: $(check_PROGRAMS) $(ref_files)
	@total=0; failed=0; \
	for i in $(USED_TESTS); \
	do \
		echo "------------------------------------------------------------------------------"; \
		./$$i$(EXEEXT) 2>/dev/null | sed -e 's/^M$$//' -e 's/\(-\?[0-9]\+\.\?[0-9]*e[+-]\)0\([0-9][0-9]\)/\1\2/g' > $$i.tmp; \
		$(DIFF) $(srcdir)/$$i.ref $$i.tmp >/dev/null; \
		if [ "x$$?" = "x0" ]; \
		then \
			echo "Test \`$$i' PASSED."; \
		else \
			echo "Test \`$$i' FAILED!!!"; \
			failed=$$(($$failed + 1)); \
		fi; \
		total=$$(($$total + 1)); \
	done; \
	echo "------------------------------------------------------------------------------"; \
	if [ $$failed -gt 0 ]; \
	then \
		echo ""; \
		echo "******************************************************************************"; \
		echo "*** ERROR: $$failed/$$total tests failed!!!"; \
		echo "******************************************************************************"; \
		echo ""; \
		exit 1; \
	fi

clean-local:
	-rm -f $(tmp_files) gf2mat_test.alist ldpc_test.codec
