include $(top_srcdir)/Makefile.common
LDADD = $(top_builddir)/itpp/libitpp.la
AM_CXXFLAGS = $(CXXFLAGS_OPT)

ref_example_sources = bpsk.ref convcode.ref interleaver.ref	\
  qpsk_simulation.ref reedsolomon.ref spread.ref timer.ref	\
  vector_and_matrix.ref

EXTRA_DIST = $(ref_example_sources)

bpsk_example_SOURCES = bpsk.cpp
convcode_example_SOURCES = convcode.cpp
interleaver_example_SOURCES = interleaver.cpp
ldpc_bersim_awgn_example_SOURCES = ldpc_bersim_awgn.cpp
pccc_bersim_awgn_example_SOURCES = pccc_bersim_awgn.cpp
exit_pccc_example_SOURCES = exit_pccc.cpp
sccc_bersim_awgn_example_SOURCES = sccc_bersim_awgn.cpp
turbo_equalizer_bersim_multipath_example_SOURCES = turbo_equalizer_bersim_multipath.cpp
ldpc_gen_codes_example_SOURCES = ldpc_gen_codes.cpp
mimoconv_example_SOURCES = mimoconv.cpp
mog_example_SOURCES = mog.cpp
qpsk_simulation_example_SOURCES = qpsk_simulation.cpp
rayleigh_example_SOURCES = rayleigh.cpp
read_it_file_example_SOURCES = read_it_file.cpp
reedsolomon_example_SOURCES = reedsolomon.cpp
spread_example_SOURCES = spread.cpp
timer_example_SOURCES = timer.cpp
vector_and_matrix_example_SOURCES = vector_and_matrix.cpp
write_it_file_example_SOURCES = write_it_file.cpp
stbicm_example_SOURCES = stbicm.cpp
dump_it_file_example_SOURCES = dump_it_file.cpp

BASE_EXAMPLES = mog_example read_it_file_example timer_example	\
  write_it_file_example dump_it_file_example
if RUN_LAPACK_TESTS
  BASE_EXAMPLES += vector_and_matrix_example
endif # RUN_LAPACK_TESTS

COMM_EXAMPLES = bpsk_example convcode_example interleaver_example		\
  ldpc_gen_codes_example ldpc_bersim_awgn_example qpsk_simulation_example	\
  rayleigh_example reedsolomon_example spread_example pccc_bersim_awgn_example \
  exit_pccc_example sccc_bersim_awgn_example turbo_equalizer_bersim_multipath_example
if RUN_LAPACK_TESTS
  COMM_EXAMPLES += mimoconv_example stbicm_example
endif # RUN_LAPACK_TESTS

USED_EXAMPLES = $(BASE_EXAMPLES)
if ENABLE_COMM
  USED_EXAMPLES += $(COMM_EXAMPLES)
endif # ENABLE_COMM

check_PROGRAMS = $(USED_EXAMPLES)
