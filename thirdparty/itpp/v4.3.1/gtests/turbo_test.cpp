/*!
 * \file
 * \brief Turbo encoder/decoder class test program
 * \author Pal Frenger and Erik G. Larsson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
void assert_mat_p(const mat &ref, const mat &act, int line, double tol = 1e-4)
{
  ASSERT_EQ(ref.rows(), act.rows()) << line;
  ASSERT_EQ(ref.cols(), act.cols()) << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol) << line;
    }
  }
}
#define assert_mat(ref, act) assert_mat_p(ref, act, __LINE__)
#define assert_mat_tol(ref, act, tol) assert_mat_p(ref, act, __LINE__, tol)

TEST(Turbo, All)
{
  RNG_reset(12345);

  Turbo_Codec turbo;
  ivec gen(2);
  gen(0) = 07;
  gen(1) = 05;
  int constraint_length = 3;
  int block_length = 400;
  ivec interleaver_sequence = wcdma_turbo_interleaver_sequence(block_length);
  int iterations = 8;
  string metric = "LOGMAX";
  double logmax_scale_factor = 0.7;
  bool adaptive_stop = true;
  turbo.set_parameters(gen, gen, constraint_length, interleaver_sequence,
                       iterations, metric, logmax_scale_factor,
                       adaptive_stop);
  int num_blocks = 50;

  bmat puncture_matrix = "1 1;1 1;1 1";
  Punctured_Turbo_Codec pturbo_ref;

  pturbo_ref.set_parameters(gen, gen, constraint_length, interleaver_sequence,
                            puncture_matrix, iterations, metric, logmax_scale_factor,
                            adaptive_stop);


  vec EbN0db = "0.0 0.5 1.0 1.5 2.0";
  double A = 1.0;
  double Ts = 1.0;
  double Ec = A * A * Ts;
  double r = block_length / (3.0 * block_length + 8.0); // 8 tailbits
  double Eb = Ec / r;
  vec EbN0 = pow(10.0, 0.1 * EbN0db);
  vec N0 = Eb * pow(EbN0, -1.0);
  vec sigma2 = N0 / 2;
  ivec nrof_used_iterations;

  vec symbols, received;
  bvec input, coded_bits, decoded_bits, transmitted;
  bvec decoded_bits_p;

  Normal_RNG noise_src;

  BPSK bpsk;
  BERC berc;

  mat err = zeros(4, EbN0db.length());
  mat cor = zeros(4, EbN0db.length());
  mat ber = zeros(4, EbN0db.length());
  mat avg_nrof_iterations = zeros(4, EbN0db.length());
  LLR_calc_unit lowresllrcalc(10, 7, 9);  // table with low resolution
  Array<Real_Timer> timer(4);
  for(int i = 0; i < 4; i++) {
    timer(i).reset();
  }

  for(int i = 0; i < EbN0db.length(); i++) {
    noise_src.setup(0.0, sigma2(i));
    turbo.set_awgn_channel_parameters(Ec, N0(i));
    pturbo_ref.set_awgn_channel_parameters(Ec, N0(i));
    input = randb(block_length * num_blocks);

    turbo.encode(input, transmitted);
    bpsk.modulate_bits(transmitted, symbols);
    received = symbols + noise_src(transmitted.length());

    // -- logmax decoding --
    turbo.set_metric("LOGMAX", 1.0);
    timer(0).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(0).stop();
    pturbo_ref.set_metric("LOGMAX", 1.0);
    pturbo_ref.decode(received, decoded_bits_p, nrof_used_iterations);
    ASSERT_TRUE(decoded_bits == decoded_bits_p);

    berc.clear();
    berc.count(input, decoded_bits);
    err(0, i) = berc.get_errors();
    cor(0, i) = berc.get_corrects();
    ber(0, i) = berc.get_errorrate();
    avg_nrof_iterations(0, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- logmap decoding --
    turbo.set_metric("LOGMAP", 1.0);
    timer(1).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(1).stop();
    pturbo_ref.set_metric("LOGMAP", 1.0);
    pturbo_ref.decode(received, decoded_bits_p, nrof_used_iterations);
    ASSERT_TRUE(decoded_bits == decoded_bits_p);

    berc.clear();
    berc.count(input, decoded_bits);
    err(1, i) = berc.get_errors();
    cor(1, i) = berc.get_corrects();
    ber(1, i) = berc.get_errorrate();
    avg_nrof_iterations(1, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, default resolution --
    turbo.set_metric("TABLE", 1.0);
    timer(2).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(2).stop();
    pturbo_ref.set_metric("TABLE", 1.0);
    pturbo_ref.decode(received, decoded_bits_p, nrof_used_iterations);
    ASSERT_TRUE(decoded_bits == decoded_bits_p);

    berc.clear();
    berc.count(input, decoded_bits);
    err(2, i) = berc.get_errors();
    cor(2, i) = berc.get_corrects();
    ber(2, i) = berc.get_errorrate();
    avg_nrof_iterations(2, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, low resolution --
    turbo.set_metric("TABLE", 1.0, lowresllrcalc);
    timer(3).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(3).stop();
    pturbo_ref.set_metric("TABLE", 1.0, lowresllrcalc);
    pturbo_ref.decode(received, decoded_bits_p, nrof_used_iterations);
    ASSERT_TRUE(decoded_bits == decoded_bits_p);

    berc.clear();
    berc.count(input, decoded_bits);
    err(3, i) = berc.get_errors();
    cor(3, i) = berc.get_corrects();
    ber(3, i) = berc.get_errorrate();
    avg_nrof_iterations(3, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

  }

  mat ref = "0.1174 0.0586 0.005 0.00065 0;"
            "0.0518 0.02005 0.0017 0 0;"
            "0.0522 0.0199 0.0017 0 0;"
            "0.05495 0.01905 0.0017 0 0";
  assert_mat(ref, ber);
  ref = "7.7 7.04 4.7 3.46 3.02;"
        "7.16 6.06 4.12 3.24 3;"
        "7.2 6.06 4.08 3.24 3;"
        "7.26 6.1 4.06 3.22 3.06";
  assert_mat(ref, avg_nrof_iterations);
  ref = "2348 1172 100 13 0;"
               "1036 401 34 0 0;"
               "1044 398 34 0 0;"
               "1099 381 34 0 0";
  assert_mat(ref, err);
  ref = "17652 18828 19900 19987 20000;"
          "18964 19599 19966 20000 20000;"
          "18956 19602 19966 20000 20000;"
          "18901 19619 19966 20000 20000";
  assert_mat(ref, cor);

  /* ######################################
   * # now test for punctured turbo codec #
   * ######################################
   */

  puncture_matrix = "1 1;1 0;0 1";
  Punctured_Turbo_Codec pturbo;

  pturbo.set_parameters(gen, gen, constraint_length, interleaver_sequence,
                        puncture_matrix, iterations, metric, logmax_scale_factor,
                        adaptive_stop);
  r = pturbo.get_rate();
  Eb = Ec / r;
  N0 = Eb * pow(EbN0, -1.0);
  sigma2 = N0 / 2;

  err.zeros();
  cor.zeros();
  ber.zeros();
  avg_nrof_iterations.zeros();

  for(int i = 0; i < 4; i++) {
    timer(i).reset();
  }

  for(int i = 0; i < EbN0db.length(); i++) {
    noise_src.setup(0.0, sigma2(i));
    pturbo.set_awgn_channel_parameters(Ec, N0(i));
    input = randb(block_length * num_blocks);

    pturbo.encode(input, transmitted);
    bpsk.modulate_bits(transmitted, symbols);
    received = symbols + noise_src(transmitted.length());

    // -- logmax decoding --
    pturbo.set_metric("LOGMAX", 1.0);
    timer(0).start();
    pturbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(0).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(0, i) = berc.get_errors();
    cor(0, i) = berc.get_corrects();
    ber(0, i) = berc.get_errorrate();
    avg_nrof_iterations(0, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- logmap decoding --
    pturbo.set_metric("LOGMAP", 1.0);
    timer(1).start();
    pturbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(1).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(1, i) = berc.get_errors();
    cor(1, i) = berc.get_corrects();
    ber(1, i) = berc.get_errorrate();
    avg_nrof_iterations(1, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, default resolution --
    pturbo.set_metric("TABLE", 1.0);
    timer(2).start();
    pturbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(2).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(2, i) = berc.get_errors();
    cor(2, i) = berc.get_corrects();
    ber(2, i) = berc.get_errorrate();
    avg_nrof_iterations(2, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, low resolution --
    pturbo.set_metric("TABLE", 1.0, lowresllrcalc);
    timer(3).start();
    pturbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(3).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(3, i) = berc.get_errors();
    cor(3, i) = berc.get_corrects();
    ber(3, i) = berc.get_errorrate();
    avg_nrof_iterations(3, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

  }

  ref = "0.1489 0.0907 0.0519 0.00355 0.00195;"
        "0.1068 0.0428 0.0242 0.0002 0;"
        "0.10685 0.04285 0.02435 0.0002 0;"
        "0.10745 0.04645 0.0244 0.00035 0";
  assert_mat(ref, ber);
  ref = "8 7.74 6.92 4.42 3.26;"
        "6.38 7.26 6.18 4.08 3.32;"
        "6.48 7.2 6.2 4.08 3.3;"
        "7.88 7.44 6.42 4.02 3.28";
  assert_mat(ref, avg_nrof_iterations);
  ref = "2978 1814 1038 71 39;"
          "2136 856 484 4 0;"
          "2137 857 487 4 0;"
          "2149 929 488 7 0";
  assert_mat_tol(ref, err, 2);
  ref = "17022 18186 18962 19929 19961;"
          "17864 19144 19516 19996 20000;"
          "17863 19143 19513 19996 20000;"
          "17851 19071 19512 19993 20000";
  assert_mat_tol(ref, cor, 2);

  /*
  // The test program cannot print this out, but on my system
  // the QLLR based decoder is about 8 times faster than logmap. -EGL
  cerr << "Timers: ";
  for (int i=0; i<4; i++) { cerr << timer(i).get_time() << "  "; }
  cerr << endl;
  */

  //test interleavers
  ivec interleaver;
  ivec block_lengths = "40 48 56 64 72 80 88 96 104 112 120 128 136 144 152 "
                       "160 168 176 184 192 200 208 216 224 232 240 248 256 264 272 280 288 296 "
                       "304 312 320 328 336 344 352 360 368 376 384 392 400 "
                       "408 416 424 432 440 448 456 464 472 480 488 496 504 512 528 544 560 576 "
                       "592 608 624 640 656 672 688 704 720 736 752 768 784 "
                       "800 816 832 848 864 880 896 912 928 944 960 976 992 1008 1024 1056 1088 "
                       "1120 1152 1184 1216 1248 1280 1312 1344 1376 1408 "
                       "1440 1472 1504 1536 1568 1600 1632 1664 1696 1728 1760 1792 1824 1856 1888 "
                       "1920 1952 1984 2016 2048 2112 2176 2240 2304 2368 "
                       "2432 2496 2560 2624 2688 2752 2816 2880 2944 3008 3072 3136 3200 3264 3328 "
                       "3392 3456 3520 3584 3648 3712 3776 3840 3904 3968 "
                       "4032 4096 4160 4224 4288 4352 4416 4480 4544 4608 4672 4736 4800 4864 4928 "
                       "4992 5056 5120 5184 5248 5312 5376 5440 5504 5568 "
                       "5632 5696 5760 5824 5888 5952 6016 6080 6144";

  for(int i = 0; i < block_lengths.length(); ++i) {
    if(5114 >= block_lengths[i]) {  //use allowed lengths
      interleaver = wcdma_turbo_interleaver_sequence(block_lengths[i]);
      sort(interleaver);
      for(int j = 0; j < block_lengths[i]; ++j) {
        ASSERT_TRUE(j == interleaver[j]) << "WCDMA: wrong value for intl length " << block_lengths[i];
      }
    }

    interleaver = lte_turbo_interleaver_sequence(block_lengths[i]);
    sort(interleaver);
    for(int j = 0; j < block_lengths[i]; ++j) {
      ASSERT_TRUE(j == interleaver[j]) << "LTE: wrong value for intl length " << block_lengths[i];
    }
  }
}
