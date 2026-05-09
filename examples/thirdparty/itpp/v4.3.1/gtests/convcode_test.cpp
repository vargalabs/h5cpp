/*!
 * \file
 * \brief Convolutional encoder/decoder class test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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


TEST (ConvCode, All)
{
  // Test of convolutional coders

  Array<ivec> spectrum, spectrum_fast, spectrum_punct, spectrum_punct_fast;
  spectrum.set_size(2);
  spectrum_fast.set_size(2);
  spectrum_punct.set_size(2);
  spectrum_punct_fast.set_size(2);

  Convolutional_Code code;
  Punctured_Convolutional_Code code_punct;
  BPSK bpsk;
  BERC berc;

  const int no_bits = 2500;
  const int packet_size = 500;

  int coded_packet_size;
  bvec bits, tail_coded_bits, tail_decoded_bits, tailbite_coded_bits,
  tailbite_decoded_bits, trunc_coded_bits, trunc_decoded_bits;
  vec symbols;
  ivec dist_profile;

  ivec G(2);
  G(0) = 0133;
  G(1) = 0171;
  int L = max(int2bits(G(0)), int2bits(G(1))); // L = 7

  code.set_generator_polynomials(G, L);
  code_punct.set_generator_polynomials(G, L);

  bmat punct_matrix = "1 0 1; 1 1 0"; // results in R = 3/4
  code_punct.set_puncture_matrix(punct_matrix);


  //1) Rate 1/2 code
  ASSERT_FALSE(code.catastrophic());
  ASSERT_DOUBLE_EQ(0.5, code.get_rate());

  code.calculate_spectrum(spectrum, 10, 10);
  code.fast(spectrum_fast, 10, 10);
  code.distance_profile(dist_profile, 10);

  // Spectrum and Spectrum, fast
  ivec expected = "0 0 0 0 0 0 0 0 0 0 11 0 38 0 193 0 1331 0 7275 0";
  ASSERT_TRUE(expected == spectrum(0));//Ad
  ASSERT_TRUE(expected == spectrum_fast(0));
  expected = "0 0 0 0 0 0 0 0 0 0 36 0 211 0 1404 0 11633 0 77433 0";
  ASSERT_TRUE(expected == spectrum(1));//Cd
  ASSERT_TRUE(expected == spectrum_fast(1));

  // Distance profile
  expected = "2 3 3 4 4 4 4";
  ASSERT_TRUE(expected == dist_profile);

  // Tail method test. Printing 30 bits starting from bit 1400
  RNG_reset(0);
  bits = randb(no_bits);
  //Input bits
  bvec expect_bits = "0 1 1 0 1 0 1 1 1 0 1 0 1 0 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 0";
  bvec expect_enc_bits = "1 0 0 1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 0 1 0 0 0";
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  tail_coded_bits = code.encode_tail(bits);
  //Coded bits
  ASSERT_TRUE(expect_enc_bits == tail_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(tail_coded_bits, symbols);
  tail_decoded_bits = code.decode_tail(symbols);
  //Decoded bits
  ASSERT_TRUE(expect_bits == tail_decoded_bits.mid(1400, 30));
  berc.count(bits, tail_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //Tailbite method test. Printing 30 bits starting from bit 1400
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  tailbite_coded_bits = code.encode_tailbite(bits);
  ASSERT_TRUE(expect_enc_bits == tailbite_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(tailbite_coded_bits, symbols);
  tailbite_decoded_bits = code.decode_tailbite(symbols);
  //Decoded bits
  ASSERT_TRUE(expect_bits == tailbite_decoded_bits.mid(1400, 30));
  berc.clear();
  berc.count(bits, tailbite_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //Trunc method test. Printing 30 bits starting from bit 1400
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  trunc_coded_bits.set_size(0);
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_coded_bits = concat(trunc_coded_bits,
                              code.encode_trunc(bits.mid(i * packet_size,
                                                         packet_size)));
  }
  ASSERT_TRUE(expect_enc_bits == trunc_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(trunc_coded_bits, symbols);
  trunc_decoded_bits.set_size(0);
  coded_packet_size = round_i(packet_size / code.get_rate());
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_decoded_bits =
      concat(trunc_decoded_bits,
             code.decode_trunc(symbols.mid(i * coded_packet_size,
                                           coded_packet_size)));
  }
  ASSERT_TRUE(expect_bits == trunc_decoded_bits.mid(1400, 30));
  berc.clear();
  berc.count(bits, trunc_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //2) Punctured code (R = 3/4)
  ASSERT_FALSE(code_punct.catastrophic());
  ASSERT_DOUBLE_EQ(0.75, code_punct.get_rate());
  bmat expect_punct_mat = "1 0 1; 1 1 0";
  ASSERT_TRUE(expect_punct_mat == code_punct.get_puncture_matrix());

  //Tail method test. Printing 30 bits starting from bit 1400
  bits = randb(no_bits);
  expect_bits = "1 1 1 0 1 1 1 1 1 0 0 1 1 0 1 0 1 1 0 1 0 1 1 1 0 0 0 0 1 1";
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  tail_coded_bits = code_punct.encode_tail(bits);
  expect_enc_bits = "0 1 1 0 1 0 0 1 0 0 0 1 1 0 1 1 0 0 1 0 1 1 0 1 0 0 1 0 1 0";
  ASSERT_TRUE(expect_enc_bits == tail_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(tail_coded_bits, symbols);
  tail_decoded_bits = code_punct.decode_tail(symbols);
  ASSERT_TRUE(expect_bits == tail_decoded_bits.mid(1400, 30));
  berc.count(bits, tail_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //Tailbite method test. Printing 30 bits starting from bit 1400
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  tailbite_coded_bits = code_punct.encode_tailbite(bits);
  ASSERT_TRUE(expect_enc_bits == tailbite_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(tailbite_coded_bits, symbols);
  tailbite_decoded_bits = code_punct.decode_tailbite(symbols);
  ASSERT_TRUE(expect_bits == tailbite_decoded_bits.mid(1400, 30));
  berc.clear();
  berc.count(bits, tailbite_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //Trunc method test. Printing 30 bits starting from bit 1400
  ASSERT_TRUE(expect_bits == bits.mid(1400, 30));
  trunc_coded_bits.set_size(0);
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_coded_bits = concat(trunc_coded_bits,
                              code_punct.encode_trunc(bits.mid(i * packet_size,
                                                               packet_size)));
  }
  expect_enc_bits = "1 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 0 0 1 1 1 1 0 0 1 0 0 1 0 1";
  ASSERT_TRUE(expect_enc_bits == trunc_coded_bits.mid(1400, 30));
  bpsk.modulate_bits(trunc_coded_bits, symbols);
  trunc_decoded_bits.set_size(0);
  coded_packet_size = round_i(packet_size / code_punct.get_rate());
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_decoded_bits =
      concat(trunc_decoded_bits,
             code_punct.decode_trunc(symbols.mid(i * coded_packet_size,
                                                 coded_packet_size)));
  }
  ASSERT_TRUE(expect_bits == trunc_decoded_bits.mid(1400, 30));
  berc.clear();
  berc.count(bits, trunc_decoded_bits);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());

  //3) Rate 1/7 code

  ivec generator(7);
  generator(0)=02;
  generator(1)=011;
  generator(2)=015;
  generator(3)=014;
  generator(4)=07;
  generator(5)=012;
  generator(6)=06;

  code.set_generator_polynomials(generator, 4);

  bvec uncoded = "1 1 1 1 1 1 1";
  bvec coded = code.encode_tail(uncoded);
  bvec decoded = code.decode_tail(to_vec((-2)*to_ivec(coded)+1));
  ASSERT_TRUE(uncoded == decoded);
  berc.clear();
  berc.count(uncoded, decoded);
  ASSERT_DOUBLE_EQ(0, berc.get_errorrate());
}
