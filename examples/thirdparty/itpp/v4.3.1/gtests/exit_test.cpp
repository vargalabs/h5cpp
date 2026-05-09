/*!
 * \file
 * \brief EXIT class test program
 * \author Bogdan Cristea
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

#include "itpp/itcomm.h"
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST(EXIT, All)
{
  //test parameters
  const double eps = 1e-12;
  const vec apriori_mutual_info_expect = "1.803346244821924e-05 0.48594415413293501 0.91282228577448177 0.99902335921776853";
  const vec extrinsic_mutual_info_expect = "0.57474253793992003 0.79081789562325711 0.98191798133410091 0.9975400908252704";
  const vec extrinsic_mutual_info_p_expect = "0.031706232459936402 0.33837811448717559 0.93275900297128034 0.99710134578550191";

  //general parameters
  vec sigmaA = "0.01 2 4 7";//standard deviation (sqrt(variance)) of the mutual a priori information
  double threshold_value = 50;
  string map_metric = "maxlogMAP";
  ivec gen = "07 05";//octal form
  int constraint_length = 3;
  int nb_blocks_lim = 10;
  int perm_len = int(itpp::pow10(3.0));//total number of bits in a block (with tail)
  double EbN0_dB = 0.8;
  double R = 1.0 / 3.0; //coding rate of PCCC
  double Ec = 1.0;//coded bit energy

  //other parameters
  vec sigma2A = sqr(sigmaA);
  int sigma2A_len = sigma2A.length();
  int nb_bits = perm_len - (constraint_length - 1); //number of bits in a block (without tail)
  double sigma2 = (0.5 * Ec / R) * pow(inv_dB(EbN0_dB), -1.0); //N0/2
  double Lc = -2 / sigma2; //normalisation factor for intrinsic information (take into account the BPSK modulation)
  bvec bits(nb_bits);
  bvec tail;
  bvec bits_tail(perm_len);
  bmat parity_bits;
  int coded_bits_len = 2 * perm_len;
  bvec coded_bits(coded_bits_len);
  vec mod_bits(coded_bits_len);
  vec rec_sig(coded_bits_len);
  vec intrinsic_coded(coded_bits_len);
  vec intrinsic_coded_p(2 * nb_bits);
  intrinsic_coded_p.zeros();
  vec apriori_data(perm_len);
  vec extrinsic_coded;
  vec extrinsic_data;
  vec apriori_mutual_info(sigma2A_len);
  vec extrinsic_mutual_info(sigma2A_len);
  vec extrinsic_mutual_info_p(sigma2A_len);
  extrinsic_mutual_info.zeros();
  extrinsic_mutual_info_p.zeros();
  register int en, n, nb_blocks;

  //Recursive Systematic Convolutional Code
  Rec_Syst_Conv_Code rsc;
  rsc.set_generator_polynomials(gen, constraint_length);//initial state should be the zero state

  //BPSK modulator
  BPSK bpsk;

  //AWGN channel
  AWGN_Channel channel;
  channel.set_noise(sigma2);

  //SISO module
  SISO siso;
  siso.set_generators(gen, constraint_length);
  siso.set_map_metric(map_metric);

  //EXIT chart
  EXIT exit;

  //Randomize generators
  RNG_reset(12345);

  //main loop
  for(en = 0; en < sigma2A_len; en++) {
    apriori_mutual_info(en) = exit.apriori_mutual_info(sigma2A(en));//a priori mutual info
    for(nb_blocks = 0; nb_blocks < nb_blocks_lim; nb_blocks++) {
      //bits generation
      bits = randb(nb_bits);

      //RSC code
      rsc.encode_tail(bits, tail, parity_bits);//tail is added

      //form coder output
      bits_tail = concat(bits, tail);
      for(n = 0; n < perm_len; n++) {
        coded_bits(2 * n) = bits_tail(n); //systematic output
        coded_bits(2 * n + 1) = parity_bits(n, 0); //parity output
      }

      //BPSK modulation (1->-1,0->+1)
      mod_bits = bpsk.modulate_bits(coded_bits);

      //AWGN channel
      rec_sig = channel(mod_bits);

      //first SISO RSC module  (tail is added)
      //intrinsic info. of coded bits
      intrinsic_coded = Lc * rec_sig;

      //a priori info. of data bits
      apriori_data = exit.generate_apriori_info(bits_tail);

      //SISO RSC module
      siso.rsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data, true);

      //threshold
      extrinsic_data = SISO::threshold(extrinsic_data, threshold_value);

      //extrinsic mutual info
      extrinsic_mutual_info(en) += exit.extrinsic_mutual_info(extrinsic_data.left(nb_bits), bits);

      //second SISO RSC module (no tail added)
      //intrinsic info. of coded bits
      for(n = 0; n < nb_bits; n++)
        intrinsic_coded_p(2 * n + 1) = Lc * rec_sig(2 * n + 1); //parity bits only

      //a priori info. of data bits
      apriori_data = exit.generate_apriori_info(bits);

      //SISO RSC module
      siso.rsc(extrinsic_coded, extrinsic_data, intrinsic_coded_p, apriori_data, false);

      //threshold
      extrinsic_data = SISO::threshold(extrinsic_data, threshold_value);

      //extrinsic mutual info
      extrinsic_mutual_info_p(en) += exit.extrinsic_mutual_info(extrinsic_data, bits);
    }//end blocks (while loop)

    //mean extrinsic mutual info over all blocks
    extrinsic_mutual_info(en) /= nb_blocks_lim;
    extrinsic_mutual_info_p(en) /= nb_blocks_lim;

    //check results
    ASSERT_NEAR(apriori_mutual_info_expect(en), apriori_mutual_info(en), eps);
    ASSERT_NEAR(extrinsic_mutual_info_expect(en), extrinsic_mutual_info(en), eps);
    ASSERT_NEAR(extrinsic_mutual_info_p_expect(en), extrinsic_mutual_info_p(en), eps);
  }
}
