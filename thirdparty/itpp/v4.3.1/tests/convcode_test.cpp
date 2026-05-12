/*!
 * \file
 * \brief Convolutional encoder/decoder class test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

using namespace itpp;
using namespace std;


int main()
{
  cout << "====================================" << endl;
  cout << "    Test of convolutional coders    " << endl;
  cout << "====================================" << endl;

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


  cout << "------------------------------------------------------------------------------" << endl;
  cout << "1) Rate 1/2 code" << endl;
  cout << "------------------------------------------------------------------------------" << endl;

  cout << "Catastrophic test = " << code.catastrophic() << endl;
  cout << "Code rate         = " << code.get_rate() << endl << endl;

  code.calculate_spectrum(spectrum, 10, 10);
  code.fast(spectrum_fast, 10, 10);
  code.distance_profile(dist_profile, 10);

  cout << "Spectrum:" << endl;
  cout << "* Ad = " << spectrum(0) << endl;
  cout << "* Cd = " << spectrum(1) << endl;

  cout << "Spectrum, fast:" << endl;
  cout << "* Ad = " << spectrum_fast(0) << endl;
  cout << "* Cd = " << spectrum_fast(1) << endl << endl;

  cout << "Distance profile  = " << dist_profile << endl << endl;

  cout << "Tail method test. Printing 30 bits starting from bit 1400:" << endl;
  bits = randb(no_bits);
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tail_coded_bits = code.encode_tail(bits);
  cout << "* Coded bits    = " << tail_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tail_coded_bits, symbols);
  tail_decoded_bits = code.decode_tail(symbols);
  cout << "* Decoded bits  = " << tail_decoded_bits.mid(1400, 30) << endl;
  berc.count(bits, tail_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Tailbite method test. Printing 30 bits starting from bit 1400:"
       << endl;
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tailbite_coded_bits = code.encode_tailbite(bits);
  cout << "* Coded bits    = " << tailbite_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tailbite_coded_bits, symbols);
  tailbite_decoded_bits = code.decode_tailbite(symbols);
  cout << "* Decoded bits  = " << tailbite_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, tailbite_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Trunc method test. Printing 30 bits starting from bit 1400:" << endl;
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  trunc_coded_bits.set_size(0);
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_coded_bits = concat(trunc_coded_bits,
                              code.encode_trunc(bits.mid(i * packet_size,
                                                         packet_size)));
  }
  cout << "* Coded bits    = " << trunc_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(trunc_coded_bits, symbols);
  trunc_decoded_bits.set_size(0);
  coded_packet_size = round_i(packet_size / code.get_rate());
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_decoded_bits =
      concat(trunc_decoded_bits,
             code.decode_trunc(symbols.mid(i * coded_packet_size,
                                           coded_packet_size)));
  }
  cout << "* Decoded bits  = " << trunc_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, trunc_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;


  cout << "------------------------------------------------------------------------------" << endl;
  cout << "2) Punctured code (R = 3/4)" << endl;
  cout << "------------------------------------------------------------------------------" << endl;

  cout << "Catastrophic test = " << code_punct.catastrophic() << endl;
  cout << "Code rate         = " << code_punct.get_rate() << endl;
  cout << "Puncture matrix   = " << code_punct.get_puncture_matrix() << endl
       << endl;

  cout << "Tail method test. Printing 30 bits starting from bit 1400:" << endl;
  bits = randb(no_bits);
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tail_coded_bits = code_punct.encode_tail(bits);
  cout << "* Coded bits    = " << tail_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tail_coded_bits, symbols);
  tail_decoded_bits = code_punct.decode_tail(symbols);
  cout << "* Decoded bits  = " << tail_decoded_bits.mid(1400, 30) << endl;
  berc.count(bits, tail_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Tailbite method test. Printing 30 bits starting from bit 1400:"
       << endl;
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  tailbite_coded_bits = code_punct.encode_tailbite(bits);
  cout << "* Coded bits    = " << tailbite_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(tailbite_coded_bits, symbols);
  tailbite_decoded_bits = code_punct.decode_tailbite(symbols);
  cout << "* Decoded bits  = " << tailbite_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, tailbite_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "Trunc method test. Printing 30 bits starting from bit 1400:" << endl;
  cout << "* Input bits    = " << bits.mid(1400, 30) << endl;
  trunc_coded_bits.set_size(0);
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_coded_bits = concat(trunc_coded_bits,
                              code_punct.encode_trunc(bits.mid(i * packet_size,
                                                               packet_size)));
  }
  cout << "* Coded bits    = " << trunc_coded_bits.mid(1400, 30) << endl;
  bpsk.modulate_bits(trunc_coded_bits, symbols);
  trunc_decoded_bits.set_size(0);
  coded_packet_size = round_i(packet_size / code_punct.get_rate());
  for (int i = 0; i < no_bits / packet_size; i++) {
    trunc_decoded_bits =
      concat(trunc_decoded_bits,
             code_punct.decode_trunc(symbols.mid(i * coded_packet_size,
                                                 coded_packet_size)));
  }
  cout << "* Decoded bits  = " << trunc_decoded_bits.mid(1400, 30) << endl;
  berc.clear();
  berc.count(bits, trunc_decoded_bits);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  cout << "------------------------------------------------------------------------------" << endl;
  cout << "3) Rate 1/7 code" << endl;
  cout << "------------------------------------------------------------------------------" << endl;

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
  cout << "* Input bits    = " << uncoded << endl;
  cout << "* Decoded bits  = " << decoded << endl;
  berc.clear();
  berc.count(uncoded, decoded);
  cout << "BER = " << berc.get_errorrate() << endl << endl;

  return 0;
}
