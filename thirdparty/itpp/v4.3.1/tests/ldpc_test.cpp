/*!
 * \file
 * \brief LDPC class test program
 * \author Erik G. Larsson
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

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  LDPC_Parity_Regular H;
  H.generate(200, 3, 6, "rand", "100 6");
  int girth = H.cycle_removal_MGW(6);
  cout << "girth=" << girth << endl;
  LDPC_Generator_Systematic G;
  G.construct(&H);
  LDPC_Code C(&H, &G);
  C.save_code("ldpc_test.codec");
  LDPC_Code C1("ldpc_test.codec", &G);
  cout << C << endl;
  bvec bitsin = randb(C.get_nvar() - C.get_ncheck());
  bvec bitsout;
  C.encode(bitsin, bitsout);
  it_assert(C.syndrome_check(bitsout), "syndrome check failed");

  double EbN0db = 1.5;
  double N0 = pow(10.0, -EbN0db / 10.0) / C.get_rate();
  double sigma = sqrt(N0 / 2.0);
  vec x = 1.0 + sigma * randn(C.get_nvar());
  QLLRvec LLRin = C.get_llrcalc().to_qllr(2.0 * x / (N0 / 2.0));
  QLLRvec LLRout(C.get_nvar());
  C.bp_decode(LLRin, LLRout);
  cout << LLRout.left(25) << endl;

  // BLDPC code
  {
    cout.precision(5);
    imat A = "0 -1 -1 0; -1 1 4 -1; -1 2 -1 6";
    imat B = "1; -1; 2";
    imat T = "0 -1 -1; 0 0 -1; -1 0 0";
    imat C = "3 -1 5 -1";
    imat D = "3";
    imat E = "-1 -1 0";
    // base matrix
    imat H_b = concat_vertical(concat_horizontal(concat_horizontal(A, B), T),
                               concat_horizontal(concat_horizontal(C, D), E));
    int Z = 4; // expansion factor
    cout << endl;
    cout << "expansion factor Z = " << Z << endl;
    cout << "base matrix H_b =\n" << H_b << endl;

    BLDPC_Parity H(H_b, Z);
    cout << "expanded parity check matrix H =\n" << H.get_H() << endl;

    BLDPC_Generator G(&H);
    bvec in_bits = randb(H.get_nvar() - H.get_ncheck());
    bvec codeword;
    G.encode(in_bits, codeword);
    cout << in_bits << endl << codeword << endl;
  }
}
