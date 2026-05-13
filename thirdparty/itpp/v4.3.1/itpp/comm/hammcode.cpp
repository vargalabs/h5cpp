/*!
 * \file
 * \brief Implementation of a Hamming code class
 * \author Tony Ottosson
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

#include <itpp/comm/hammcode.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/converters.h>


namespace itpp
{

Hamming_Code::Hamming_Code(int m)
{
  n = pow2i(m) - 1;
  k = pow2i(m) - m - 1;
  H.set_size(n - k, n);
  G.set_size(k, n);
  generate_H(); // generate_H must be run before generate_G
  generate_G();
}

void Hamming_Code::generate_H(void)
{
  int i, j, NextPos;
  char NotUsed;
  bvec temp;
  ivec indexes(n);
  indexes.zeros();

  for (i = 1; i <= n - k; i++) { indexes(i - 1) = pow2i(n - k - i); }
  NextPos = n - k;
  for (i = 1; i <= n; i++) {
    NotUsed = 1;
    for (j = 0; j < n; j++)
      if (i == indexes(j)) { NotUsed = 0; }
    if (NotUsed) { indexes(NextPos) = i; NextPos = NextPos + 1; }
  }

  for (i = 0; i < n; i++) {
    temp = dec2bin(n - k, indexes(i)); //<-CHECK THIS OUT!!!!
    for (j = 0; j < (n - k); j++) {
      H(j, i) = temp(j);
    }
  }
}

void Hamming_Code::generate_G(void)
{
  int i, j;
  for (i = 0; i < k; i++) {
    for (j = 0; j < n - k; j++)
      G(i, j) = H(j, i + n - k);
  }

  for (i = 0; i < k; i++) {
    for (j = n - k; j < n; j++)
      G(i, j) = 0;
  }

  for (i = 0; i < k; i++)
    G(i, i + n - k) = 1;
}

void Hamming_Code::encode(const bvec &uncoded_bits, bvec &coded_bits)
{
  int length = uncoded_bits.length();
  int Itterations = floor_i(static_cast<double>(length) / k);
  bmat Gt = G.T();
  int i;

  coded_bits.set_size(Itterations * n, false);
  //Code all codewords
  for (i = 0; i < Itterations; i++)
    coded_bits.replace_mid(n*i, Gt * uncoded_bits.mid(i*k, k));
}

bvec Hamming_Code::encode(const bvec &uncoded_bits)
{
  bvec coded_bits;
  encode(uncoded_bits, coded_bits);
  return coded_bits;
}

void Hamming_Code::decode(const bvec &coded_bits, bvec &decoded_bits)
{
  int length = coded_bits.length();
  int Itterations = floor_i(static_cast<double>(length) / n);
  ivec Hindexes(n);
  bvec temp(n - k);
  bvec coded(n), syndrome(n - k);
  int isynd, errorpos = 0;
  int i, j;

  decoded_bits.set_size(Itterations*k, false);

  for (i = 0; i < n; i++) {
    for (j = 0; j < n - k; j++)
      temp(j) = H(j, i);
    Hindexes(i) = bin2dec(temp);
  }

  //Decode all codewords
  for (i = 0; i < Itterations; i++) {
    coded = coded_bits.mid(i * n, n);
    syndrome = H * coded;
    isynd = bin2dec(syndrome);
    if (isynd != 0) {
      for (j = 0; j < n; j++)
        if (Hindexes(j) == isynd) { errorpos = j; };
      coded(errorpos) += 1;
    }
    decoded_bits.replace_mid(k*i, coded.right(k));
  }
}

bvec Hamming_Code::decode(const bvec &coded_bits)
{
  bvec decoded_bits;
  decode(coded_bits, decoded_bits);
  return decoded_bits;
}


// -------------- Soft-decision decoding is not implemented ----------------
void Hamming_Code::decode(const vec &, bvec &)
{
  it_error("Hamming_Code::decode(vec, bvec); soft-decision decoding is not implemented");
}

bvec Hamming_Code::decode(const vec &)
{
  it_error("Hamming_Code::decode(vec, bvec); soft-decision decoding is not implemented");
  return bvec();
}


} // namespace itpp
