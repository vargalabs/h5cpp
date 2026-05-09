/*!
 * \file
 * \brief Implementation of the Extended Golay Code (24, 12, 8)
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

#include <itpp/comm/egolay.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/converters.h>

namespace itpp
{

Extended_Golay::Extended_Golay(void)
{
  B = "0 1 1 1 1 1 1 1 1 1 1 1;1 1 1 0 1 1 1 0 0 0 1 0;1 1 0 1 1 1 0 0 0 1 0 1;1 0 1 1 1 0 0 0 1 0 1 1;1 1 1 1 0 0 0 1 0 1 1 0;1 1 1 0 0 0 1 0 1 1 0 1;1 1 0 0 0 1 0 1 1 0 1 1;1 0 0 0 1 0 1 1 0 1 1 1;1 0 0 1 0 1 1 0 1 1 1 0;1 0 1 0 1 1 0 1 1 1 0 0;1 1 0 1 1 0 1 1 1 0 0 0;1 0 1 1 0 1 1 1 0 0 0 1";

  G = concat_horizontal(eye_b(12), B);
}

void Extended_Golay::encode(const bvec &uncoded_bits, bvec &coded_bits)
{
  int no_bits = uncoded_bits.length();
  int no_blocks = floor_i(no_bits / 12.0);

  coded_bits.set_size(24*no_blocks, false);
  bmat Gt = G.T();
  int i;

  for (i = 0; i < no_blocks; i++)
    coded_bits.replace_mid(24*i, Gt * uncoded_bits.mid(i*12, 12));
}

bvec Extended_Golay::encode(const bvec &uncoded_bits)
{
  bvec coded_bits;
  encode(uncoded_bits, coded_bits);
  return coded_bits;
}

void Extended_Golay::decode(const bvec &coded_bits, bvec &decoded_bits)
{
  int no_bits = coded_bits.length();
  int no_blocks = floor_i(no_bits / 24.0);

  decoded_bits.set_size(12*no_blocks, false);
  int i, j;
  bvec S(12), BS(12), r(12), temp(12), e(24), c(24);
  bmat eyetemp = eye_b(12);

  for (i = 0; i < no_blocks; i++) {
    r = coded_bits.mid(i * 24, 24);
    // Step 1. Compute S=G*r.
    S = G * r;
    // Step 2. w(S)<=3. e=(S,0). Goto 8.
    if (weight(S) <= 3) {
      e = concat(S, zeros_b(12));
      goto Step8;
    }

    // Step 3. w(S+Ii)<=2. e=(S+Ii,yi). Goto 8.
    for (j = 0; j < 12; j++) {

      temp = S + B.get_col(j);
      if (weight(temp) <= 2) {
        e = concat(temp, eyetemp.get_row(j));
        goto Step8;
      }
    }

    // STEP 4. Compute B*S
    BS = B * S;

    // Step 5. w(B*S)<=3. e=(0,BS). Goto8.
    if (weight(BS) <= 3) {
      e = concat(zeros_b(12), BS);
      goto Step8;
    }

    // Step 6. w(BS+Ri)<=2. e=(xi,BS+Ri). Goto 8.
    for (j = 0; j < 12; j++) {
      temp = BS + B.get_row(j);
      if (weight(temp) <= 2) {
        e = concat(eyetemp.get_row(j), temp);
        goto Step8;
      }
    }

    // Step 7. Uncorrectable erreor pattern. Choose the first 12 bits.
    e = zeros_b(24);
    goto Step8;

  Step8: // Step 8. c=r+e. STOP
    c = r + e;
    decoded_bits.replace_mid(i*12, c.left(12));
  }
}

bvec Extended_Golay::decode(const bvec &coded_bits)
{
  bvec decoded_bits;
  decode(coded_bits, decoded_bits);
  return decoded_bits;
}


// -------------- Soft-decision decoding is not implemented ------------------
void Extended_Golay::decode(const vec &, bvec &)
{
  it_error("Extended_Golay::decode(vec, bvec); soft-decision decoding is not implemented");
}

bvec Extended_Golay::decode(const vec &)
{
  it_error("Extended_Golay::decode(vec, bvec); soft-decision decoding is not implemented");
  return bvec();
}


} // namespace itpp
