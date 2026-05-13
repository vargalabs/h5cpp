/*!
 * \file
 * \brief Implementation of a CRC code class
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

#include <itpp/comm/crc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>


namespace itpp
{

void CRC_Code::set_generator(const bvec &poly)
{
  //it_assert(poly(0) == 1 && poly(poly.size()-1) == 1, "CRC_Code::set_polynomial: not a valid polynomial");
  it_assert(poly(0) == 1, "CRC_Code::set_polynomial: not a valid polynomial");
  polynomial = poly;
  no_parity = polynomial.size() - 1;
}

//! \cond

std::string crccode[18][2] = {
  {"CRC-4", "1 1 1 1 1"},
  {"CRC-7", "1 1 0 1 0 0 0 1"},
  {"CRC-8", "1 1 1 0 1 0 1 0 1"},
  {"CRC-12", "1 1 0 0 0 0 0 0 0 1 1 1 1"},
  {"CRC-24", "1 1 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 1"},
  {"CRC-32", "1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1 1 0 0 0 1 1 1 0 0 0 1 0"},
  {"CCITT-4", "1 0 0 1 1"},
  {"CCITT-5", "1 1 0 1 0 1"},
  {"CCITT-6", "1 0 0 0 0 1 1"},
  {"CCITT-16", "1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1"},
  {"CCITT-32", "1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1"},
  {"WCDMA-8", "1 1 0 0 1 1 0 1 1"},
  {"WCDMA-12", "1 1 0 0 0 0 0 0 0 1 1 1 1"},
  {"WCDMA-16", "1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1"},
  {"WCDMA-24", "1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1"},
  {"ATM-8", "1 0 0 0 0 0 1 1 1"},
  {"ANSI-16", "1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1"},
  {"SDLC-16", "1 1 0 1 0 0 0 0 0 1 0 0 1 0 1 1 1"},
};

//! \endcond

void CRC_Code::set_code(const std::string &code)
{
  bvec poly;
  for (int i = 0; i < 18;i++) {
    if (crccode[i][0] == code)
      poly = bvec(crccode[i][1]);
  }

  if ((code == "WCDMA-8") || (code == "WCDMA-12") || (code == "WCDMA-16") || (code == "WCDMA-24")) {
    reverse_parity = true;
  }

  it_assert(poly.size() > 0, "This CRC code doesn't exist in the tables");
  set_generator(poly);
}

// Not optimized for speed!
void CRC_Code::parity(const bvec &in_bits, bvec &out) const
{
  bvec temp = concat(in_bits, zeros_b(no_parity));

  for (int i = 0; i < temp.size() - polynomial.size() + 1; i++) {
    if (temp(i) == 1) {
      temp.set_subvector(i, temp(i, i + no_parity) + polynomial);
    }
  }

  out = temp(temp.size() - no_parity, temp.size() - 1);

  if (reverse_parity) {
    out = reverse(out);
  }

}

// Not optimized for speed
bool CRC_Code::check_parity(const bvec &coded_bits) const
{
  int n = coded_bits.size();
  bvec temp;

  if (reverse_parity) {
    temp = concat(coded_bits.left(n - no_parity), reverse(coded_bits.right(no_parity)));
  }
  else {
    temp = coded_bits;
  }

  for (int i = 0; i < temp.size() - polynomial.size() + 1; i++) {
    if (temp(i) == 1) {
      temp.set_subvector(i, temp(i, i + no_parity) + polynomial);
    }
  }

  if (temp(temp.size() - no_parity, temp.size() - 1) == zeros_b(no_parity))
    return true;
  else
    return false;
}

void CRC_Code::encode(const bvec &in_bits, bvec &out) const
{
  bvec p;
  parity(in_bits, p);
  out = concat(in_bits, p);
}

bvec CRC_Code::encode(const bvec &in_bits) const
{
  bvec temp;
  encode(in_bits, temp);
  return temp;
}

bool CRC_Code::decode(const bvec &coded_bits, bvec &out) const
{
  out = coded_bits(0, coded_bits.size() - no_parity - 1);
  if (check_parity(coded_bits)) {
    return true;
  }
  else
    return false;
}

bool CRC_Code::decode(bvec &coded_bits) const
{
  //coded_bits = coded_bits(0, coded_bits.size()-no_parity-1); <-- OLD CODE
  if (check_parity(coded_bits)) {
    return true;
  }
  else
    return false;
}

} // namespace itpp
