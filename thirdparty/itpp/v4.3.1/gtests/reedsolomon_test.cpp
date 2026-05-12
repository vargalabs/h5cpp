/*!
 * \file
 * \brief Reed-Solomon encoder/decoder class test program
 * \author Steve Peters and Adam Piatyszek
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

TEST(ReedSolomon, All)
{
  // Test of Reed-Solomon encoder/decoder
  RNG_reset(0);

  bmat u = randb(8, 15);
  bmat c(u.rows(), 21);
  bmat y(u.rows(), 21);
  bvec codeword, errorword;
  bmat decoded(u.rows(), u.cols());
  Reed_Solomon rs(3, 1);
  Reed_Solomon rs_sys(3, 1, true);

  bmat f = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0";

  // Non-systematic case
  for (int i = 0; i < u.rows(); i++) {
    codeword = rs.encode(u.get_row(i));
    ASSERT_TRUE(c.cols() == length(codeword)) << "Error 1";
    c.set_row(i, rs.encode(u.get_row(i)));
    errorword = f.get_row(i) + c.get_row(i);
    ASSERT_TRUE(y.cols() == length(errorword)) << "Error 2";
    y.set_row(i, f.get_row(i) + c.get_row(i));
    decoded.set_row(i, rs.decode(y.get_row(i)));
    ASSERT_TRUE(u.get_row(i) == decoded.get_row(i)) << i;
  }

  // Systematic case
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    y.set_row(i, f.get_row(i) + c.get_row(i));
    decoded.set_row(i, rs_sys.decode(y.get_row(i)));
    ASSERT_TRUE(u.get_row(i) == decoded.get_row(i)) << i;
  }

  bmat g = "1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1; 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0";

  // Systematic case with decoder failure indicator
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    y.set_row(i, g.get_row(i) + c.get_row(i));
    bvec message, cw_is_valid;
    rs_sys.decode(y.get_row(i), message, cw_is_valid);
    decoded.set_row(i, message);
    //decoded bits are wrong without a decoding failure
    switch (i) {
      case 0:
      case 6:
      ASSERT_TRUE(u.get_row(i) == decoded.get_row(i));
      break;
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 7:
      ASSERT_FALSE(u.get_row(i) == decoded.get_row(i));
      break;
      default:
      void(0);
    }
  }

  rs = Reed_Solomon(3, 1, false, 0);
  rs_sys = Reed_Solomon(3, 1, true, 0);

  // Erasure decoding (1 erasure) and non-narrow-sense
  // Non-systematic case
  for (int i = 0; i < u.rows(); i++) {
    codeword = rs.encode(u.get_row(i));
    ASSERT_TRUE(c.cols() == length(codeword)) << "Error 1";
    c.set_row(i, rs.encode(u.get_row(i)));
    ASSERT_TRUE(y.cols() == length(errorword)) << "Error 2";
    y.set_row(i, c.get_row(i));
    bvec decoded_tmp;
    bvec cw_isvalid;
    ivec erasure(1);
    erasure(0) = i/y.cols();
    rs.decode(y.get_row(i), erasure, decoded_tmp, cw_isvalid);
    decoded.set_row(i, decoded_tmp);
    ASSERT_TRUE(u.get_row(i) == decoded.get_row(i));
  }

  // Systematic case
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    y.set_row(i, c.get_row(i));
    bvec decoded_tmp;
    bvec cw_isvalid;
    ivec erasure(1);
    erasure(0) = 0;//i/y.cols();
    rs_sys.decode(y.get_row(i), erasure, decoded_tmp, cw_isvalid);
    decoded.set_row(i, decoded_tmp);
    ASSERT_TRUE(u.get_row(i) == decoded.get_row(i));
  }
}
