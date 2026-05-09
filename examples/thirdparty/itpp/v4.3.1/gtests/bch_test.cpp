/*!
 * \file
 * \brief BCH encoder/decoder class test program
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
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


bvec set_errors(const bvec &input, const ivec errpos)
{
  bvec output = input;
  for (int i = 0; i < errpos.length(); i++)
    output(errpos(i)) ^= 1;
  return output ;
}


TEST(BCH, codec)
{
  ostringstream ss(ostringstream::out);
  const string ref[] = {"encoded = [1 0 1 1 1 0 0 1 1 1 1 1 0 0 1 0 0 1 0 1 1 0 1 0 0 0 0 1 1 0 0]",
    "err =     [1 1 0 1 1 0 0 1 1 1 1 1 0 0 1 0 0 1 0 1 1 0 1 0 0 0 0 1 1 0 0]",
    "input =   [0 1 0 0 0 1 1 1 0 1 0 1 1 1 0 1 0 0 1 1 0]",
    "encoded = [0 1 1 1 0 0 1 1 1 1 1 1 0 1 0 0 1 0 1 1 0 0 1 1 1 0 1 0 1 1 0]",
    "err =     [0 0 0 1 0 0 1 1 1 1 1 1 0 1 0 0 1 0 1 1 0 0 1 1 1 0 1 1 1 1 0]",
    "decoded = [0 0 0 1 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 1 0]",
    "input =   [1 0 1 1 0 0 0 1 1 0 1 0 0 1 1 1 1 1 1 0 0]",
    "encoded = [1 0 1 1 0 0 0 1 1 0 1 0 0 1 1 1 0 1 0 1 0 1 0 0 1 0 0 0 1 0 1]",
    "err =     [1 1 0 1 0 0 0 1 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 0 1 0 0 1 1 0 1]",
    "decoded = [1 1 0 1 0 0 0 1 1 0 1 0 0 1 0 1]"};
  int i = 0;

  {
    BCH bch(31, 2);
    RNG_reset(0);

    bvec input = randb(21);
    bvec encoded = bch.encode(input);
    bvec err = set_errors(encoded, (ivec) "1 2"); // error positions
    bvec decoded = bch.decode(err);

    // A two error case (should be corrected)
    ASSERT_TRUE(input == decoded);
    ss << "encoded = " << encoded;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "err =     " << err;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");

    input = randb(21);
    encoded = bch.encode(input);
    err = set_errors(encoded, (ivec) "1 2 27"); // error positions
    decoded = bch.decode(err);

    // A three error case (will cause decoding errors);
    ss << "input =   " << input;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "encoded = " << encoded;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "err =     " << err;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "decoded = " << decoded;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
  }

  // Systematic vs. non-systematic test

  {
    bmat u = "0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; 0 1 1 0";
    bmat c(u.rows(), 7);
    bmat y(u.rows(), 7);
    bmat decoded(u.rows(), u.cols());
    BCH bch_nsys(7, 1);
    BCH bch_sys(7, 1, true);

    bmat f = "1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1";

    // Non-systematic case
    for (int i = 0; i < u.rows(); i++) {
      c.set_row(i, bch_nsys.encode(u.get_row(i)));
      y.set_row(i, f.get_row(i) + c.get_row(i));
      decoded.set_row(i, bch_nsys.decode(y.get_row(i)));
    }
    ASSERT_TRUE(u == decoded);

    // Systematic case
    for (int i = 0; i < u.rows(); i++) {
      c.set_row(i, bch_sys.encode(u.get_row(i)));
      y.set_row(i, f.get_row(i) + c.get_row(i));
      decoded.set_row(i, bch_sys.decode(y.get_row(i)));
    }
    ASSERT_TRUE(u == decoded);
  }

  // Systematic decoding failure test

  {
    BCH bch(31, 3, true);

    bvec input = randb(21);
    bvec encoded = bch.encode(input);
    bvec err = set_errors(encoded, (ivec) "1 2 14 27"); // error positions

    bvec decoded;
    bvec is_valid_cw;    // test the new decoding procedure for the systematic case (should extract the systematics)
    // all codewords valid?
    ASSERT_FALSE(bch.decode(err, decoded, is_valid_cw));
    // valid codeword?
    ASSERT_TRUE(bvec("0") == is_valid_cw);

    // A four error case (will cause decoding failure)
    ss << "input =   " << input;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "encoded = " << encoded;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "err =     " << err;
    ASSERT_TRUE(ss.str() == ref[i++]);
    ss.str("");
    ss << "decoded = " << decoded;
    ASSERT_TRUE(ss.str() == ref[i++]);
  }
}
