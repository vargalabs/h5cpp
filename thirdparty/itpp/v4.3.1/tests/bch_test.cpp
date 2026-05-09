/*!
 * \file
 * \brief BCH encoder/decoder class test program
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
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


bvec set_errors(const bvec &input, const ivec errpos)
{
  bvec output = input;
  for (int i = 0; i < errpos.length(); i++)
    output(errpos(i)) ^= 1;
  return output ;
}


int main()
{
  cout << "===================================" << endl;
  cout << "    Test of BCH encoder/decoder    " << endl;
  cout << "===================================" << endl;

  {
    BCH bch(31, 2);

    bvec input = randb(21);
    bvec encoded = bch.encode(input);
    bvec err = set_errors(encoded, (ivec) "1 2"); // error positions
    bvec decoded = bch.decode(err);

    cout << "A two error case (should be corrected)" << endl;
    cout << "input =   " << input << endl ;
    cout << "encoded = " << encoded << endl ;
    cout << "err =     " << err <<  endl ;
    cout << "decoded = " << decoded << endl ;

    input = randb(21);
    encoded = bch.encode(input);
    err = set_errors(encoded, (ivec) "1 2 27"); // error positions
    decoded = bch.decode(err);

    cout << "A three error case (will cause decoding errors)" << endl;
    cout << "input =   " << input << endl;
    cout << "encoded = " << encoded << endl;
    cout << "err =     " << err <<  endl;
    cout << "decoded = " << decoded << endl << endl;

  }

  cout << "========================================" << endl;
  cout << "   Systematic vs. non-systematic test   " << endl;
  cout << "========================================" << endl;

  {
    bmat u = "0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; 0 1 1 0";
    bmat c(u.rows(), 7);
    bmat y(u.rows(), 7);
    bmat decoded(u.rows(), u.cols());
    BCH bch_nsys(7, 1);
    BCH bch_sys(7, 1, true);

    bmat f = "1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1";

    cout << "Non-systematic case" << endl;
    cout << "-------------------" << endl;
    for (int i = 0; i < u.rows(); i++) {
      c.set_row(i, bch_nsys.encode(u.get_row(i)));
      cout << "Encoded " << u.get_row(i) << " to " << c.get_row(i) << endl;
      y.set_row(i, f.get_row(i) + c.get_row(i));
      cout << "One error added: " << y.get_row(i) << endl;
      decoded.set_row(i, bch_nsys.decode(y.get_row(i)));
      cout << "Decoded to:" << decoded.get_row(i) << endl << endl;
    }

    cout << "Systematic case" << endl;
    cout << "---------------" << endl;
    for (int i = 0; i < u.rows(); i++) {
      c.set_row(i, bch_sys.encode(u.get_row(i)));
      cout << "Encoded " << u.get_row(i) << " to " << c.get_row(i) << endl;
      y.set_row(i, f.get_row(i) + c.get_row(i));
      cout << "One error added: " << y.get_row(i) << endl;
      decoded.set_row(i, bch_sys.decode(y.get_row(i)));
      cout << "Decoded to:" << decoded.get_row(i) << endl << endl;
    }
  }

  cout << "========================================" << endl;
  cout << "   Systematic decoding failure test     " << endl;
  cout << "========================================" << endl;

  {
    BCH bch(31, 3, true);

    bvec input = randb(21);
    bvec encoded = bch.encode(input);
    bvec err = set_errors(encoded, (ivec) "1 2 14 27"); // error positions

    bvec decoded;
    bvec is_valid_cw;    // test the new decoding procedure for the systematic case (should extract the systematics)
    cout << "all codewords valid? " << bch.decode(err, decoded, is_valid_cw) << endl;
    cout << "valid codeword? = " << is_valid_cw << endl;

    cout << "A four error case (will cause decoding failure)" << endl;
    cout << "input =   " << input << endl;
    cout << "encoded = " << encoded << endl;
    cout << "err =     " << err <<  endl;
    cout << "decoded = " << decoded << endl << endl << endl;
  }

  return 0;
}
