/*!
 * \file
 * \brief Reed-Solomon encoder/decoder class test program
 * \author Steve Peters and Adam Piatyszek
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
using std::cout;
using std::endl;


int main()
{
  cout << "==========================================" << endl;
  cout << "   Test of Reed-Solomon encoder/decoder   " << endl;
  cout << "==========================================" << endl;

  bmat u = randb(8, 15);
  bmat c(u.rows(), 21);
  bmat y(u.rows(), 21);
  bvec codeword, errorword;
  bmat decoded(u.rows(), u.cols());
  Reed_Solomon rs(3, 1);
  Reed_Solomon rs_sys(3, 1, true);

  bmat f = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0";

  cout << "Non-systematic case" << endl;
  cout << "-------------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    cout << "Info word:       " << u.get_row(i) << endl;
    codeword = rs.encode(u.get_row(i));
    it_assert(c.cols() == length(codeword), "Error 1");
    c.set_row(i, rs.encode(u.get_row(i)));
    cout << "Encoded:         " << c.get_row(i) << endl;
    errorword = f.get_row(i) + c.get_row(i);
    it_assert(y.cols() == length(errorword), "Error 2");
    y.set_row(i, f.get_row(i) + c.get_row(i));
    cout << "One error added: " << y.get_row(i) << endl;
    decoded.set_row(i, rs.decode(y.get_row(i)));
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  cout << "Systematic case" << endl;
  cout << "---------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    cout << "Info word:       " << u.get_row(i) << endl;
    cout << "Encoded:         " << c.get_row(i) << endl;
    y.set_row(i, f.get_row(i) + c.get_row(i));
    cout << "One error added: " << y.get_row(i) << endl;
    decoded.set_row(i, rs_sys.decode(y.get_row(i)));
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  bmat g = "1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1; 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0";

  cout << "Systematic case with decoder failure indicator" << endl;
  cout << "----------------------------------------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    cout << "Info word:       " << u.get_row(i) << endl;
    cout << "Encoded:         " << c.get_row(i) << endl;
    y.set_row(i, g.get_row(i) + c.get_row(i));
    cout << "Two error added: " << y.get_row(i) << endl;
    bvec message, cw_is_valid;
    rs_sys.decode(y.get_row(i), message, cw_is_valid);
    decoded.set_row(i, message);
    cout << "Decoded to:      " << decoded.get_row(i) << endl;
    cout << "Codeword valid:  " << cw_is_valid << endl << endl;
  }


  rs = Reed_Solomon(3, 1, false, 0);
  rs_sys = Reed_Solomon(3, 1, true, 0);

  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "Erasure decoding (1 erasure) and non-narrow-sense" << endl;
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  cout << "Non-systematic case" << endl;
  cout << "-------------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    cout << "Info word:       " << u.get_row(i) << endl;
    codeword = rs.encode(u.get_row(i));
    it_assert(c.cols() == length(codeword), "Error 1");
    c.set_row(i, rs.encode(u.get_row(i)));
    cout << "Encoded:         " << c.get_row(i) << endl;
    it_assert(y.cols() == length(errorword), "Error 2");
    y.set_row(i, c.get_row(i));
    bvec decoded_tmp;
    bvec cw_isvalid;
    ivec erasure(1);
    erasure(0) = i/y.cols();
    rs.decode(y.get_row(i), erasure, decoded_tmp, cw_isvalid);
    decoded.set_row(i, decoded_tmp);
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  cout << "Systematic case" << endl;
  cout << "---------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    cout << "Info word:       " << u.get_row(i) << endl;
    cout << "Encoded:         " << c.get_row(i) << endl;
    y.set_row(i, c.get_row(i));
    bvec decoded_tmp;
    bvec cw_isvalid;
    ivec erasure(1);
    erasure(0) = 0;//i/y.cols();
    rs_sys.decode(y.get_row(i), erasure, decoded_tmp, cw_isvalid);
    decoded.set_row(i, decoded_tmp);
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  return 0;
}
