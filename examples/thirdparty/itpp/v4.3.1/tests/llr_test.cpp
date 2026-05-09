/*!
 * \file
 * \brief Test program for the LLR class
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
#include <iomanip>

using namespace itpp;
using namespace std;


int main()
{
  // This is a hack for improper rounding under MinGW
  cout.precision(13);

  LLR_calc_unit lcu1; // standard table resolution
  LLR_calc_unit lcu2(10, 7, 9);  // low table resolution
  LLR_calc_unit lcu3(2, 15, 0);  // low table resolution and low LLR granuality
  LLR_calc_unit lcu4(10, 0, 0);  // this gives logexp=logmax
  cout << lcu1 << endl;
  cout << lcu2 << endl;
  cout << lcu3 << endl;
  cout << lcu4 << endl;

  cout << "Testing Jacobian logarithm with four different resolutions." << endl;
  for (double x = 0.0; x < 10; x += 0.1) {
    cout << "JacLog(" << x << ") = "
         << lcu1.to_double(lcu1.logexp(lcu1.to_qllr(x))) << " ; "
         << lcu2.to_double(lcu2.logexp(lcu2.to_qllr(x))) << " ; "
         << lcu3.to_double(lcu3.logexp(lcu3.to_qllr(x))) << " ; "
         << lcu4.to_double(lcu4.logexp(lcu4.to_qllr(x)))
         << endl;
  }

  cout << "-------------------" << endl;
  cout << "Some special cases:" << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0), lcu1.to_qllr(0.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0), lcu1.to_qllr(-0.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0), lcu1.to_qllr(0.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0), lcu1.to_qllr(-0.75)))
       << endl;

  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0), lcu1.to_qllr(0.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0), lcu1.to_qllr(-0.75)))
       << endl;

  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75), lcu1.to_qllr(-1.25)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-1.25), lcu1.to_qllr(3.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-3.75), lcu1.to_qllr(1.25)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25), lcu1.to_qllr(-3.75)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75), lcu1.to_qllr(1.25)))
       << endl;
  cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25), lcu1.to_qllr(3.75)))
       << endl;

  return 0;
}

