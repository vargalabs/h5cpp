/*!
 * \file
 * \brief Bessel test program
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

#include <itpp/itbase.h>
#include <iomanip>

using namespace std;
using namespace itpp;

int main(void)
{

  cout << "================================" << endl;
  cout << "    Test of bessel functions " << endl;
  cout << "================================" << endl;

  vec x = linspace(0.01, 10, 20);

  cout << "x = " << x << endl;

  cout << "besselj(0, x) = " << fixed << besselj(0, x) << endl;
  cout << "besselj(1, x) = " << besselj(1, x) << endl;
  cout << "besselj(5, x) = " << round_to_infty(besselj(5, x)) << endl;
  cout << "besselj(0.3, x) = " << besselj(0.3, x) << endl;
  cout << "besselj(1.7, x) = " << besselj(1.7, x) << endl;
  cout << "besselj(5.3, x) = " << round_to_infty(besselj(5.3, x)) << endl;

  cout << "bessely(0, x) = " << bessely(0, x) << endl;
  cout << "bessely(1, x) = " << bessely(1, x) << endl;
  cout << "bessely(5, x) = " << round_to_infty(bessely(5, x)) << endl;
  cout << "bessely(0.3, x) = " << bessely(0.3, x) << endl;
  cout << "bessely(1.7, x) = " << bessely(1.7, x) << endl;
  cout << "bessely(5.3, x) = " << round_to_infty(bessely(5.3, x)) << endl;

  cout << "besseli(0, x) = " << besseli(0, x) << endl;
  cout << "besseli(1, x) = " << besseli(1, x) << endl;
  cout << "besseli(5, x) = " << besseli(5, x) << endl;
  cout << "besseli(0.3, x) = " << besseli(0.3, x) << endl;
  cout << "besseli(1.7, x) = " << besseli(1.7, x) << endl;
  cout << "besseli(5.3, x) = " << besseli(5.3, x) << endl;

  cout << "besselk(0, x) = " << besselk(0, x) << endl;
  cout << "besselk(1, x) = " << besselk(1, x) << endl;
  cout << "besselk(5, x) = " << round_to_infty(besselk(5, x)) << endl;

  return 0;

}
