/*!
 * \file
 * \brief Random number generator test program
 * \author Adam Piatyszek
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

using namespace itpp;
using namespace std;

int main()
{
  cout.setf(ios::fixed);
  cout.precision(2);

  RNG_reset(4357U);
  Bernoulli_RNG b_rng;
  cout << "Bernoulli_RNG:\n" << b_rng() << endl;
  cout << b_rng(10) << endl;
  cout << b_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Uniform_RNG u_rng;
  cout << "Uniform_RNG:\n" << u_rng() << endl;
  cout << u_rng(10) << endl;
  cout << u_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  I_Uniform_RNG iu_rng(0, 9);
  cout << "I_Uniform_RNG [0..9]:\n" << iu_rng() << endl;
  cout << iu_rng(10) << endl;
  cout << iu_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Normal_RNG n_rng;
  cout << "Normal_RNG:\n" << n_rng() << endl;
  cout << n_rng(10) << endl;
  cout << n_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Complex_Normal_RNG cn_rng;
  cout << "Complex_Normal_RNG:\n" << cn_rng() << endl;
  cout << cn_rng(10) << endl;
  cout << cn_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Exponential_RNG e_rng;
  cout << "Exponential_RNG:\n" << e_rng() << endl;
  cout << e_rng(10) << endl;
  cout << e_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Gamma_RNG g_rng;
  cout << "Gamma_RNG:\n" << g_rng() << endl;
  cout << g_rng(10) << endl;
  cout << g_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Laplace_RNG lap_rng;
  cout << "Laplace_RNG:\n" << lap_rng() << endl;
  cout << lap_rng(10) << endl;
  cout << lap_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  AR1_Normal_RNG ar1_rng;
  cout << "AR1_Normal_RNG:\n" << ar1_rng() << endl;
  cout << ar1_rng(10) << endl;
  cout << ar1_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Weibull_RNG w_rng;
  cout << "Weibull_RNG:\n" << w_rng() << endl;
  cout << w_rng(10) << endl;
  cout << w_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Rayleigh_RNG ray_rng;
  cout << "Rayleigh_RNG:\n" << ray_rng() << endl;
  cout << ray_rng(10) << endl;
  cout << ray_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Rice_RNG ric_rng;
  cout << "Rice_RNG:\n" << ric_rng() << endl;
  cout << ric_rng(10) << endl;
  cout << ric_rng(3, 5) << endl << endl;

  return 0;
}
