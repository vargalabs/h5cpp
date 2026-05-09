/*!
 * \file
 * \brief Channel classes test program
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

#include <itpp/itcomm.h>
#include <iomanip>

using namespace std;
using namespace itpp;

int main()
{
  cout.setf(ios::fixed);
  cout.precision(3);

  double Ts = 100e-9; // channel sampling time Ts = 100ns
  double fd_norm = 100e3 * Ts; // normalised Doppler fd_norm = 10kHz * Ts

  Channel_Specification cost207ra(COST207_RA);
  Channel_Specification cost207bu(COST207_BU);

  TDL_Channel tdl_207ra(cost207ra, Ts);
  TDL_Channel tdl_207bu(cost207bu, Ts);

  cmat ch_coeffs;

  tdl_207ra.set_fading_type(Independent);
  tdl_207ra.generate(10, ch_coeffs);
  cout << "Independent Fading Generator\n"
       << "----------------------------\n"
       << ch_coeffs << "\n\n";

  tdl_207bu.set_fading_type(Static);
  tdl_207bu.generate(5, ch_coeffs);
  cout << "Static Fading Generator\n"
       << "-----------------------\n"
       << ch_coeffs << "\n\n";

  tdl_207bu.set_fading_type(Correlated);
  tdl_207bu.set_norm_doppler(fd_norm);
  tdl_207bu.generate(10, ch_coeffs);
  cout << "Correlated Fading Generator (Rice method)\n"
       << "-----------------------------------------\n"
       << round_to_zero(ch_coeffs) << "\n\n";

  tdl_207ra.set_correlated_method(FIR);
  tdl_207ra.set_norm_doppler(fd_norm);
  tdl_207ra.generate(10, ch_coeffs);
  cout << "Correlated Fading Generator (FIR method)\n"
       << "----------------------------------------\n"
       << round_to_zero(ch_coeffs) << "\n\n";

  tdl_207ra.set_correlated_method(IFFT);
  tdl_207ra.generate(200, ch_coeffs);
  cout << "Correlated Fading Generator (IFFT method)\n"
       << "-----------------------------------------\n"
       << round_to_zero(ch_coeffs.get_rows(0, 9)) << "\n\n";

  return 0;
}
