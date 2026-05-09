/*!
 * \file
 * \brief Channel classes test program
 * \author Adam Piatyszek
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
#include <iomanip>
#include "gtest/gtest.h"

using namespace std;
using namespace itpp;

TEST(Channel, All)
{
  ostringstream ss(ostringstream::out);
  const string ref[] = {"[[0.572-0.281i 0.163-0.260i -0.243-0.025i -0.047+0.004i]\n\
 [0.917-0.028i 0.149+0.552i 0.029+0.382i -0.084-0.073i]\n\
 [0.611-0.193i -0.083-0.199i 0.194+0.068i 0.120+0.034i]\n\
 [0.975-0.149i 0.273+0.616i 0.218+0.335i -0.102+0.008i]\n\
 [0.537+0.179i 0.249+0.108i 0.194-0.144i 0.030+0.051i]\n\
 [0.574-0.822i 0.289+0.025i -0.208+0.181i 0.003-0.042i]\n\
 [0.666-0.207i 0.234+0.546i 0.175-0.032i -0.016-0.058i]\n\
 [0.746-0.281i -0.624-0.064i 0.077+0.155i 0.034-0.058i]\n\
 [0.776-0.085i -0.273+0.379i 0.129-0.407i 0.003+0.034i]\n\
 [0.694+0.193i -0.217-0.746i 0.198-0.074i -0.009-0.042i]]",
 "[[0.064+0.155i 0.074+0.540i 0.065-0.507i -0.416+0.264i -0.455-0.071i -0.203+0.058i]\n\
 [0.064+0.155i 0.074+0.540i 0.065-0.507i -0.416+0.264i -0.455-0.071i -0.203+0.058i]\n\
 [0.064+0.155i 0.074+0.540i 0.065-0.507i -0.416+0.264i -0.455-0.071i -0.203+0.058i]\n\
 [0.064+0.155i 0.074+0.540i 0.065-0.507i -0.416+0.264i -0.455-0.071i -0.203+0.058i]\n\
 [0.064+0.155i 0.074+0.540i 0.065-0.507i -0.416+0.264i -0.455-0.071i -0.203+0.058i]]",
 "[[0.063+0.143i -0.203-0.354i 0.309+0.000i -0.396+0.000i -0.285+0.000i 0.312+0.000i]\n\
 [0.062+0.142i -0.204-0.374i 0.309-0.025i -0.395+0.007i -0.285-0.008i 0.312+0.007i]\n\
 [0.062+0.140i -0.205-0.393i 0.308-0.050i -0.393+0.014i -0.285-0.016i 0.312+0.014i]\n\
 [0.061+0.138i -0.207-0.413i 0.306-0.075i -0.390+0.021i -0.285-0.024i 0.311+0.020i]\n\
 [0.060+0.135i -0.208-0.433i 0.303-0.100i -0.387+0.028i -0.284-0.032i 0.310+0.027i]\n\
 [0.060+0.132i -0.209-0.452i 0.299-0.124i -0.383+0.035i -0.283-0.039i 0.308+0.034i]\n\
 [0.059+0.128i -0.210-0.471i 0.294-0.148i -0.378+0.041i -0.281-0.048i 0.306+0.041i]\n\
 [0.059+0.123i -0.211-0.490i 0.288-0.172i -0.373+0.047i -0.279-0.056i 0.303+0.048i]\n\
 [0.058+0.118i -0.212-0.509i 0.281-0.195i -0.367+0.053i -0.276-0.063i 0.300+0.055i]\n\
 [0.058+0.112i -0.213-0.527i 0.273-0.217i -0.360+0.059i -0.273-0.071i 0.296+0.061i]]",
 "[[0.526-0.156i -0.671-0.072i 0.050-0.042i -0.106+0.044i]\n\
 [0.528-0.133i -0.685-0.062i 0.048-0.043i -0.103+0.041i]\n\
 [0.528-0.109i -0.700-0.053i 0.045-0.043i -0.100+0.038i]\n\
 [0.527-0.086i -0.714-0.044i 0.042-0.044i -0.098+0.036i]\n\
 [0.525-0.063i -0.728-0.034i 0.040-0.045i -0.095+0.033i]\n\
 [0.521-0.040i -0.743-0.025i 0.037-0.046i -0.092+0.030i]\n\
 [0.516-0.017i -0.757-0.015i 0.034-0.046i -0.089+0.028i]\n\
 [0.510+0.005i -0.771-0.006i 0.031-0.047i -0.086+0.025i]\n\
 [0.502+0.027i -0.786+0.003i 0.029-0.048i -0.083+0.022i]\n\
 [0.494+0.048i -0.800+0.013i 0.026-0.048i -0.081+0.019i]]",
 "[[0.874+0.095i -0.851+0.319i 0.178+0.255i 0.016+0.000i]\n\
 [0.873+0.112i -0.828+0.274i 0.170+0.250i 0.016-0.003i]\n\
 [0.870+0.128i -0.802+0.228i 0.162+0.245i 0.016-0.007i]\n\
 [0.865+0.144i -0.774+0.182i 0.153+0.239i 0.016-0.010i]\n\
 [0.859+0.159i -0.743+0.137i 0.143+0.233i 0.016-0.014i]\n\
 [0.852+0.173i -0.709+0.091i 0.132+0.225i 0.016-0.017i]\n\
 [0.843+0.187i -0.673+0.045i 0.120+0.217i 0.016-0.021i]\n\
 [0.833+0.200i -0.635+0.000i 0.108+0.209i 0.016-0.025i]\n\
 [0.822+0.213i -0.595-0.045i 0.096+0.200i 0.017-0.029i]\n\
 [0.809+0.226i -0.552-0.089i 0.082+0.190i 0.017-0.032i]]"};
  int i = 0;

  ss.setf(ios::fixed);
  ss.precision(3);

  double Ts = 100e-9; // channel sampling time Ts = 100ns
  double fd_norm = 100e3 * Ts; // normalised Doppler fd_norm = 10kHz * Ts

  RNG_reset(0);

  Channel_Specification cost207ra(COST207_RA);
  Channel_Specification cost207bu(COST207_BU);

  TDL_Channel tdl_207ra(cost207ra, Ts);
  TDL_Channel tdl_207bu(cost207bu, Ts);

  cmat ch_coeffs;

  tdl_207ra.set_fading_type(Independent);
  tdl_207ra.generate(10, ch_coeffs);
  ss << ch_coeffs;
  ASSERT_TRUE(ss.str() == ref[i++]);
  ss.str("");

  tdl_207bu.set_fading_type(Static);
  tdl_207bu.generate(5, ch_coeffs);
  ss << ch_coeffs;
  ASSERT_TRUE(ss.str() == ref[i++]);
  ss.str("");

  tdl_207bu.set_fading_type(Correlated);
  tdl_207bu.set_norm_doppler(fd_norm);
  tdl_207bu.generate(10, ch_coeffs);
  ss << round_to_zero(ch_coeffs);
  ASSERT_TRUE(ss.str() == ref[i++]);
  ss.str("");

  tdl_207ra.set_correlated_method(FIR);
  tdl_207ra.set_norm_doppler(fd_norm);
  tdl_207ra.generate(10, ch_coeffs);
  ss << round_to_zero(ch_coeffs);
  ASSERT_TRUE(ss.str() == ref[i++]);
  ss.str("");

  tdl_207ra.set_correlated_method(IFFT);
  tdl_207ra.generate(200, ch_coeffs);
  ss << round_to_zero(ch_coeffs.get_rows(0, 9));
  ASSERT_TRUE(ss.str() == ref[i++]);
}
