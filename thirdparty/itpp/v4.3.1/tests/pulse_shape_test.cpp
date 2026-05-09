/*!
 * \file
 * \brief Pulse shaping classes test program
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

using namespace std;
using namespace itpp;

int main()
{
  cout << "======================================" << endl;
  cout << "    Test of pulse shaping routines    " << endl;
  cout << "======================================" << endl << endl;

  BPSK bpsk;
  vec symbols, samples, rsymbols;
  Root_Raised_Cosine<double> rrc_tx(0.5), rrc_rx(0.5);
  Raised_Cosine<double> rc_tx(0.5);

  bpsk.modulate_bits(randb(20), symbols);
  samples = rrc_tx.shape_symbols(symbols);
  rsymbols = rrc_rx.shape_samples(samples);

  cout << "*** Root Raised Cosine; real input ***" << endl << endl;
  cout << "pulse, RRC = " << round_to_zero(rrc_tx.get_pulse_shape())
       << endl << endl;
  cout << "symbols = " << round_to_zero(symbols) << endl << endl;
  cout << "samples = " << round_to_zero(samples) << endl << endl;
  cout << "received symbols =" << fixed << round_to_zero(rsymbols) << endl << endl;

  samples = rc_tx.shape_symbols(symbols);

  cout << "*** Raised Cosine; real input ***" << endl << endl;
  cout << "pulse, RC = " << fixed << round_to_zero(rc_tx.get_pulse_shape())
       << endl << endl;
  cout << "symbols = " << fixed << round_to_zero(symbols) << endl << endl;
  cout << "samples = " << fixed << round_to_zero(samples) << endl << endl;

  QPSK qpsk;
  cvec csymbols, csamples, crsymbols;
  Root_Raised_Cosine<complex<double> > crrc_tx(0.5), crrc_rx(0.5);
  Raised_Cosine<complex<double> > crc_tx(0.5);

  qpsk.modulate_bits(randb(40), csymbols);
  csamples = crrc_tx.shape_symbols(csymbols);
  crsymbols = crrc_rx.shape_samples(csamples);

  cout << "*** Root Raised Cosine; complex input ***" << endl << endl;
  cout << "pulse, RRC = " << fixed << round_to_zero(crrc_tx.get_pulse_shape())
       << endl << endl;
  cout << "symbols = " << fixed << round_to_zero(csymbols) << endl << endl;
  cout << "samples = " << fixed << round_to_zero(csamples) << endl << endl;
  cout << "received symbols = " << fixed << round_to_zero(crsymbols) << endl << endl;

  csamples = crc_tx.shape_symbols(csymbols);

  cout << "*** Raised Cosine; complex input ***" << endl << endl;
  cout << "pulse, RC = " << fixed << round_to_zero(crc_tx.get_pulse_shape())
       << endl << endl;
  cout << "symbols = " << fixed << round_to_zero(csymbols) << endl << endl;
  cout << "samples = " << fixed << round_to_zero(csamples) << endl << endl;

  return 0;
}
