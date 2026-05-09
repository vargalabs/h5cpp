/*!
 * \file
 * \brief Interface of an Orthogonal Frequency Division Multiplex
 * (OFDM) class
 * \author Pal Frenger, Anders Persson and Tony Ottosson
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

#ifndef OFDM_H
#define OFDM_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \ingroup modulators
  \brief Class for modulating and demodulation of OFDM signals using the FFT

  The modulated signal is normalized taking into account the cyclic prefix
*/
class ITPP_EXPORT OFDM
{
public:
  //! Empty constructor
  OFDM(void) { setup_done = false; }
  //! Constructor \a Nfft is the size of the FFT. \a Ncp is the length of the cyclic prefix. \a Nupsample is the upsampling factor (default=1)
  OFDM(int inNfft, int inNcp, int inNupsample = 1);
  //! Return the number of carriers
  int no_carriers() {return Nfft;}
  //! Set parameters
  void set_parameters(const int Nfft, const int Ncp, const int inNupsample = 1);
  //! Modulate complex data symbols. Length of \c input must be an integer multiple of \c Nfft
  cvec modulate(const cvec &input);
  //! Modulate complex data symbols. Length of \c input must be an integer multiple of \c Nfft
  void modulate(const cvec &input, cvec &output);
  //! Demodulate to complex valued symbols. Length of \c input must be an integer multiple of \c Nfft+Ncp
  cvec demodulate(const cvec &input);
  //! Demodulate to complex valued symbols. Length of \c input must be an integer multiple of \c Nfft+Ncp
  void demodulate(const cvec &input, cvec &output);
private:
  double norm_factor;
  bool setup_done;
  int Nfft, Ncp, Nupsample;
};

} // namespace itpp

#endif // #ifndef OFDM_H
