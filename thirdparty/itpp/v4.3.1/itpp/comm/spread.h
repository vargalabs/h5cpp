/*!
 * \file
 * \brief Definition of spread spectrum classes and functions
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

#ifndef SPREAD_H
#define SPREAD_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup modulators
  \brief Spreading of float symbols to float output

  Spreading block for generation of 1-dimensional DS-CDMA signals
  Useful in the simulation of DS-CDMA systems on chip level or upsampled chip level.

  Obeserve that the spreading is
  normalized so that the energy per bit is preserved before and after spreading
  (that is each symbol is multiplied with \c 1/sqrt(N)).
  Hence, for the multicode case the energy is normalized for one symbol
  (code) but the transmitted signal consist of a sum of several signals.

  Four different classes exist:

  <ul>
  <li> Spread_1d for 1-dimensional symbols (e.g. BPSK).</li>
  <li> Spread_2d for 2-dimensional symbols (e.g. QPSK).</li>
  </ul>

  For multicode transmission, that is each user uses several codes in parallel to transmit data, there exist

  <ul>
  <li> Multicode_Spread_1d for 1-dimensional symbols. </li>
  <li> Multicode_Spread_2d for 2-dimensional symbols. </li>
  </ul>

  Example:
  \code
  #include "itpp/itcomm.h"

  int main() {

  //Generate the spreading code you want to use.
  vec spreading_code = "-1 1 -1 1";

  //Initiate th Spreading class
  Spread_1d spread_1d(spreading_code);

  //Generate the symbols to transmitt
  bvec transmitted_bits = randb(10);
  BPSK bpsk;
  vec transmitted_symbols = bpsk.modulate_bits(transmitted_bits);

  //Spread the symbols
  vec transmitted_signal = spread_1d.spread(transmitted_symbols);

  //Generate the received signal
  vec received_signal = transmitted_signal;

  //Despread the received signal
  vec received_symbols  = spread_1d.despread(received_signal,0);

  //demodulate the bits
  bvec received_bits = bpsk.demodulate_bits(received_symbols);

  }
  \endcode

*/
class ITPP_EXPORT Spread_1d
{
public:
  //! Constructor
  Spread_1d() { }
  //! Constructor
  Spread_1d(const vec &incode);
  //! Spreading of signal return i out.
  void spread(const vec &symbols, vec &out);
  //! Spreading of signal.
  vec spread(const vec &symbols) { vec out; spread(symbols, out); return out; }
  /*!
    \brief Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  */
  void despread(const vec &rec_signal, vec &out, int timing);
  /*!
    \brief Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  */
  vec despread(const vec &rec_signal, int timing)
  { vec out; despread(rec_signal, out, timing); return out; }
  //! Set the spreading code used for spreading
  void set_code(const vec &incode);
  //! Returns the spreading code used
  vec get_code();
  //! Get the period of the code (length of code vector).
  int get_period() { return N; }
protected:
  //! The spreading code
  vec code;
  //! The spreading factor
  int N;
};

/*!
  \ingroup modulators
  \brief Spreading of complex symbols to complex output.

  The spreading are done independently for the I and Q phases. That is
  real(symbols) are spread by the incodeI and imag(symbols) are spread
  by incodeQ.

  Before despreading the phase should be corrected, that is the
  complex baseband signal should be multiplied by exp(j*PHIk), where
  PHIk is the phase of that user (and path).

  Obeserve that the spreading is
  normalized so that the energy per bit is preserved before and after spreading
  (that is each symbol is multiplied with \c 1/sqrt(N)).
  Hence, for the multicode case the energy is normalized for one symbol
  (code) but the transmitted signal consist of a sum of several signals.

  Example: See Spread_1d
*/
class ITPP_EXPORT Spread_2d
{
public:
  //! Constructor
  Spread_2d() { }
  //! Constructor
  Spread_2d(const vec &incodeI, const vec &incodeQ);
  //! Spreading of signal
  void spread(const cvec &symbols, cvec &out);
  //! Spreading of signal
  cvec spread(const cvec &symbols) { cvec out; spread(symbols, out); return out; }
  /*!
    \brief Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  */
  void despread(const cvec &rec_signal, cvec &out, int timing);
  /*!
    \brief Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  */
  cvec despread(const cvec &rec_signal, int timing)
  { cvec out; despread(rec_signal, out, timing); return out; }
  //! Set the in-phase and the quadrature-phase spreading codes
  void set_code(const vec &incodeI, const vec &incodeQ);
  //! Returns the in-phase spreading code
  vec get_codeI();
  //! Returns the quadrature-phase spreading code
  vec get_codeQ();
  //! Get the period of the code (length of code vector).
  int get_period() { return spreadI.get_period(); }
protected:
  /*! The spreaders for the I and Q channels respectively
   * @{ */
  Spread_1d spreadI, spreadQ;
  /*! @} */
};

/*!
  \ingroup modulators
  \brief Multicode spreading of float symbols

  Obeserve that the spreading is
  normalized so that the energy per bit is preserved before and after spreading
  (that is each symbol is multiplied with \c 1/sqrt(N)).
  Hence, for the multicode case the energy is normalized for one symbol
  (code) but the transmitted signal consist of a sum of several signals.

  Example: See Spread_1d
*/
class ITPP_EXPORT Multicode_Spread_1d
{
public:
  //! Constructor
  Multicode_Spread_1d() { }
  //! Constructor
  Multicode_Spread_1d(const mat &incodes);
  //! Spreading function
  vec spread(const vec &symbols);
  //! Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  vec despread(const vec &receivedsignal, int timing);
  //! Set the spreading codes. Each row represent one spreading code. The spreading factor equals the number of columns
  void set_codes(const mat &incodes);
  //! Returns the matrix containing the spreading codes used as rows in the matrix
  mat get_codes();
  //! Returns the spreading factor
  int get_period() { return N; }
  //! Returns the number of multi-codes used
  int get_nocodes() { return L; }
protected:
  //! The spreading codes used size (\f$L \times N\f$)
  mat codes;
  //! The number of multi-codes
  int L;
  //! The spreading factor
  int N;
};

/*!
  \ingroup modulators
  \brief Multicode spreading of complex symbols to complex output.

  The spreading are done independently for the I and Q phases. That is
  real(symbols) are spread by the incodeI and imag(symbols) are spread
  by incodeQ.

  Before despreading the phase should be corrected, that is the
  complex baseband signal should be multiplied by exp(j*PHIk), where
  PHIk is the phase of that user (and path).

  Obeserve that the spreading is
  normalized so that the energy per bit is preserved before and after spreading
  (that is each symbol is multiplied with \c 1/sqrt(N)).
  Hence, for the multicode case the energy is normalized for one symbol
  (code) but the transmitted signal consist of a sum of several signals.

  Example: See Spread_1d
*/
class ITPP_EXPORT Multicode_Spread_2d
{
public:
  //! Constructor
  Multicode_Spread_2d() { }
  //! Constructor
  Multicode_Spread_2d(const mat &incodesI, const mat &incodesQ);
  //! Spreading of signal
  cvec spread(const cvec &symbols);
  //! Despreading of signal. \a timing is the start position of the first symbol, given in number of samples.
  cvec despread(const cvec &receivedsignal, int timing);
  /*!
    \brief Set the spreading codes

    The codes are given as rows in the matricies \a incodesI and \a incodesQ. The number of rows
    shall equal the number of multiple spreading codes
  */
  void set_codes(const mat &incodesI, const mat &incodesQ);
  //! Return the matrix containing the in-phase codes (as rows)
  mat get_codesI();
  //! Return the matrix containing the quadrature-phase codes (as rows)
  mat get_codesQ();
  //! Returns the spreading factor
  int get_period() { return mcspreadI.get_period(); }
protected:
  /*! The multicode spreaders for the I and Q channels respectively
   * @{ */
  Multicode_Spread_1d mcspreadI, mcspreadQ;
  /*! @} */
};

} // namespace itpp

#endif // #ifndef SPREAD_H
