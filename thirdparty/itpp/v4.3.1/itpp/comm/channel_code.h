/*!
 * \file
 * \brief Channel Code class virtual interface
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

#ifndef CHANNEL_CODE_H
#define CHANNEL_CODE_H

#include <itpp/base/vec.h>
#include <itpp/comm/modulator.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup fec
*/

//---------------------- BCH --------------------------------------

/*!
  \ingroup fec
  \brief Generic Channel Code class

*/
class ITPP_EXPORT Channel_Code
{
public:
  //! Default constructor
  Channel_Code() {}
  //! Destructor
  virtual ~Channel_Code() {}

  //! Encode a bvec of input
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits) = 0;
  //! Encode a bvec of input
  virtual bvec encode(const bvec &uncoded_bits) = 0;

  //! Decode a bvec of coded data.
  virtual void decode(const bvec &codedbits, bvec &decoded_bits) = 0;
  //! Decode a bvec of coded data
  //virtual bvec decode(const bvec &coded_bits);
  virtual bvec decode(const bvec &coded_bits) = 0;

  //! Decode a vec of received data.
  virtual void decode(const vec &received_signal, bvec &decoded_bits) = 0;
  //! Decode a vec of received data
  //virtual bvec decode(const vec &received_signal);
  virtual bvec decode(const vec &received_signal) = 0;

  //! Get the code rate
  virtual double get_rate() const = 0;
};


/*!
  \ingroup fec
  \brief Dummy Channel Code class

  A dummy code class. Uncoded output.
*/
class ITPP_EXPORT Dummy_Code : public Channel_Code
{
public:
  //! Default constructor
  Dummy_Code() {}
  //! Destructor
  virtual ~Dummy_Code() {}

  //! Encode a bvec of input
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits) { coded_bits = uncoded_bits; }
  //! Encode a bvec of input
  virtual bvec encode(const bvec &uncoded_bits) { return uncoded_bits; }

  //! Decode a bvec of coded data.
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits) { decoded_bits = coded_bits; }
  //! Decode a bvec of coded data.
  virtual bvec decode(const bvec &coded_bits) { return coded_bits; }

  //! Decode a vec of received data. Assumes soft input (BPSK modulated)
  virtual void decode(const vec &received_signal, bvec &decoded_bits) { BPSK bpsk; bpsk.demodulate_bits(received_signal, decoded_bits); }
  //! Decode a vec of received data. Assumes soft input (BPSK modulated)
  virtual bvec decode(const vec &received_signal) { bvec out; decode(received_signal, out); return out; }

  //! Get the code rate
  virtual double get_rate() const { return 1.0; }
};



} // namespace itpp

#endif // #ifndef CHANNEL_CODE_H
