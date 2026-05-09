/*!
 * \file
 * \brief Definitions of a Reed-Solomon codec class
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
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

#ifndef REEDSOLOMON_H
#define REEDSOLOMON_H

#include <itpp/base/vec.h>
#include <itpp/comm/galois.h>
#include <itpp/comm/channel_code.h>
#include <itpp/itexports.h>

namespace itpp
{

//---------------------- Reed-Solomon --------------------------------------

/*!
  \ingroup fec
  \brief Reed-Solomon Codes.

  Uses the Berlkamp-Massey algorithm for decoding as described in: S. B. Wicker,
  "Error Control Systems for digital communication and storage," Prentice Hall.

  The code is \f$2^m\f$ - ary of length \f$2^m-1\f$ capable of correcting \f$t\f$ errors.
*/
class ITPP_EXPORT Reed_Solomon : public Channel_Code
{
public:
  /*! Class constructor for the \f$2^m\f$ - ary, \f$t\f$ error correcting RS-code.
	* The generator polynomial will be $g(x)=\prod_{i=0}^{2t-1}(x-\alpha^{b+i})$,
	* where $\alpha$ is a root of the primitive polynomial of \c itpp::GF.
	*/
  Reed_Solomon(int in_m, int in_t, bool sys = false, int in_b = 1);
  //! Destructor
  virtual ~Reed_Solomon() { }

  //! Encoder function.
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
  //! Encoder function.
  virtual bvec encode(const bvec &uncoded_bits);

  /*!
   * \brief Decode the RS code bits. Return false if there has been any decoding failure.
   *
   * The \c coded_bits are block-wise decoded into \c decoded_message messages. If
   * there has been any decoding failure, the function will return \c false.
   * If this happened in the n-th block, then \c cw_isvalid(n) will be set
   * to \c false (zero-indexed). In case of a systematic code the systematic bits will be
   * extracted and presented in the corresponding block of \c decoded_message.
   * This is better than just presenting zeros, which is done in case of a
   * decoding failure of non-systematic codes.
   * 
   * For erasure decoding the indices of erased positions need to be passed in \c erasure_positions 
   * as indices to the erased \fsymbol\f (not bit!).
   * 
   */
  virtual bool decode(const bvec &coded_bits, const ivec &erasure_positions, bvec &decoded_message, bvec &cw_isvalid);

  /*!
   * \brief Decode the RS code bits. Return false if there has been any decoding failure.
   *
   * The \c coded_bits are block-wise decoded into \c decoded_message messages. If
   * there has been any decoding failure, the function will return \c false.
   * If this happened in the n-th block, then \c cw_isvalid(n) will be set
   * to \c false (zero-indexed). In case of a systematic code the systematic bits will be
   * extracted and presented in the corresponding block of \c decoded_message.
   * This is better than just presenting zeros, which is done in case of a
   * decoding failure of non-systematic codes.
   */
  virtual bool decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid);

  //! Decoder a \c bvec of coded data. This function is kept for backward compatibility. Better use \code bool decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid)\endcode.
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
  //! Decoder a \c bvec of coded data. This function is kept for backward compatibility. Better use \code bool decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid)\endcode.
  virtual bvec decode(const bvec &coded_bits);

  // Soft-decision decoding is not implemented
  virtual void decode(const vec &received_signal, bvec &output);
  virtual bvec decode(const vec &received_signal);

  //! Gets the rate of the RS-code.
  virtual double get_rate() const { return static_cast<double>(k) / n; }

  //! Dummy assignment operator - MSVC++ warning C4512
  Reed_Solomon & operator=(const Reed_Solomon &) { return *this; }

protected:
  /*! Internal encoder/decoder parameters
   * @{ */
  int m, t, k, n, q, b;
  /*! @} */
  //! The generator polynomial of the RS code
  GFX g;
  //! Whether or not the code is systematic
  const bool systematic;
};

} // namespace itpp

#endif // #ifndef REEDSOLOMON_H
