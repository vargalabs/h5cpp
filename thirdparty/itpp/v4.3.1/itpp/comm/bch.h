/*!
 * \file
 * \brief Definition of a BCH encoder/decoder class
 * \author Pal Frenger, Steve Peters, Adam Piatyszek and Stephan Ludwig
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

#ifndef BCH_H
#define BCH_H

#include <itpp/comm/galois.h>
#include <itpp/comm/channel_code.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup fec
*/

//---------------------- BCH --------------------------------------

/*!
  \ingroup fec
  \brief Class for binary, narrow-sense BCH codes.

  The notation used is found in S. B. Wicker, "Error control systems for
  digital communication and storage", Appendix E, Prentice-Hall, 1995.

  Example:
  \code BCH bch(31,21,2,ivec("3 5 5 1"))
  \endcode
  uses the generator polynomial
  \f$g(x) = x^{10} + x^9 + x^8 + x^6 + x^5 + x^3 + 1\f$, and is capable of
  correcting 2 errors with \a n = 31 and \a k = 21.
*/
class ITPP_EXPORT BCH : public Channel_Code
{
public:
  /*!
   * \brief Initialize a (n,k)-code that can correct t errors
   *
   * \note Do not call this constructor with e.g. BCH bch(31, 21, 2, "3 5 5 1")
   * but with BCH bch(31, 21, 2, ivec("3 5 5 1")) instead. Otherwise the
   * complier will complain.
   */
  BCH(int in_n, int in_k, int in_t, const ivec &genpolynom, bool sys = false);

  /*!
   * \brief Initialize a (n,k)-code that can correct t errors
   * \author Stephan Ludwig
   *
   * The generator polynomial is automatically generated from the (n, t)
   * parameters of the BCH code. The constructor generates the generator
   * polynomial (and determines k) according to the method described in:
   *
   * [Wic95] S. B. Wicker, "Error control systems for digital communication
   * and storage", Prentice-Hall, 1995
   */
  BCH(int in_n, int in_t, bool sys = false);

  //! Destructor
  virtual ~BCH() { }

  //! Encode a bvec of indata
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
  //! Encode a bvec of indata
  virtual bvec encode(const bvec &uncoded_bits);

  /*!
  * \brief Decode the BCH code bits. Return false if there has been any decoding failure.
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

  //! Decode a bvec of coded data. This function is kept for backward compatibility. Better use \code bool decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid) \endcode.
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);

  //! Decode a bvec of coded data. This function is kept for backward compatibility. Better use \code bool decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid)\endcode.
  virtual bvec decode(const bvec &coded_bits);

  // Soft-decision decoding is not implemented
  virtual void decode(const vec &received_signal, bvec &output);
  virtual bvec decode(const vec &received_signal);

  //! Get the code rate
  virtual double get_rate() const {
    return static_cast<double>(k) / n;
  }

  //! Get cardinality of code k
  virtual int get_k() const {
    return k;
  }

  //! Dummy assignment operator - MSVC++ warning C4512
  BCH & operator= (const BCH &) {
    return *this;
  }

private:
  int n, k, t;
  GFX g;
  const bool systematic;
};

} // namespace itpp

#endif // #ifndef BCH_H
