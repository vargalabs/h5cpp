/*!
 * \file
 * \brief Definition of a CRC code class
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

#ifndef CRC_H
#define CRC_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief Cyclic Redundancy Check Codes

  This class will add the CRC bits after each input word. With b(i)
denoting the i-th input bit and p(i) the i-th parity check bit,
the order of the outbut bits will be:
  \code [b(1), b(2), ..., b(k), p(1), p(2), ..., p(n-k)] \endcode

  When the WCDMA CRC polynomials are used, this class will reverse
  the order of the parity check bits in order to comply to the WCDMA
  standard. Thus for the polynomials WCDMA-8, WCDMA-12, WCDMA-16,
  and WCDMA-24 the output will be:
  \code [b(1), b(2), ..., b(k), p(n-k), ..., p(2), p(1)] \endcode

  Usage:
  \code
  CRC_Code crc(string("CRC-4"));
  bvec bits = randb(10), coded_bits, decoded_bits;
  bool error;

  coded_bits = crc.encode(bits);
  error = crc.decode(rec_bits, decoded_bits);
  \endcode
*/
class ITPP_EXPORT CRC_Code
{
public:

  //! Default Constructor
  CRC_Code() { reverse_parity = false; }

  /*!
    \brief Set CRC code to one of the standardpolynomials using the
  string value.
    \param code Possible values: CRC-4, CRC-7, CRC-8, CRC-12,
    CRC-24, CRC-32, CCITT-4, CCITT-5, CCITT-6, CCITT-16, CCITT-32,
    WCDMA-8, WCDMA-12, WCDMA-16, WCDMA-24, ATM-8, ANSI-16, SDLC-16
  */
  CRC_Code(const std::string &code) { reverse_parity = false; set_code(code); }

  //! Set an arbitary polynomial in bvec form. Start with highest order terms.
  void set_generator(const bvec &poly);

  //! Set CRC code to one of the standardpolynomials using the string value.
  void set_code(const std::string &code);

  //! Calulate the parity bits
  void parity(const bvec &in_bits, bvec &out) const;

  //! Return true if parity checks OK otherwise flase
  bool check_parity(const bvec &coded_bits) const;

  //! Calculate and add parity to the in_bits.
  void encode(const bvec &in_bits, bvec &out) const;

  //! Returns the in_bits vector with parity added
  bvec encode(const bvec &in_bits) const;

  //! Return true if parity checks OK otherwise flase. Also returns the message part in out.
  bool decode(const bvec &coded_bits, bvec &out) const;

  //! Return true if parity checks OK otherwise flase. Also returns the message part in bits.
  bool decode(bvec &bits) const;

private:
  bool reverse_parity;
  bvec polynomial;
  int no_parity;
};

} // namespace itpp

#endif // #ifndef CRC_H
