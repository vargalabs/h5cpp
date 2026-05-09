/*!
 * \file
 * \brief Binary class definition
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

#ifndef BINARY_H
#define BINARY_H

#include <itpp/base/itassert.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \brief Binary arithmetic (boolean) class
  \author Tony Ottosson

  This class creates a binary arithmetic class, following the ordinary
  rules for binary (GF(2)) fields.

  Examples:
  \code
  bin a;         // Creation of variable
  bin a = 0;     // Creating a variable and assigning it value 0
  bin b = 1;     // Creating a variable and assigning it value 1
  bin c = a + b; // XOR operation
  c = !a;        // NOT
  c = a * b;     // AND
  c = a / b;     // OR
  \endcode
*/
class bin
{
public:
  //! Default constructor
  bin(): b(0) {}

  //! Set the binary object equal to \c value. Either "0" or "1".
  bin(const int &value): b(static_cast<char>(value)) {
    it_assert_debug((value == 0) || (value == 1),
                    "bin::bin(): value must be 0 or 1");
  }

  //! Copy constructor
  bin(const bin &inbin): b(inbin.b) {}

  //! Assign a value
  void operator=(const int &value) {
    it_assert_debug((value == 0) || (value == 1),
                    "bin::operator=(): value must be 0 or 1");
    b = static_cast<char>(value);
  }

  //! Assign a value
  void operator=(const bin &inbin) { b = inbin.b; }

  //! OR
  void operator/=(const bin &inbin) { b |= inbin.b; }

  //! OR
  void operator|=(const bin &inbin) { b |= inbin.b; }
  //! OR
  bin operator/(const bin &inbin) const { return bin(b | inbin.b); }
  //! OR
  bin operator|(const bin &inbin) const { return bin(b | inbin.b); }

  //! XOR
  void operator+=(const bin &inbin) { b ^= inbin.b; }
  //! XOR
  void operator^=(const bin &inbin) { b ^= inbin.b; }
  //! XOR
  bin operator+(const bin &inbin) const { return bin(b ^ inbin.b); }
  //! XOR
  bin operator^(const bin &inbin) const { return bin(b ^ inbin.b); }
  //! XOR
  void operator-=(const bin &inbin) { b ^= inbin.b; }
  //! XOR
  bin operator-(const bin &inbin) const { return bin(b ^ inbin.b); }
  //! Dummy definition to be able to use vec<bin>
  bin operator-() const { return bin(b); }

  //! AND
  void operator*=(const bin &inbin) { b &= inbin.b; }
  //! AND
  void operator&=(const bin &inbin) { b &= inbin.b; }
  //! AND
  bin operator*(const bin &inbin) const { return bin(b & inbin.b); }
  //! AND
  bin operator&(const bin &inbin) const { return bin(b & inbin.b); }

  //! NOT
  bin operator!(void) const { return bin(b ^ 1); }
  //! NOT
  bin operator~(void) const { return bin(b ^ 1); }

  //! Check if equal
  bool operator==(const bin &inbin) const { return b == inbin.b; }
  //! Check if equal
  bool operator==(const int &i) const { return b == i; }

  //! Check if not equal
  bool operator!=(const bin &inbin) const { return b != inbin.b; }
  //! Check if not equal
  bool operator!=(const int &i) const { return b != i; }

  //! Less than (interpret the binary values {0,1} as integers)
  bool operator<(const bin &inbin) const  { return b < inbin.b; }
  //! Less than equal (interpret the binary values {0,1} as integers)
  bool operator<=(const bin &inbin) const { return b <= inbin.b; }

  //! Greater than (interpret the binary values {0,1} as integers)
  bool operator>(const bin &inbin) const  { return b > inbin.b; }
  //! Greater than equal (interpret the binary values {0,1} as integers)
  bool operator>=(const bin &inbin) const { return b >= inbin.b; }

  //! Convert \c bin to \c short
  operator short() const  { return static_cast<short>(b); }
  //! Convert \c bin to \c int
  operator int() const    { return static_cast<int>(b); }
  //! Convert \c bin to \c bool
  operator bool() const   { return b != 0; }
  //! Convert \c bin to \c float
  operator float() const  { return static_cast<float>(b); }
  //! Convert \c bin to \c double
  operator double() const { return static_cast<double>(b); }

  //! Output the binary value of the object
  char value() const { return b; }

private:
  char b;
};

/*!
  \relatesalso bin
  \brief Output stream of bin
*/
ITPP_EXPORT std::ostream &operator<<(std::ostream &output, const bin &inbin);

/*!
  \relatesalso bin
  \brief Input stream of bin
*/
ITPP_EXPORT std::istream &operator>>(std::istream &input, bin &outbin);

/*!
  \relatesalso bin
  \brief absolute value of bin
*/
inline bin abs(const bin &inbin) { return inbin; }

} // namespace itpp


namespace std   // added 11/2005, EGL
{

/*!
  \relatesalso itpp::bin
  \brief absolute value of bin
*/
inline int abs(const itpp::bin &inbin) { return inbin; }

} // namespace std

#endif // #ifndef BINARY_H

