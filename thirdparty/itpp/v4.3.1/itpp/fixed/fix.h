/*!
 * \file
 * \brief Definitions of a fixed-point data type Fix
 * \author Johan Bergman
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

#ifndef FIX_H
#define FIX_H

#include <itpp/fixed/fix_base.h>
#include <itpp/fixed/fix_factory.h>
#include <itpp/itexports.h>

namespace itpp
{

// Forward declarations
template<class Num_T> class Vec;
template<class Num_T> class Mat;

//! \addtogroup fixed
//!@{

/*!
  \brief Fixed-point data type

  See the Detailed Description in the \ref fixed module.
*/
class ITPP_EXPORT Fix : public Fix_Base
{
  friend class CFix;
  template<int, e_mode, o_mode, q_mode> friend class Fixed;
  template<int, e_mode, o_mode, q_mode> friend class CFixed;
public:
  //! Default constructor
  Fix(double x = 0.0, int s = 0, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(s, w, e, o, q, ptr), re(scale_and_apply_modes(x)) {}
  //! Constructor
  explicit Fix(const Fix_Factory &f)
      : Fix_Base(0, f.wordlen, f.emode, f.omode, f.qmode, f.stat_ptr), re(0) {}
  //! Constructor for internal use. No restrictions are applied. The dummies help to avoid ambiguities
  Fix(fixrep r, int s, int, int)
      : Fix_Base(s), re(r) {}
  //! Copy constructor
  Fix(const Fix &x, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(x.shift, w, e, o, q, ptr), re(x.re) {}
  //! Destructor
  virtual ~Fix() {}

  //! Assignment from Fix
  Fix& operator=(const Fix &x);
  //! Assignment from int
  Fix& operator=(const int x);
  //! Addition of Fix
  Fix& operator+=(const Fix &x);
  //! Addition of int
  Fix& operator+=(const int x);
  //! Subtraction of Fix
  Fix& operator-=(const Fix &x);
  //! Subtraction of int
  Fix& operator-=(const int x);
  //! Multiplication with Fix
  Fix& operator*=(const Fix &x);
  //! Multiplication with int
  Fix& operator*=(const int x);
  //! Division with Fix using quantization mode \c TRN
  Fix& operator/=(const Fix &x);
  //! Division with int using quantization mode \c TRN
  Fix& operator/=(const int x);
  //! Unary negative of Fix
  Fix operator-() const;
  //! Left shift \c n bits
  Fix& operator<<=(const int n);
  //! Right shift \c n bits using quantization mode \c qmode (constructor argument)
  Fix& operator>>=(const int n);

  //! Set to <tt>x * pow2(n)</tt> using quantization mode \c qmode (constructor argument)
  void set(double x, int n);
  //! Set to <tt>x * pow2(n)</tt> using quantization mode \c q (function argument)
  void set(double x, int n, q_mode q);
  //! Set data representation (mainly for internal use since it reveals the representation type)
  void set_re(fixrep x) {re = apply_o_mode(x);}

  //! Left shift \c n bits
  void lshift(int n);
  //! Right shift \c n bits using quantization mode \c qmode (constructor argument)
  void rshift(int n);
  //! Right shift \c n bits using quantization mode \c q (function argument)
  void rshift(int n, q_mode q);

  //! Print restrictions
  virtual void print() const;
  //! Get data representation (mainly for internal use since it reveals the representation type)
  fixrep get_re() const {return re;}
  //! Conversion to double
  double unfix() const;

#ifndef NO_IMPLICIT_FIX_CONVERSION
  //! Conversion to double
  operator double() const {
    it_assert_debug(shift>=-63 && shift <= 64, "Fix::operator double: Illegal shift!");
    return double(re)*DOUBLE_POW2[64 - shift];
  }
#endif

  //! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
  friend ITPP_EXPORT int assert_shifts(const CFix &x, const Fix &y);
  //! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
  friend ITPP_EXPORT int assert_shifts(const Fix &x, const Fix &y);
  //! Check that x.shift==0 OR x==0 OR y==0 and return x.shift
  friend ITPP_EXPORT int assert_shifts(const Fix &x, int y);

protected:
  //! Data representation
  fixrep re;
};

//! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
ITPP_EXPORT int assert_shifts(const Fix &x, const Fix &y);
//! Check that x.shift==0 OR x==0 OR y==0 and return x.shift
ITPP_EXPORT int assert_shifts(const Fix &x, int y);

//! Input bit representation and, optionally, the shift
ITPP_EXPORT std::istream &operator>>(std::istream &is, Fix &x);
//! Output bit representation and, optionally, the shift
ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const Fix &x);

//! Typedef for fixed-point vector type
typedef Vec<Fix> fixvec;
//! Typedef for fixed-point matrix type
typedef Mat<Fix> fixmat;

// Specialization of template definition in vec.cpp
template<> void fixvec::set(const char *values);
// Specialization of template definition in mat.cpp
template<> void fixmat::set(const char *values);

//!@}

} // namespace itpp

#endif // #ifndef FIX_H
