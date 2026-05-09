/*!
 * \file
 * \brief Definitions of a complex fixed-point data type CFix
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

#ifndef CFIX_H
#define CFIX_H

#include <itpp/fixed/fix.h>
#include <itpp/itexports.h>


namespace itpp
{

// Forward declarations
template<class Num_T> class Vec;
template<class Num_T> class Mat;

//! \addtogroup fixed
//!@{

/*!
  \brief Complex fixed-point data type.

  See the Detailed Description in the \ref fixed module.
*/
class ITPP_EXPORT CFix : public Fix_Base
{
  template<int, e_mode, o_mode, q_mode> friend class CFixed;
public:
  //! Default constructor
  CFix(double r = 0.0, double i = 0.0, int s = 0, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(s, w, e, o, q, ptr), re(scale_and_apply_modes(r)), im(scale_and_apply_modes(i)) {}
  //! Constructor
  CFix(std::complex<double> x, double, int s = 0, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(s, w, e, o, q, ptr), re(scale_and_apply_modes(std::real(x))), im(scale_and_apply_modes(std::imag(x))) {}
  //! Constructor
  explicit CFix(const Fix_Factory &f)
      : Fix_Base(0, f.wordlen, f.emode, f.omode, f.qmode, f.stat_ptr), re(0), im(0) {}
  //! Constructor for internal use. No restrictions are applied. The dummies help to avoid ambiguities
  CFix(fixrep r, fixrep i, int s, int, int)
      : Fix_Base(s), re(r), im(i) {}
  //! Constructor
  CFix(const Fix &r, const Fix &i = 0.0, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(assert_shifts(r, i), w, e, o, q, ptr), re(r.re), im(i.re) {}
  //! Copy constructor
  CFix(const CFix &x, double, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : Fix_Base(x.shift, w, e, o, q, ptr), re(x.re), im(x.im) {}
  //! Destructor
  virtual ~CFix() {}

  //! Assignment from CFix
  CFix& operator=(const CFix &x);
  //! Assignment from Fix
  CFix& operator=(const Fix &x);
  //! Assignment from std::complex<double>. Fractional part is truncated
  CFix& operator=(const std::complex<double> &x);
  //! Assignment from int
  CFix& operator=(const int x);
  //! Addition of CFix
  CFix& operator+=(const CFix &x);
  //! Addition of Fix
  CFix& operator+=(const Fix &x);
  //! Addition of int
  CFix& operator+=(const int x);
  //! Subtraction of CFix
  CFix& operator-=(const CFix &x);
  //! Subtraction of Fix
  CFix& operator-=(const Fix &x);
  //! Subtraction of int
  CFix& operator-=(const int x);
  //! Multiplication with CFix. Temporary variables use the maximum word length (64 bits)
  CFix& operator*=(const CFix &x);
  //! Multiplication with Fix. Temporary variables use the maximum word length (64 bits)
  CFix& operator*=(const Fix &x);
  //! Multiplication with int. Temporary variables use the maximum word length (64 bits)
  CFix& operator*=(const int x);
  //! Division with CFix using quantization mode \c TRN. Temporary variables use the maximum word length (64 bits)
  CFix& operator/=(const CFix &x);
  //! Division with Fix using quantization mode \c TRN. Temporary variables use the maximum word length (64 bits)
  CFix& operator/=(const Fix &x);
  //! Division with int using quantization mode \c TRN. Temporary variables use the maximum word length (64 bits)
  CFix& operator/=(const int x);
  //! Unary negative of CFix
  CFix operator-() const;
  //! Left shift \c n bits
  CFix& operator<<=(const int n);
  //! Right shift \c n bits using quantization mode \c qmode (constructor argument)
  CFix& operator>>=(const int n);

  //! Set to <tt>(real + i*imag) * pow2(n)</tt> using quantization mode \c qmode (constructor argument)
  void set(double real, double imag, int n);
  //! Set to <tt>(real + i*imag) * pow2(n)</tt> using quantization mode \c q (function argument)
  void set(double real, double imag, int n, q_mode q);
  //! Set to <tt>x * pow2(n)</tt> using quantization mode \c qmode (constructor argument)
  void set(const std::complex<double> &x, int n);
  //! Set to <tt>x * pow2(n)</tt> using quantization mode \c q (function argument)
  void set(const std::complex<double> &x, int n, q_mode q);
  //! Set data representation for real part (mainly for internal use since it reveals the representation type)
  void set_re(fixrep x) {re = apply_o_mode(x);}
  //! Set data representation for imaginary part (mainly for internal use since it reveals the representation type)
  void set_im(fixrep x) {im = apply_o_mode(x);}

  //! Left shift \c n bits
  void lshift(int n);
  //! Right shift \c n bits using quantization mode \c qmode (constructor argument)
  void rshift(int n);
  //! Right shift \c n bits using quantization mode \c q (function argument)
  void rshift(int n, q_mode q);

  //! Print restrictions
  virtual void print() const;
  //! Get data representation for real part (mainly for internal use since it reveals the representation type)
  fixrep get_re() const {return re;}
  //! Get data representation for imaginary part (mainly for internal use since it reveals the representation type)
  fixrep get_im() const {return im;}
  //! Conversion to std::complex<double>
  std::complex<double> unfix() const;

#ifndef NO_IMPLICIT_FIX_CONVERSION
  //! Conversion to std::complex<double>
  operator std::complex<double>() const {
    it_assert_debug(shift >= -63 && shift <= 64, "CFix::operator complex<double>: Illegal shift!");
    return std::complex<double>(double(re)*DOUBLE_POW2[64 - shift],
                                double(im)*DOUBLE_POW2[64 - shift]);
  }
#endif

  //! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
  friend ITPP_EXPORT int assert_shifts(const CFix &x, const CFix &y);
  //! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
  friend ITPP_EXPORT int assert_shifts(const CFix &x, const Fix &y);
  //! Check that x.shift==0 OR x==0 OR y==0 and return x.shift
  friend ITPP_EXPORT int assert_shifts(const CFix &x, int y);

protected:
  fixrep re;   //!< Real data part
  fixrep im;   //!< Imaginary data part
};

//! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
ITPP_EXPORT int assert_shifts(const CFix &x, const CFix &y);
//! Check that x.shift==y.shift OR x==0 OR y==0 and return the shift (for the non-zero argument)
ITPP_EXPORT int assert_shifts(const CFix &x, const Fix &y);
//! Check that x.shift==0 OR x==0 OR y==0 and return x.shift
ITPP_EXPORT int assert_shifts(const CFix &x, int y);

//! Input bit representation and, optionally, the shift
ITPP_EXPORT std::istream &operator>>(std::istream &is, CFix &x);
//! Output bit representation and, optionally, the shift
ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const CFix &x);

//! Typedef for complex fixed-point vector type
typedef Vec<CFix> cfixvec;
//! Typedef for complex fixed-point matrix type
typedef Mat<CFix> cfixmat;

// Specialization of template definition in vec.cpp
template<> void cfixvec::set(const char *values);
// Specialization of template definition in mat.cpp
template<> void cfixmat::set(const char *values);

//!@}

} // namespace itpp

#endif // #ifndef CFIX_H
