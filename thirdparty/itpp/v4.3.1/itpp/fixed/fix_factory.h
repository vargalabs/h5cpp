/*!
 * \file
 * \brief Definitions of a class factory for fixed-point data types Fix
 * and CFix
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

#ifndef FIX_FACTORY_H
#define FIX_FACTORY_H

#include <itpp/base/factory.h>
#include <itpp/fixed/fix_base.h>
#include <itpp/itexports.h>


namespace itpp
{

// Forward declarations
class Fix;
class CFix;

//! \addtogroup fixed
//!@{

/*!
  \brief Class factory for fixed-point data types Fix and CFix

  For an introduction to factories, see the Detailed Description for Factory.
  For more information on the fixed-point data types, see the Detailed
  Description in the \ref fixed module.

  This example shows how to declare a Fix_Factory:
  \code
  // Declare UFIX32, a factory for 32-bit unsigned Fix/CFix with wrap-around
  // i.e. a factory for Fix(0.0, 0, 32, US, WRAP) and CFix(0.0, 0, 32, US, WRAP)
  Fix_Factory UFIX32(32, US, WRAP);
  \endcode

  However, the user does not need to declare \c UFIX32 since it is one of the
  already declared factories in fix_factory.h (which is included by itbase.h):
  \code
  Fix_Factory FIX1(1, TC, WRAP);  // for Fix/CFix with 1 bit
  ...
  Fix_Factory FIX64(64, TC, WRAP);  // for Fix/CFix with 64 bits

  Fix_Factory UFIX1(1, US, WRAP);  // for Unsigned Fix/CFix with 1 bit
  ...
  Fix_Factory UFIX64(64, US, WRAP);  // for Unsigned Fix/CFix with 64 bits

  Fix_Factory SFIX1(1, TC, SAT);  // for Saturated Fix/CFix with 1 bit
  ...
  Fix_Factory SFIX64(64, TC, SAT);  // for Saturated Fix/CFix with 64 bits

  Fix_Factory SUFIX1(1, US, SAT);  // for Saturated Unsigned Fix/CFix with 1 bit
  ...
  Fix_Factory SUFIX64(64, US, SAT);  // for Saturated Unsigned Fix/CFix with 64 bits
  \endcode
  This means that it is only necessary for the user to declare a Fix_Factory if
  it is desired to have some other overflow mode than \c WRAP or \c SAT, or some
  other quantization mode than \c TRN, or a non-zero statistics object pointer.

  \note U stands for Unsigned but S stands for Saturated, NOT for Signed.

  The Array/Vec/Mat constructors can take a Fix_Factory as an argument:
  \code
  // Declare a Vec<Fix> with size 10 that will use
  // Fix(0.0, 0, 32, US, WRAP) for element creation
  Vec<Fix> vf(10, UFIX32);

  // Declare an Array<Array<Mat<CFix> > > with size 10 that will use
  // CFix(0.0, 0, 32, US, WRAP) for element creation
  Array<Array<Mat<CFix> > > aamcf(10, UFIX32);
  \endcode

  Even a Fix/CFix declaration can take a Fix_Factory as a constructor argument:
  \code
  // Equivalent to
  // Fix f(0.0, 0, 32, US, WRAP);
  Fix f(UFIX32);
  \endcode

  This syntax is also legal if Fix is replaced with \c double and CFix is
  replaced with <tt>complex<double></tt>, i.e.
  \code
  // The factory will be ignored
  Vec<double> vd(10, UFIX32);

  // The factory will be ignored
  Array<Array<Mat<complex<double> > > > aamcd(10, UFIX32);

  // The factory will be converted to double(0.0) i.e. innocent initialization
  double d(UFIX32);
  \endcode
  which can be useful in templated code, e.g. when the same code should support
  both floating- and fixed-point data types.
*/
class ITPP_EXPORT Fix_Factory : public Factory
{
  friend class Fix;
  friend class CFix;
public:
  //! Constructor
  explicit Fix_Factory(int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : wordlen(w), emode(e), omode(o), qmode(q), stat_ptr(ptr) {}
  //! Destructor
  virtual ~Fix_Factory() {}
  //! Conversion operator. Useful in templated code
  operator double() const {return 0.0;}
  //! Create an n-length array of Fix
  virtual void create(Fix* &ptr, const int n) const;
  //! Create an n-length array of CFix
  virtual void create(CFix* &ptr, const int n) const;
protected:
  //! Word length
  int wordlen;
  //! Sign encoding mode
  e_mode emode;
  //! Overflow mode
  o_mode omode;
  //! Quantization mode
  q_mode qmode;
  //! Pointer to statistics object
  Stat *stat_ptr;
};

//! Create an n-length array of Fix using Fix_Factory \c f
template<>
ITPP_EXPORT void create_elements<Fix>(Fix* &ptr, const int n, const Factory &f);

//! Create an n-length array of CFix using Fix_Factory \c f
template<>
ITPP_EXPORT void create_elements<CFix>(CFix* &ptr, const int n, const Factory &f);

//!@}


//! Fix_Factories for signed Fix/CFix with wrap-around (FIX1, FIX2, ..., FIX64)
const Fix_Factory FIX1(1, TC, WRAP);
//! \cond
const Fix_Factory FIX2(2, TC, WRAP);
const Fix_Factory FIX3(3, TC, WRAP);
const Fix_Factory FIX4(4, TC, WRAP);
const Fix_Factory FIX5(5, TC, WRAP);
const Fix_Factory FIX6(6, TC, WRAP);
const Fix_Factory FIX7(7, TC, WRAP);
const Fix_Factory FIX8(8, TC, WRAP);
const Fix_Factory FIX9(9, TC, WRAP);
const Fix_Factory FIX10(10, TC, WRAP);
const Fix_Factory FIX11(11, TC, WRAP);
const Fix_Factory FIX12(12, TC, WRAP);
const Fix_Factory FIX13(13, TC, WRAP);
const Fix_Factory FIX14(14, TC, WRAP);
const Fix_Factory FIX15(15, TC, WRAP);
const Fix_Factory FIX16(16, TC, WRAP);
const Fix_Factory FIX17(17, TC, WRAP);
const Fix_Factory FIX18(18, TC, WRAP);
const Fix_Factory FIX19(19, TC, WRAP);
const Fix_Factory FIX20(20, TC, WRAP);
const Fix_Factory FIX21(21, TC, WRAP);
const Fix_Factory FIX22(22, TC, WRAP);
const Fix_Factory FIX23(23, TC, WRAP);
const Fix_Factory FIX24(24, TC, WRAP);
const Fix_Factory FIX25(25, TC, WRAP);
const Fix_Factory FIX26(26, TC, WRAP);
const Fix_Factory FIX27(27, TC, WRAP);
const Fix_Factory FIX28(28, TC, WRAP);
const Fix_Factory FIX29(29, TC, WRAP);
const Fix_Factory FIX30(30, TC, WRAP);
const Fix_Factory FIX31(31, TC, WRAP);
const Fix_Factory FIX32(32, TC, WRAP);
const Fix_Factory FIX33(33, TC, WRAP);
const Fix_Factory FIX34(34, TC, WRAP);
const Fix_Factory FIX35(35, TC, WRAP);
const Fix_Factory FIX36(36, TC, WRAP);
const Fix_Factory FIX37(37, TC, WRAP);
const Fix_Factory FIX38(38, TC, WRAP);
const Fix_Factory FIX39(39, TC, WRAP);
const Fix_Factory FIX40(40, TC, WRAP);
const Fix_Factory FIX41(41, TC, WRAP);
const Fix_Factory FIX42(42, TC, WRAP);
const Fix_Factory FIX43(43, TC, WRAP);
const Fix_Factory FIX44(44, TC, WRAP);
const Fix_Factory FIX45(45, TC, WRAP);
const Fix_Factory FIX46(46, TC, WRAP);
const Fix_Factory FIX47(47, TC, WRAP);
const Fix_Factory FIX48(48, TC, WRAP);
const Fix_Factory FIX49(49, TC, WRAP);
const Fix_Factory FIX50(50, TC, WRAP);
const Fix_Factory FIX51(51, TC, WRAP);
const Fix_Factory FIX52(52, TC, WRAP);
const Fix_Factory FIX53(53, TC, WRAP);
const Fix_Factory FIX54(54, TC, WRAP);
const Fix_Factory FIX55(55, TC, WRAP);
const Fix_Factory FIX56(56, TC, WRAP);
const Fix_Factory FIX57(57, TC, WRAP);
const Fix_Factory FIX58(58, TC, WRAP);
const Fix_Factory FIX59(59, TC, WRAP);
const Fix_Factory FIX60(60, TC, WRAP);
const Fix_Factory FIX61(61, TC, WRAP);
const Fix_Factory FIX62(62, TC, WRAP);
const Fix_Factory FIX63(63, TC, WRAP);
const Fix_Factory FIX64(64, TC, WRAP);
//! \endcond

//! Fix_Factories for unsigned Fix/CFix with wrap-around (UFIX1, UFIX2, ..., UFIX64)
const Fix_Factory UFIX1(1, US, WRAP);
//! \cond
const Fix_Factory UFIX2(2, US, WRAP);
const Fix_Factory UFIX3(3, US, WRAP);
const Fix_Factory UFIX4(4, US, WRAP);
const Fix_Factory UFIX5(5, US, WRAP);
const Fix_Factory UFIX6(6, US, WRAP);
const Fix_Factory UFIX7(7, US, WRAP);
const Fix_Factory UFIX8(8, US, WRAP);
const Fix_Factory UFIX9(9, US, WRAP);
const Fix_Factory UFIX10(10, US, WRAP);
const Fix_Factory UFIX11(11, US, WRAP);
const Fix_Factory UFIX12(12, US, WRAP);
const Fix_Factory UFIX13(13, US, WRAP);
const Fix_Factory UFIX14(14, US, WRAP);
const Fix_Factory UFIX15(15, US, WRAP);
const Fix_Factory UFIX16(16, US, WRAP);
const Fix_Factory UFIX17(17, US, WRAP);
const Fix_Factory UFIX18(18, US, WRAP);
const Fix_Factory UFIX19(19, US, WRAP);
const Fix_Factory UFIX20(20, US, WRAP);
const Fix_Factory UFIX21(21, US, WRAP);
const Fix_Factory UFIX22(22, US, WRAP);
const Fix_Factory UFIX23(23, US, WRAP);
const Fix_Factory UFIX24(24, US, WRAP);
const Fix_Factory UFIX25(25, US, WRAP);
const Fix_Factory UFIX26(26, US, WRAP);
const Fix_Factory UFIX27(27, US, WRAP);
const Fix_Factory UFIX28(28, US, WRAP);
const Fix_Factory UFIX29(29, US, WRAP);
const Fix_Factory UFIX30(30, US, WRAP);
const Fix_Factory UFIX31(31, US, WRAP);
const Fix_Factory UFIX32(32, US, WRAP);
const Fix_Factory UFIX33(33, US, WRAP);
const Fix_Factory UFIX34(34, US, WRAP);
const Fix_Factory UFIX35(35, US, WRAP);
const Fix_Factory UFIX36(36, US, WRAP);
const Fix_Factory UFIX37(37, US, WRAP);
const Fix_Factory UFIX38(38, US, WRAP);
const Fix_Factory UFIX39(39, US, WRAP);
const Fix_Factory UFIX40(40, US, WRAP);
const Fix_Factory UFIX41(41, US, WRAP);
const Fix_Factory UFIX42(42, US, WRAP);
const Fix_Factory UFIX43(43, US, WRAP);
const Fix_Factory UFIX44(44, US, WRAP);
const Fix_Factory UFIX45(45, US, WRAP);
const Fix_Factory UFIX46(46, US, WRAP);
const Fix_Factory UFIX47(47, US, WRAP);
const Fix_Factory UFIX48(48, US, WRAP);
const Fix_Factory UFIX49(49, US, WRAP);
const Fix_Factory UFIX50(50, US, WRAP);
const Fix_Factory UFIX51(51, US, WRAP);
const Fix_Factory UFIX52(52, US, WRAP);
const Fix_Factory UFIX53(53, US, WRAP);
const Fix_Factory UFIX54(54, US, WRAP);
const Fix_Factory UFIX55(55, US, WRAP);
const Fix_Factory UFIX56(56, US, WRAP);
const Fix_Factory UFIX57(57, US, WRAP);
const Fix_Factory UFIX58(58, US, WRAP);
const Fix_Factory UFIX59(59, US, WRAP);
const Fix_Factory UFIX60(60, US, WRAP);
const Fix_Factory UFIX61(61, US, WRAP);
const Fix_Factory UFIX62(62, US, WRAP);
const Fix_Factory UFIX63(63, US, WRAP);
const Fix_Factory UFIX64(64, US, WRAP);
//! \endcond

//! Fix_Factories for unsigned Fix/CFix with wrap-around (SFIX1, SFIX2, ..., SFIX64)
const Fix_Factory SFIX1(1, TC, SAT);
//! \cond
const Fix_Factory SFIX2(2, TC, SAT);
const Fix_Factory SFIX3(3, TC, SAT);
const Fix_Factory SFIX4(4, TC, SAT);
const Fix_Factory SFIX5(5, TC, SAT);
const Fix_Factory SFIX6(6, TC, SAT);
const Fix_Factory SFIX7(7, TC, SAT);
const Fix_Factory SFIX8(8, TC, SAT);
const Fix_Factory SFIX9(9, TC, SAT);
const Fix_Factory SFIX10(10, TC, SAT);
const Fix_Factory SFIX11(11, TC, SAT);
const Fix_Factory SFIX12(12, TC, SAT);
const Fix_Factory SFIX13(13, TC, SAT);
const Fix_Factory SFIX14(14, TC, SAT);
const Fix_Factory SFIX15(15, TC, SAT);
const Fix_Factory SFIX16(16, TC, SAT);
const Fix_Factory SFIX17(17, TC, SAT);
const Fix_Factory SFIX18(18, TC, SAT);
const Fix_Factory SFIX19(19, TC, SAT);
const Fix_Factory SFIX20(20, TC, SAT);
const Fix_Factory SFIX21(21, TC, SAT);
const Fix_Factory SFIX22(22, TC, SAT);
const Fix_Factory SFIX23(23, TC, SAT);
const Fix_Factory SFIX24(24, TC, SAT);
const Fix_Factory SFIX25(25, TC, SAT);
const Fix_Factory SFIX26(26, TC, SAT);
const Fix_Factory SFIX27(27, TC, SAT);
const Fix_Factory SFIX28(28, TC, SAT);
const Fix_Factory SFIX29(29, TC, SAT);
const Fix_Factory SFIX30(30, TC, SAT);
const Fix_Factory SFIX31(31, TC, SAT);
const Fix_Factory SFIX32(32, TC, SAT);
const Fix_Factory SFIX33(33, TC, SAT);
const Fix_Factory SFIX34(34, TC, SAT);
const Fix_Factory SFIX35(35, TC, SAT);
const Fix_Factory SFIX36(36, TC, SAT);
const Fix_Factory SFIX37(37, TC, SAT);
const Fix_Factory SFIX38(38, TC, SAT);
const Fix_Factory SFIX39(39, TC, SAT);
const Fix_Factory SFIX40(40, TC, SAT);
const Fix_Factory SFIX41(41, TC, SAT);
const Fix_Factory SFIX42(42, TC, SAT);
const Fix_Factory SFIX43(43, TC, SAT);
const Fix_Factory SFIX44(44, TC, SAT);
const Fix_Factory SFIX45(45, TC, SAT);
const Fix_Factory SFIX46(46, TC, SAT);
const Fix_Factory SFIX47(47, TC, SAT);
const Fix_Factory SFIX48(48, TC, SAT);
const Fix_Factory SFIX49(49, TC, SAT);
const Fix_Factory SFIX50(50, TC, SAT);
const Fix_Factory SFIX51(51, TC, SAT);
const Fix_Factory SFIX52(52, TC, SAT);
const Fix_Factory SFIX53(53, TC, SAT);
const Fix_Factory SFIX54(54, TC, SAT);
const Fix_Factory SFIX55(55, TC, SAT);
const Fix_Factory SFIX56(56, TC, SAT);
const Fix_Factory SFIX57(57, TC, SAT);
const Fix_Factory SFIX58(58, TC, SAT);
const Fix_Factory SFIX59(59, TC, SAT);
const Fix_Factory SFIX60(60, TC, SAT);
const Fix_Factory SFIX61(61, TC, SAT);
const Fix_Factory SFIX62(62, TC, SAT);
const Fix_Factory SFIX63(63, TC, SAT);
const Fix_Factory SFIX64(64, TC, SAT);
//! \endcond

//! Fix_Factories for unsigned Fix/CFix with saturation (SUFIX1, SUFIX2, ..., SUFIX64)
const Fix_Factory SUFIX1(1, US, SAT);
//! \cond
const Fix_Factory SUFIX2(2, US, SAT);
const Fix_Factory SUFIX3(3, US, SAT);
const Fix_Factory SUFIX4(4, US, SAT);
const Fix_Factory SUFIX5(5, US, SAT);
const Fix_Factory SUFIX6(6, US, SAT);
const Fix_Factory SUFIX7(7, US, SAT);
const Fix_Factory SUFIX8(8, US, SAT);
const Fix_Factory SUFIX9(9, US, SAT);
const Fix_Factory SUFIX10(10, US, SAT);
const Fix_Factory SUFIX11(11, US, SAT);
const Fix_Factory SUFIX12(12, US, SAT);
const Fix_Factory SUFIX13(13, US, SAT);
const Fix_Factory SUFIX14(14, US, SAT);
const Fix_Factory SUFIX15(15, US, SAT);
const Fix_Factory SUFIX16(16, US, SAT);
const Fix_Factory SUFIX17(17, US, SAT);
const Fix_Factory SUFIX18(18, US, SAT);
const Fix_Factory SUFIX19(19, US, SAT);
const Fix_Factory SUFIX20(20, US, SAT);
const Fix_Factory SUFIX21(21, US, SAT);
const Fix_Factory SUFIX22(22, US, SAT);
const Fix_Factory SUFIX23(23, US, SAT);
const Fix_Factory SUFIX24(24, US, SAT);
const Fix_Factory SUFIX25(25, US, SAT);
const Fix_Factory SUFIX26(26, US, SAT);
const Fix_Factory SUFIX27(27, US, SAT);
const Fix_Factory SUFIX28(28, US, SAT);
const Fix_Factory SUFIX29(29, US, SAT);
const Fix_Factory SUFIX30(30, US, SAT);
const Fix_Factory SUFIX31(31, US, SAT);
const Fix_Factory SUFIX32(32, US, SAT);
const Fix_Factory SUFIX33(33, US, SAT);
const Fix_Factory SUFIX34(34, US, SAT);
const Fix_Factory SUFIX35(35, US, SAT);
const Fix_Factory SUFIX36(36, US, SAT);
const Fix_Factory SUFIX37(37, US, SAT);
const Fix_Factory SUFIX38(38, US, SAT);
const Fix_Factory SUFIX39(39, US, SAT);
const Fix_Factory SUFIX40(40, US, SAT);
const Fix_Factory SUFIX41(41, US, SAT);
const Fix_Factory SUFIX42(42, US, SAT);
const Fix_Factory SUFIX43(43, US, SAT);
const Fix_Factory SUFIX44(44, US, SAT);
const Fix_Factory SUFIX45(45, US, SAT);
const Fix_Factory SUFIX46(46, US, SAT);
const Fix_Factory SUFIX47(47, US, SAT);
const Fix_Factory SUFIX48(48, US, SAT);
const Fix_Factory SUFIX49(49, US, SAT);
const Fix_Factory SUFIX50(50, US, SAT);
const Fix_Factory SUFIX51(51, US, SAT);
const Fix_Factory SUFIX52(52, US, SAT);
const Fix_Factory SUFIX53(53, US, SAT);
const Fix_Factory SUFIX54(54, US, SAT);
const Fix_Factory SUFIX55(55, US, SAT);
const Fix_Factory SUFIX56(56, US, SAT);
const Fix_Factory SUFIX57(57, US, SAT);
const Fix_Factory SUFIX58(58, US, SAT);
const Fix_Factory SUFIX59(59, US, SAT);
const Fix_Factory SUFIX60(60, US, SAT);
const Fix_Factory SUFIX61(61, US, SAT);
const Fix_Factory SUFIX62(62, US, SAT);
const Fix_Factory SUFIX63(63, US, SAT);
const Fix_Factory SUFIX64(64, US, SAT);
//! \endcond

} // namespace itpp

#endif // #ifndef FIX_FACTORY_H
