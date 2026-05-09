/*!
 * \file
 * \brief Definitions of a complex fixed-point data type CFixed
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

#ifndef CFIXED_H
#define CFIXED_H

#include <itpp/fixed/cfix.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
 * \addtogroup fixed
 * @{
 */

/*!
  \brief Templated complex fixed-point data type

  See the Detailed Description in the \ref fixed module.
*/
template < int w, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN >
class CFixed : public CFix
{
public:
  //! Default constructor
  CFixed(double r = 0.0, double i = 0.0, int s = 0, Stat *ptr = 0)
      : CFix(r, i, s, w, e, o, q, ptr) {}
  //! Constructor
  CFixed(std::complex<double> x, double, int s = 0, Stat *ptr = 0)
      : CFix(x, 0.0, s, w, e, o, q, ptr) {}
  //! Constructor
  explicit CFixed(Stat *ptr)
      : CFix(0.0, 0.0, 0, w, e, o, q, ptr) {}
  //! Constructor
  CFixed(const Fix &r, const Fix &i = 0.0, Stat *ptr = 0)
      : CFix(r, i, w, e, o, q, ptr) {}
  //! Constructor
  CFixed(const CFix &x, double, Stat *ptr = 0)
      : CFix(x, 0.0, w, e, o, q, ptr) {}
  //! Destructor
  virtual ~CFixed() {}

  //! Assignment from CFix
  CFixed& operator=(const CFix &x) {
    shift = x.shift;
    re = apply_o_mode(x.re);
    im = apply_o_mode(x.im);
    return *this;
  }
  //! Assignment from Fix
  CFixed& operator=(const Fix &x) {
    shift = x.shift;
    re = apply_o_mode(x.re);
    im = 0;
    return *this;
  }
  //! Assignment from complex<double>. Fractional part is truncated
  CFixed& operator=(const std::complex<double> &x) {
    shift = 0;
    re = apply_o_mode(fixrep(real(x)));
    im = apply_o_mode(fixrep(imag(x)));
    return *this;
  }
  //! Assignment from int
  CFixed& operator=(int x) {
    shift = 0;
    re = apply_o_mode(x);
    im = 0;
    return *this;
  }
protected:
};

/*! @} */

//! Typedefs for CFixed (cfixed1, cfixed2, ..., cfixed64)
typedef CFixed<1, TC, WRAP> cfixed1;
//! \cond
typedef CFixed<2, TC, WRAP> cfixed2;
typedef CFixed<3, TC, WRAP> cfixed3;
typedef CFixed<4, TC, WRAP> cfixed4;
typedef CFixed<5, TC, WRAP> cfixed5;
typedef CFixed<6, TC, WRAP> cfixed6;
typedef CFixed<7, TC, WRAP> cfixed7;
typedef CFixed<8, TC, WRAP> cfixed8;
typedef CFixed<9, TC, WRAP> cfixed9;
typedef CFixed<10, TC, WRAP> cfixed10;
typedef CFixed<11, TC, WRAP> cfixed11;
typedef CFixed<12, TC, WRAP> cfixed12;
typedef CFixed<13, TC, WRAP> cfixed13;
typedef CFixed<14, TC, WRAP> cfixed14;
typedef CFixed<15, TC, WRAP> cfixed15;
typedef CFixed<16, TC, WRAP> cfixed16;
typedef CFixed<17, TC, WRAP> cfixed17;
typedef CFixed<18, TC, WRAP> cfixed18;
typedef CFixed<19, TC, WRAP> cfixed19;
typedef CFixed<20, TC, WRAP> cfixed20;
typedef CFixed<21, TC, WRAP> cfixed21;
typedef CFixed<22, TC, WRAP> cfixed22;
typedef CFixed<23, TC, WRAP> cfixed23;
typedef CFixed<24, TC, WRAP> cfixed24;
typedef CFixed<25, TC, WRAP> cfixed25;
typedef CFixed<26, TC, WRAP> cfixed26;
typedef CFixed<27, TC, WRAP> cfixed27;
typedef CFixed<28, TC, WRAP> cfixed28;
typedef CFixed<29, TC, WRAP> cfixed29;
typedef CFixed<30, TC, WRAP> cfixed30;
typedef CFixed<31, TC, WRAP> cfixed31;
typedef CFixed<32, TC, WRAP> cfixed32;
typedef CFixed<33, TC, WRAP> cfixed33;
typedef CFixed<34, TC, WRAP> cfixed34;
typedef CFixed<35, TC, WRAP> cfixed35;
typedef CFixed<36, TC, WRAP> cfixed36;
typedef CFixed<37, TC, WRAP> cfixed37;
typedef CFixed<38, TC, WRAP> cfixed38;
typedef CFixed<39, TC, WRAP> cfixed39;
typedef CFixed<40, TC, WRAP> cfixed40;
typedef CFixed<41, TC, WRAP> cfixed41;
typedef CFixed<42, TC, WRAP> cfixed42;
typedef CFixed<43, TC, WRAP> cfixed43;
typedef CFixed<44, TC, WRAP> cfixed44;
typedef CFixed<45, TC, WRAP> cfixed45;
typedef CFixed<46, TC, WRAP> cfixed46;
typedef CFixed<47, TC, WRAP> cfixed47;
typedef CFixed<48, TC, WRAP> cfixed48;
typedef CFixed<49, TC, WRAP> cfixed49;
typedef CFixed<50, TC, WRAP> cfixed50;
typedef CFixed<51, TC, WRAP> cfixed51;
typedef CFixed<52, TC, WRAP> cfixed52;
typedef CFixed<53, TC, WRAP> cfixed53;
typedef CFixed<54, TC, WRAP> cfixed54;
typedef CFixed<55, TC, WRAP> cfixed55;
typedef CFixed<56, TC, WRAP> cfixed56;
typedef CFixed<57, TC, WRAP> cfixed57;
typedef CFixed<58, TC, WRAP> cfixed58;
typedef CFixed<59, TC, WRAP> cfixed59;
typedef CFixed<60, TC, WRAP> cfixed60;
typedef CFixed<61, TC, WRAP> cfixed61;
typedef CFixed<62, TC, WRAP> cfixed62;
typedef CFixed<63, TC, WRAP> cfixed63;
typedef CFixed<64, TC, WRAP> cfixed64;
//! \endcond

//! Typedefs for saturated CFixed (scfixed1, scfixed2, ..., scfixed64)
typedef CFixed<1, TC, WRAP> cfixed1;
//! \cond
typedef CFixed<1, TC, SAT> scfixed1;
typedef CFixed<2, TC, SAT> scfixed2;
typedef CFixed<3, TC, SAT> scfixed3;
typedef CFixed<4, TC, SAT> scfixed4;
typedef CFixed<5, TC, SAT> scfixed5;
typedef CFixed<6, TC, SAT> scfixed6;
typedef CFixed<7, TC, SAT> scfixed7;
typedef CFixed<8, TC, SAT> scfixed8;
typedef CFixed<9, TC, SAT> scfixed9;
typedef CFixed<10, TC, SAT> scfixed10;
typedef CFixed<11, TC, SAT> scfixed11;
typedef CFixed<12, TC, SAT> scfixed12;
typedef CFixed<13, TC, SAT> scfixed13;
typedef CFixed<14, TC, SAT> scfixed14;
typedef CFixed<15, TC, SAT> scfixed15;
typedef CFixed<16, TC, SAT> scfixed16;
typedef CFixed<17, TC, SAT> scfixed17;
typedef CFixed<18, TC, SAT> scfixed18;
typedef CFixed<19, TC, SAT> scfixed19;
typedef CFixed<20, TC, SAT> scfixed20;
typedef CFixed<21, TC, SAT> scfixed21;
typedef CFixed<22, TC, SAT> scfixed22;
typedef CFixed<23, TC, SAT> scfixed23;
typedef CFixed<24, TC, SAT> scfixed24;
typedef CFixed<25, TC, SAT> scfixed25;
typedef CFixed<26, TC, SAT> scfixed26;
typedef CFixed<27, TC, SAT> scfixed27;
typedef CFixed<28, TC, SAT> scfixed28;
typedef CFixed<29, TC, SAT> scfixed29;
typedef CFixed<30, TC, SAT> scfixed30;
typedef CFixed<31, TC, SAT> scfixed31;
typedef CFixed<32, TC, SAT> scfixed32;
typedef CFixed<33, TC, SAT> scfixed33;
typedef CFixed<34, TC, SAT> scfixed34;
typedef CFixed<35, TC, SAT> scfixed35;
typedef CFixed<36, TC, SAT> scfixed36;
typedef CFixed<37, TC, SAT> scfixed37;
typedef CFixed<38, TC, SAT> scfixed38;
typedef CFixed<39, TC, SAT> scfixed39;
typedef CFixed<40, TC, SAT> scfixed40;
typedef CFixed<41, TC, SAT> scfixed41;
typedef CFixed<42, TC, SAT> scfixed42;
typedef CFixed<43, TC, SAT> scfixed43;
typedef CFixed<44, TC, SAT> scfixed44;
typedef CFixed<45, TC, SAT> scfixed45;
typedef CFixed<46, TC, SAT> scfixed46;
typedef CFixed<47, TC, SAT> scfixed47;
typedef CFixed<48, TC, SAT> scfixed48;
typedef CFixed<49, TC, SAT> scfixed49;
typedef CFixed<50, TC, SAT> scfixed50;
typedef CFixed<51, TC, SAT> scfixed51;
typedef CFixed<52, TC, SAT> scfixed52;
typedef CFixed<53, TC, SAT> scfixed53;
typedef CFixed<54, TC, SAT> scfixed54;
typedef CFixed<55, TC, SAT> scfixed55;
typedef CFixed<56, TC, SAT> scfixed56;
typedef CFixed<57, TC, SAT> scfixed57;
typedef CFixed<58, TC, SAT> scfixed58;
typedef CFixed<59, TC, SAT> scfixed59;
typedef CFixed<60, TC, SAT> scfixed60;
typedef CFixed<61, TC, SAT> scfixed61;
typedef CFixed<62, TC, SAT> scfixed62;
typedef CFixed<63, TC, SAT> scfixed63;
typedef CFixed<64, TC, SAT> scfixed64;

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT CFixed<64, TC, WRAP>;

//! \endcond

} // namespace itpp

#endif // #ifndef CFIXED_H
