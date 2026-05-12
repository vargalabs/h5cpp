/*!
 * \file
 * \brief Definitions of a fixed-point data type Fixed
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

#ifndef FIXED_H
#define FIXED_H

#include <itpp/fixed/fix.h>
#include <itpp/itexports.h>


namespace itpp
{

//! \addtogroup fixed
//!@{

/*!
  \brief Templated fixed-point data type

  See the Detailed Description in the \ref fixed module.
*/
template < int w, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN >
class Fixed : public Fix
{
public:
  //! Default constructor
  Fixed(double x = 0.0, int s = 0, Stat *ptr = 0)
      : Fix(x, s, w, e, o, q, ptr) {}
  //! Constructor
  explicit Fixed(Stat *ptr)
      : Fix(0.0, 0, w, e, o, q, ptr) {}
  //! Constructor
  Fixed(const Fix &x, Stat *ptr = 0)
      : Fix(x, w, e, o, q, ptr) {}
  //! Destructor
  virtual ~Fixed() {}

  //! Assignment from Fix
  Fixed& operator=(const Fix &x) {
    shift = x.shift;
    re = apply_o_mode(x.re);
    return *this;
  }
  //! Assignment from int
  Fixed& operator=(int x) {
    shift = 0;
    re = apply_o_mode(x);
    return *this;
  }
protected:
};

//!@}

//! Typedefs for Fixed (fixed1, fixed2, ..., fixed64)
typedef Fixed<1, TC, WRAP> fixed1;
//! \cond
typedef Fixed<2, TC, WRAP> fixed2;
typedef Fixed<3, TC, WRAP> fixed3;
typedef Fixed<4, TC, WRAP> fixed4;
typedef Fixed<5, TC, WRAP> fixed5;
typedef Fixed<6, TC, WRAP> fixed6;
typedef Fixed<7, TC, WRAP> fixed7;
typedef Fixed<8, TC, WRAP> fixed8;
typedef Fixed<9, TC, WRAP> fixed9;
typedef Fixed<10, TC, WRAP> fixed10;
typedef Fixed<11, TC, WRAP> fixed11;
typedef Fixed<12, TC, WRAP> fixed12;
typedef Fixed<13, TC, WRAP> fixed13;
typedef Fixed<14, TC, WRAP> fixed14;
typedef Fixed<15, TC, WRAP> fixed15;
typedef Fixed<16, TC, WRAP> fixed16;
typedef Fixed<17, TC, WRAP> fixed17;
typedef Fixed<18, TC, WRAP> fixed18;
typedef Fixed<19, TC, WRAP> fixed19;
typedef Fixed<20, TC, WRAP> fixed20;
typedef Fixed<21, TC, WRAP> fixed21;
typedef Fixed<22, TC, WRAP> fixed22;
typedef Fixed<23, TC, WRAP> fixed23;
typedef Fixed<24, TC, WRAP> fixed24;
typedef Fixed<25, TC, WRAP> fixed25;
typedef Fixed<26, TC, WRAP> fixed26;
typedef Fixed<27, TC, WRAP> fixed27;
typedef Fixed<28, TC, WRAP> fixed28;
typedef Fixed<29, TC, WRAP> fixed29;
typedef Fixed<30, TC, WRAP> fixed30;
typedef Fixed<31, TC, WRAP> fixed31;
typedef Fixed<32, TC, WRAP> fixed32;
typedef Fixed<33, TC, WRAP> fixed33;
typedef Fixed<34, TC, WRAP> fixed34;
typedef Fixed<35, TC, WRAP> fixed35;
typedef Fixed<36, TC, WRAP> fixed36;
typedef Fixed<37, TC, WRAP> fixed37;
typedef Fixed<38, TC, WRAP> fixed38;
typedef Fixed<39, TC, WRAP> fixed39;
typedef Fixed<40, TC, WRAP> fixed40;
typedef Fixed<41, TC, WRAP> fixed41;
typedef Fixed<42, TC, WRAP> fixed42;
typedef Fixed<43, TC, WRAP> fixed43;
typedef Fixed<44, TC, WRAP> fixed44;
typedef Fixed<45, TC, WRAP> fixed45;
typedef Fixed<46, TC, WRAP> fixed46;
typedef Fixed<47, TC, WRAP> fixed47;
typedef Fixed<48, TC, WRAP> fixed48;
typedef Fixed<49, TC, WRAP> fixed49;
typedef Fixed<50, TC, WRAP> fixed50;
typedef Fixed<51, TC, WRAP> fixed51;
typedef Fixed<52, TC, WRAP> fixed52;
typedef Fixed<53, TC, WRAP> fixed53;
typedef Fixed<54, TC, WRAP> fixed54;
typedef Fixed<55, TC, WRAP> fixed55;
typedef Fixed<56, TC, WRAP> fixed56;
typedef Fixed<57, TC, WRAP> fixed57;
typedef Fixed<58, TC, WRAP> fixed58;
typedef Fixed<59, TC, WRAP> fixed59;
typedef Fixed<60, TC, WRAP> fixed60;
typedef Fixed<61, TC, WRAP> fixed61;
typedef Fixed<62, TC, WRAP> fixed62;
typedef Fixed<63, TC, WRAP> fixed63;
typedef Fixed<64, TC, WRAP> fixed64;
//! \endcond

//! Typedefs for unsigned Fixed (ufixed1, ufixed2, ..., ufixed64)
typedef Fixed<1, US, WRAP> ufixed1;
//! \cond
typedef Fixed<2, US, WRAP> ufixed2;
typedef Fixed<3, US, WRAP> ufixed3;
typedef Fixed<4, US, WRAP> ufixed4;
typedef Fixed<5, US, WRAP> ufixed5;
typedef Fixed<6, US, WRAP> ufixed6;
typedef Fixed<7, US, WRAP> ufixed7;
typedef Fixed<8, US, WRAP> ufixed8;
typedef Fixed<9, US, WRAP> ufixed9;
typedef Fixed<10, US, WRAP> ufixed10;
typedef Fixed<11, US, WRAP> ufixed11;
typedef Fixed<12, US, WRAP> ufixed12;
typedef Fixed<13, US, WRAP> ufixed13;
typedef Fixed<14, US, WRAP> ufixed14;
typedef Fixed<15, US, WRAP> ufixed15;
typedef Fixed<16, US, WRAP> ufixed16;
typedef Fixed<17, US, WRAP> ufixed17;
typedef Fixed<18, US, WRAP> ufixed18;
typedef Fixed<19, US, WRAP> ufixed19;
typedef Fixed<20, US, WRAP> ufixed20;
typedef Fixed<21, US, WRAP> ufixed21;
typedef Fixed<22, US, WRAP> ufixed22;
typedef Fixed<23, US, WRAP> ufixed23;
typedef Fixed<24, US, WRAP> ufixed24;
typedef Fixed<25, US, WRAP> ufixed25;
typedef Fixed<26, US, WRAP> ufixed26;
typedef Fixed<27, US, WRAP> ufixed27;
typedef Fixed<28, US, WRAP> ufixed28;
typedef Fixed<29, US, WRAP> ufixed29;
typedef Fixed<30, US, WRAP> ufixed30;
typedef Fixed<31, US, WRAP> ufixed31;
typedef Fixed<32, US, WRAP> ufixed32;
typedef Fixed<33, US, WRAP> ufixed33;
typedef Fixed<34, US, WRAP> ufixed34;
typedef Fixed<35, US, WRAP> ufixed35;
typedef Fixed<36, US, WRAP> ufixed36;
typedef Fixed<37, US, WRAP> ufixed37;
typedef Fixed<38, US, WRAP> ufixed38;
typedef Fixed<39, US, WRAP> ufixed39;
typedef Fixed<40, US, WRAP> ufixed40;
typedef Fixed<41, US, WRAP> ufixed41;
typedef Fixed<42, US, WRAP> ufixed42;
typedef Fixed<43, US, WRAP> ufixed43;
typedef Fixed<44, US, WRAP> ufixed44;
typedef Fixed<45, US, WRAP> ufixed45;
typedef Fixed<46, US, WRAP> ufixed46;
typedef Fixed<47, US, WRAP> ufixed47;
typedef Fixed<48, US, WRAP> ufixed48;
typedef Fixed<49, US, WRAP> ufixed49;
typedef Fixed<50, US, WRAP> ufixed50;
typedef Fixed<51, US, WRAP> ufixed51;
typedef Fixed<52, US, WRAP> ufixed52;
typedef Fixed<53, US, WRAP> ufixed53;
typedef Fixed<54, US, WRAP> ufixed54;
typedef Fixed<55, US, WRAP> ufixed55;
typedef Fixed<56, US, WRAP> ufixed56;
typedef Fixed<57, US, WRAP> ufixed57;
typedef Fixed<58, US, WRAP> ufixed58;
typedef Fixed<59, US, WRAP> ufixed59;
typedef Fixed<60, US, WRAP> ufixed60;
typedef Fixed<61, US, WRAP> ufixed61;
typedef Fixed<62, US, WRAP> ufixed62;
typedef Fixed<63, US, WRAP> ufixed63;
typedef Fixed<64, US, WRAP> ufixed64;
//! \endcond

//! Typedefs for saturated Fixed (sfixed1, sfixed2, ..., sfixed64)
typedef Fixed<1, TC, SAT> sfixed1;
//! \cond
typedef Fixed<2, TC, SAT> sfixed2;
typedef Fixed<3, TC, SAT> sfixed3;
typedef Fixed<4, TC, SAT> sfixed4;
typedef Fixed<5, TC, SAT> sfixed5;
typedef Fixed<6, TC, SAT> sfixed6;
typedef Fixed<7, TC, SAT> sfixed7;
typedef Fixed<8, TC, SAT> sfixed8;
typedef Fixed<9, TC, SAT> sfixed9;
typedef Fixed<10, TC, SAT> sfixed10;
typedef Fixed<11, TC, SAT> sfixed11;
typedef Fixed<12, TC, SAT> sfixed12;
typedef Fixed<13, TC, SAT> sfixed13;
typedef Fixed<14, TC, SAT> sfixed14;
typedef Fixed<15, TC, SAT> sfixed15;
typedef Fixed<16, TC, SAT> sfixed16;
typedef Fixed<17, TC, SAT> sfixed17;
typedef Fixed<18, TC, SAT> sfixed18;
typedef Fixed<19, TC, SAT> sfixed19;
typedef Fixed<20, TC, SAT> sfixed20;
typedef Fixed<21, TC, SAT> sfixed21;
typedef Fixed<22, TC, SAT> sfixed22;
typedef Fixed<23, TC, SAT> sfixed23;
typedef Fixed<24, TC, SAT> sfixed24;
typedef Fixed<25, TC, SAT> sfixed25;
typedef Fixed<26, TC, SAT> sfixed26;
typedef Fixed<27, TC, SAT> sfixed27;
typedef Fixed<28, TC, SAT> sfixed28;
typedef Fixed<29, TC, SAT> sfixed29;
typedef Fixed<30, TC, SAT> sfixed30;
typedef Fixed<31, TC, SAT> sfixed31;
typedef Fixed<32, TC, SAT> sfixed32;
typedef Fixed<33, TC, SAT> sfixed33;
typedef Fixed<34, TC, SAT> sfixed34;
typedef Fixed<35, TC, SAT> sfixed35;
typedef Fixed<36, TC, SAT> sfixed36;
typedef Fixed<37, TC, SAT> sfixed37;
typedef Fixed<38, TC, SAT> sfixed38;
typedef Fixed<39, TC, SAT> sfixed39;
typedef Fixed<40, TC, SAT> sfixed40;
typedef Fixed<41, TC, SAT> sfixed41;
typedef Fixed<42, TC, SAT> sfixed42;
typedef Fixed<43, TC, SAT> sfixed43;
typedef Fixed<44, TC, SAT> sfixed44;
typedef Fixed<45, TC, SAT> sfixed45;
typedef Fixed<46, TC, SAT> sfixed46;
typedef Fixed<47, TC, SAT> sfixed47;
typedef Fixed<48, TC, SAT> sfixed48;
typedef Fixed<49, TC, SAT> sfixed49;
typedef Fixed<50, TC, SAT> sfixed50;
typedef Fixed<51, TC, SAT> sfixed51;
typedef Fixed<52, TC, SAT> sfixed52;
typedef Fixed<53, TC, SAT> sfixed53;
typedef Fixed<54, TC, SAT> sfixed54;
typedef Fixed<55, TC, SAT> sfixed55;
typedef Fixed<56, TC, SAT> sfixed56;
typedef Fixed<57, TC, SAT> sfixed57;
typedef Fixed<58, TC, SAT> sfixed58;
typedef Fixed<59, TC, SAT> sfixed59;
typedef Fixed<60, TC, SAT> sfixed60;
typedef Fixed<61, TC, SAT> sfixed61;
typedef Fixed<62, TC, SAT> sfixed62;
typedef Fixed<63, TC, SAT> sfixed63;
typedef Fixed<64, TC, SAT> sfixed64;
//! \endcond

//! Typedefs for saturated unsigned Fixed (sufixed1, sufixed2, ..., sufixed64)
typedef Fixed<1, US, SAT> sufixed1;
//! \cond
typedef Fixed<2, US, SAT> sufixed2;
typedef Fixed<3, US, SAT> sufixed3;
typedef Fixed<4, US, SAT> sufixed4;
typedef Fixed<5, US, SAT> sufixed5;
typedef Fixed<6, US, SAT> sufixed6;
typedef Fixed<7, US, SAT> sufixed7;
typedef Fixed<8, US, SAT> sufixed8;
typedef Fixed<9, US, SAT> sufixed9;
typedef Fixed<10, US, SAT> sufixed10;
typedef Fixed<11, US, SAT> sufixed11;
typedef Fixed<12, US, SAT> sufixed12;
typedef Fixed<13, US, SAT> sufixed13;
typedef Fixed<14, US, SAT> sufixed14;
typedef Fixed<15, US, SAT> sufixed15;
typedef Fixed<16, US, SAT> sufixed16;
typedef Fixed<17, US, SAT> sufixed17;
typedef Fixed<18, US, SAT> sufixed18;
typedef Fixed<19, US, SAT> sufixed19;
typedef Fixed<20, US, SAT> sufixed20;
typedef Fixed<21, US, SAT> sufixed21;
typedef Fixed<22, US, SAT> sufixed22;
typedef Fixed<23, US, SAT> sufixed23;
typedef Fixed<24, US, SAT> sufixed24;
typedef Fixed<25, US, SAT> sufixed25;
typedef Fixed<26, US, SAT> sufixed26;
typedef Fixed<27, US, SAT> sufixed27;
typedef Fixed<28, US, SAT> sufixed28;
typedef Fixed<29, US, SAT> sufixed29;
typedef Fixed<30, US, SAT> sufixed30;
typedef Fixed<31, US, SAT> sufixed31;
typedef Fixed<32, US, SAT> sufixed32;
typedef Fixed<33, US, SAT> sufixed33;
typedef Fixed<34, US, SAT> sufixed34;
typedef Fixed<35, US, SAT> sufixed35;
typedef Fixed<36, US, SAT> sufixed36;
typedef Fixed<37, US, SAT> sufixed37;
typedef Fixed<38, US, SAT> sufixed38;
typedef Fixed<39, US, SAT> sufixed39;
typedef Fixed<40, US, SAT> sufixed40;
typedef Fixed<41, US, SAT> sufixed41;
typedef Fixed<42, US, SAT> sufixed42;
typedef Fixed<43, US, SAT> sufixed43;
typedef Fixed<44, US, SAT> sufixed44;
typedef Fixed<45, US, SAT> sufixed45;
typedef Fixed<46, US, SAT> sufixed46;
typedef Fixed<47, US, SAT> sufixed47;
typedef Fixed<48, US, SAT> sufixed48;
typedef Fixed<49, US, SAT> sufixed49;
typedef Fixed<50, US, SAT> sufixed50;
typedef Fixed<51, US, SAT> sufixed51;
typedef Fixed<52, US, SAT> sufixed52;
typedef Fixed<53, US, SAT> sufixed53;
typedef Fixed<54, US, SAT> sufixed54;
typedef Fixed<55, US, SAT> sufixed55;
typedef Fixed<56, US, SAT> sufixed56;
typedef Fixed<57, US, SAT> sufixed57;
typedef Fixed<58, US, SAT> sufixed58;
typedef Fixed<59, US, SAT> sufixed59;
typedef Fixed<60, US, SAT> sufixed60;
typedef Fixed<61, US, SAT> sufixed61;
typedef Fixed<62, US, SAT> sufixed62;
typedef Fixed<63, US, SAT> sufixed63;
typedef Fixed<64, US, SAT> sufixed64;

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Fixed<64, TC, WRAP>;

//! \endcond

} // namespace itpp

#endif // #ifndef FIXED_H
