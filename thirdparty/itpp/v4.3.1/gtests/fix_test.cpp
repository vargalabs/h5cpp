/*!
 * \file
 * \brief Fixed-point classes test program
 * \author Johan Bergman and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#include <itpp/itfixed.h>
#include <iomanip>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST (Fix, All)
{
  // This is a hack for improper rounding under MinGW
  static const double eps = 1e-4;

  // Testing declaration, initialization and conversion
  int shift(10);       // -64...+63 (0 is default)
  int wordlen(20);     // 1...64 (64 is default)
  e_mode emode(TC);    // TC or US (TC is default)
  o_mode omode(WRAP);  // WRAP or SAT (WRAP is default)
  q_mode qmode(TRN);   // RND or TRN (TRN is default)
  Stat *stat_ptr(0);   // 0 or Stat* value (0 is default)

  // For double and complex<double>
  double real_value(3.14159265358979323846);
  complex<double> complex_value(100.0 / 3.0, 200.0 / 3.0);

  // For Fix and CFix
  Fix the_fix(real_value, shift, wordlen, emode, omode, qmode, stat_ptr);
  ASSERT_NEAR(3.14062, double(the_fix), eps);
  CFix the_cfix(complex_value, 0.0, shift, wordlen, emode, omode, qmode,
  stat_ptr);
  complex<double> actual = complex<double>(the_cfix);
  ASSERT_NEAR(33.333, actual.real(), eps);
  ASSERT_NEAR(66.666, actual.imag(), eps);

  // For Fixed and CFixed
  Fixed<20, TC, WRAP, TRN> the_fixed(real_value, shift, stat_ptr);
  ASSERT_NEAR(3.14062, double(the_fixed), eps);
  CFixed<20, TC, WRAP, TRN> the_cfixed(complex_value, 0.0, shift, stat_ptr);
  actual = complex<double>(the_cfixed);
  ASSERT_NEAR(33.333, actual.real(), eps);
  ASSERT_NEAR(66.666, actual.imag(), eps);

  // For Fixed and CFixed declared using a typedef
  fixed20 the_fixed20(real_value, shift, stat_ptr);
  ASSERT_NEAR(3.14062, double(the_fixed20), eps);
  cfixed20 the_cfixed20(complex_value, 0.0, shift, stat_ptr);
  actual = complex<double>(the_cfixed20);
  ASSERT_NEAR(33.333, actual.real(), eps);
  ASSERT_NEAR(66.666, actual.imag(), eps);

  // For Fix and CFix declared using a factory
  Fix the_fix20(FIX20);
  the_fix20.set(real_value, shift);
  ASSERT_NEAR(3.14062, double(the_fix20), eps);
  CFix the_cfix20(FIX20);
  the_cfix20.set(complex_value, shift);
  actual = complex<double>(the_cfix20);
  ASSERT_NEAR(33.333, actual.real(), eps);
  ASSERT_NEAR(66.666, actual.imag(), eps);

  // Testing Array/Vec/Mat declarations and operations
  int vec_length(2);

  // For Vec<Fix> and Vec<CFix>
  fixvec the_fixvec(vec_length, FIX20);
  the_fixvec = Fix(real_value, shift);
  vec v_actual = to_vec(the_fixvec);
  vec v_expect = "3.14062 3.14062";
  int i;
  for (i = 0; i < v_actual.length(); ++i) {
    ASSERT_NEAR(v_actual[i], v_expect[i], eps);
  }

  cfixvec the_cfixvec(vec_length, FIX20);
  the_cfixvec = CFix(complex_value, 0.0, shift);
  cvec cv_actual = to_cvec(the_cfixvec);
  cvec cv_expect = "33.333+66.666i 33.333+66.666i";
  for (i = 0; i < cv_actual.length(); ++i) {
    ASSERT_NEAR(cv_actual[i].real(), cv_expect[i].real(), eps);
    ASSERT_NEAR(cv_actual[i].imag(), cv_expect[i].imag(), eps);
  }
  cv_actual = to_cvec(the_cfixvec + the_fixvec);
  cv_expect = "36.4736+66.666i 36.4736+66.666i";
  for (i = 0; i < cv_actual.length(); ++i) {
    ASSERT_NEAR(cv_actual[i].real(), cv_expect[i].real(), eps);
    ASSERT_NEAR(cv_actual[i].imag(), cv_expect[i].imag(), eps);
  }
  cv_actual = to_cvec(the_cfixvec - the_fixvec);
  cv_expect = "30.1924+66.666i 30.1924+66.666i";
  for (i = 0; i < cv_actual.length(); ++i) {
    ASSERT_NEAR(cv_actual[i].real(), cv_expect[i].real(), eps);
    ASSERT_NEAR(cv_actual[i].imag(), cv_expect[i].imag(), eps);
  }
  actual = complex<double>(the_cfixvec * the_fixvec);
  cv_expect = "209.373+418.746i";
  ASSERT_NEAR(209.373, actual.real(), eps);
  ASSERT_NEAR(418.746, actual.imag(), eps);
  cv_actual = to_cvec(the_cfixvec / the_fix);
  cv_expect = "10+21i 10+21i";
  for (i = 0; i < cv_actual.length(); ++i) {
    ASSERT_NEAR(cv_actual[i].real(), cv_expect[i].real(), eps);
    ASSERT_NEAR(cv_actual[i].imag(), cv_expect[i].imag(), eps);
  }

  // Testing functions

  // Function is_fix
  Array<Array<fixvec> > the_array2d_fixvec;
  ASSERT_TRUE(is_fix(the_array2d_fixvec));

  // Function set_fix
  vec original_float = "0:7";
  fixvec resulting_fix(FIX3);
  set_fix(resulting_fix, original_float, 0);
  v_actual = to_vec(resulting_fix);
  v_expect = "0 1 2 3 -4 -3 -2 -1";
  for (i = 0; i < v_actual.length(); ++i) {
    ASSERT_NEAR(v_actual[i], v_expect[i], eps);
  }
  vec resulting_float(FIX3);
  set_fix(resulting_float, original_float, 0);
  for (i = 0; i < resulting_float.length(); ++i) {
    ASSERT_NEAR(resulting_float[i], original_float[i], eps);
  }

  // Function lshift_fix
  Fix fix_to_be_lshifted(FIX16);
  fix_to_be_lshifted = 77;
  lshift_fix(fix_to_be_lshifted, 1);
  ASSERT_NEAR(77, double(fix_to_be_lshifted), eps);

  // Function rshift_fix
  Fix fix_to_be_rshifted(FIX16);
  fix_to_be_rshifted = Fix(3.14, 8);
  rshift_fix(fix_to_be_rshifted, 6, RND);
  ASSERT_NEAR(3.25, double(fix_to_be_rshifted), eps);
}
