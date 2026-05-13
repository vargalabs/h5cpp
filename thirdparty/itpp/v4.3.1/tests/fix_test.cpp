/*!
 * \file
 * \brief Fixed-point classes test program
 * \author Johan Bergman and Adam Piatyszek
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

#include <itpp/itfixed.h>
#include <iomanip>

using namespace itpp;
using namespace std;

int main()
{
  // This is a hack for improper rounding under MinGW
  cout.precision(8);

  cout << "Testing declaration, initialization and conversion" << endl;
  cout << "==================================================" << endl;
  int shift(10);       // -64...+63 (0 is default)
  int wordlen(20);     // 1...64 (64 is default)
  e_mode emode(TC);    // TC or US (TC is default)
  o_mode omode(WRAP);  // WRAP or SAT (WRAP is default)
  q_mode qmode(TRN);   // RND or TRN (TRN is default)
  Stat *stat_ptr(0);   // 0 or Stat* value (0 is default)

  cout << "For double and complex<double>:" << endl;
  double real_value(3.14159265358979323846);
  cout << "  real_value = " << real_value << endl;
  complex<double> complex_value(100.0 / 3.0, 200.0 / 3.0);
  cout << "  complex_value = " << complex_value << endl;

  cout << "For Fix and CFix:" << endl;
  Fix the_fix(real_value, shift, wordlen, emode, omode, qmode, stat_ptr);
  cout << "  the_fix = " << double(the_fix) << endl;
  CFix the_cfix(complex_value, 0.0, shift, wordlen, emode, omode, qmode,
                stat_ptr);
  cout << "  the_cfix = " << complex<double>(the_cfix) << endl;

  cout << "For Fixed and CFixed:" << endl;
  Fixed<20, TC, WRAP, TRN> the_fixed(real_value, shift, stat_ptr);
  cout << "  the_fixed = " << double(the_fixed) << endl;
  CFixed<20, TC, WRAP, TRN> the_cfixed(complex_value, 0.0, shift, stat_ptr);
  cout << "  the_cfixed = " << complex<double>(the_cfixed) << endl;

  cout << "For Fixed and CFixed declared using a typedef:" << endl;
  fixed20 the_fixed20(real_value, shift, stat_ptr);
  cout << "  the_fixed20 = " << double(the_fixed20) << endl;
  cfixed20 the_cfixed20(complex_value, 0.0, shift, stat_ptr);
  cout << "  the_cfixed20 = " << complex<double>(the_cfixed20) << endl;

  cout << "For Fix and CFix declared using a factory:" << endl;
  Fix the_fix20(FIX20);
  the_fix20.set(real_value, shift);
  cout << "  the_fix20 = " << double(the_fix20) << endl;
  CFix the_cfix20(FIX20);
  the_cfix20.set(complex_value, shift);
  cout << "  the_cfix20 = " << complex<double>(the_cfix20) << endl << endl;

  cout << "Testing Array/Vec/Mat declarations and operations" << endl;
  cout << "=================================================" << endl;
  int vec_length(2);

  cout << "For Vec<Fix> and Vec<CFix>:" << endl;
  fixvec the_fixvec(vec_length, FIX20);
  the_fixvec = Fix(real_value, shift);
  cout << "  the_fixvec = " << to_vec(the_fixvec) << endl;
  cfixvec the_cfixvec(vec_length, FIX20);
  the_cfixvec = CFix(complex_value, 0.0, shift);
  cout << "  the_cfixvec = " << to_cvec(the_cfixvec) << endl;
  cout << "  the_cfixvec + the_fixvec = " << to_cvec(the_cfixvec + the_fixvec)
       << endl;
  cout << "  the_cfixvec - the_fixvec = " << to_cvec(the_cfixvec - the_fixvec)
       << endl;
  cout << "  the_cfixvec * the_fixvec = "
       << complex<double>(the_cfixvec * the_fixvec) << endl;
  cout << "  the_cfixvec / the_fix = " << to_cvec(the_cfixvec / the_fix)
       << endl << endl;

  cout << "Testing functions" << endl;
  cout << "=================" << endl;

  cout << "Function is_fix:" << endl;
  Array<Array<fixvec> > the_array2d_fixvec;
  cout << "  is_fix(the_array2d_fixvec) = " << is_fix(the_array2d_fixvec)
       << endl;

  cout << "Function set_fix:" << endl;
  vec original_float = "0:7";
  fixvec resulting_fix(FIX3);
  set_fix(resulting_fix, original_float, 0);
  cout << "  original_float = " << original_float << " => resulting_fix = "
       << resulting_fix << endl;
  vec resulting_float(FIX3);
  set_fix(resulting_float, original_float, 0);
  cout << "  original_float = " << original_float << " => resulting_float = "
       << resulting_float << endl;

  cout << "Function lshift_fix:" << endl;
  Fix fix_to_be_lshifted(FIX16);
  fix_to_be_lshifted = 77;
  cout << "  before lshift: " << fix_to_be_lshifted << " , rep: "
       << fix_to_be_lshifted.get_re() << endl;
  lshift_fix(fix_to_be_lshifted, 1);
  cout << "  after lshift: " << fix_to_be_lshifted << " , rep: "
       << fix_to_be_lshifted.get_re() << endl;

  cout << "Function rshift_fix:" << endl;
  Fix fix_to_be_rshifted(FIX16);
  fix_to_be_rshifted = Fix(3.14, 8);
  cout << "  before rshift: " << fix_to_be_rshifted << " , rep: "
       << fix_to_be_rshifted.get_re() << endl;
  rshift_fix(fix_to_be_rshifted, 6, RND);
  cout << "  after rshift: " << fix_to_be_rshifted << " , rep: "
       << fix_to_be_rshifted.get_re() << endl;

  return 0;
}
