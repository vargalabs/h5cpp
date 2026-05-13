/*!
 * \file
 * \brief Matrix class test program
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/itbase.h>
#include <iomanip>

using namespace itpp;
using namespace std;

template <typename T>
void common_operators(const Mat<T> &A, const Mat<T> &B, const Mat<T> &C,
                      const Vec<T> &u, const Vec<T> &v, T c)
{
  cout.setf(ios::fixed);
  cout.precision(4);

  cout << "A =\n" << A << endl;
  cout << "B =\n" << B << endl;
  cout << "C =\n" << C << endl;
  cout << "v = " << v << endl;
  cout << "u = " << u << endl;
  cout << "c = " << c << endl << endl;

  // indexing
  cout << "A(1,2) = " << A(1, 2) << endl;
  cout << "A(2,3) = " << A(2, 3) << endl;
  cout << "A(6) = " << A(6) << endl;
  cout << "A(0,2,1,3) =\n" << A(0, 2, 1, 3) << endl;
  cout << "A.get_row(1) = " << A.get_row(1) << endl;
  cout << "A.get_rows(1,2) =\n" << A.get_rows(1, 2) << endl;
  cout << "A.get_col(2) = " << A.get_col(2) << endl;
  cout << "A.get_cols(2,3) =\n" << A.get_cols(2, 3) << endl << endl;

  // setting, copying, swapping
  Mat<T> Mv(v);
  cout << "Mv(A) =\n" << Mv << endl;
  Mat<T> D(A);
  cout << "D(A) =\n" << D << endl;
  D.set_size(6, 5, true);
  cout << "D.set_size(6,5, true); D =\n" << D << endl;
  D.set_size(3, 2, true);
  cout << "D.set_size(3,2, true); D =\n" << D << endl;
  D.zeros();
  cout << "D.zeros(); D =\n" << D << endl;
  D.ones();
  cout << "D.ones(); D =\n" << D << endl;
  D = A;
  cout << "D = A; D =\n" << D << endl;
  D(2, 2) = c;
  cout << "D(2,2) = c; D =\n" << D << endl;
  D(9) = c;
  cout << "D(7) = c; D =\n" << D << endl;
  D.set(0, 1, c);
  cout << "D.set(0,1, c); D =\n" << D << endl;
  D.set_row(1, v);
  cout << "D.set_row(1, v); D =\n" << D << endl;
  D.set_col(2, u);
  cout << "D.set_col(2, u): D =\n" << D << endl;
  D.set_rows(0, B.get_rows(1, 2));
  cout << "D.set_rows(0, B.get_rows(1,2)); D =\n" << D << endl;
  D.set_cols(2, B.get_cols(0, 1));
  cout << "D.set_cols(2, B.get_cols(0,1)); D =\n" << D << endl;
  D.copy_row(1, 2);
  cout << "D.copy_row(1, 2); D =\n" << D << endl;
  D.copy_col(2, 3);
  cout << "D.copy_col(2, 3); D =\n" << D << endl;
  D.swap_rows(0, 2);
  cout << "D.swap_rows(0, 2); D =\n" << D << endl;
  D.swap_cols(0, 3);
  cout << "D.swap_cols(0, 3); D =\n" << D << endl;
  D.set_submatrix(1, 2, A(0, 1, 0, 1));
  cout << "D.set_submatrix(1,2, A(0,1,0,1); D =\n" << D << endl;
  D.set_submatrix(0, 0, A(0, 1, 0, 1));
  cout << "D.set_submatrix(0,0, A(0,1,0,1); D =\n" << D << endl;
  D.set_submatrix(1, 2, 2, 3, c);
  cout << "D.set_submatrix(1,2,2,3, c); D =\n" << D << endl << endl;

  // transposition
  cout << "A.T() =\n" << A.T() << endl;
  cout << "A.T().T() =\n" << A.T().T() << endl;
  cout << "A.H() =\n" << A.H() << endl << endl;

  // concatenation
  D = concat_horizontal(A, B);
  cout << "D = concat_horizontal(A,B); D =\n" << D << endl;
  D = concat_vertical(A, B);
  cout << "D = concat_vertical(A,B); D =\n" << D << endl << endl;

  // deleting rows, cols
  D.del_row(2);
  cout << "D.del_row(2); D =\n" << D << endl;
  D.del_rows(0, 2);
  cout << "D.del_rows(0,2); D =\n" << D << endl;
  D.del_col(3);
  cout << "D.del_col(3); D =\n" << D << endl;
  D.del_cols(0, 1);
  cout << "D.del_cols(0,1); D =\n" << D << endl << endl;

  // inserting, appending rows cols
  Mat<T> A2 = A;
  A2.ins_row(1, v);
  cout << "A.ins_row(1, v); A =\n" << A2 << endl;
  A2.ins_col(0, v);
  cout << "A.ins_col(0, v); A =\n" << A2 << endl;
  A2.append_col(A2.get_col(3));
  cout << "A.append_col(A2.get_col(3)); A =\n" << A2 << endl;
  A2.append_row(A2.get_row(0));
  cout << "A.append_row(A2.get_row(0)); A =\n" << A2 << endl << endl;

  // addition
  cout << "A+B =\n" << A + B << endl;
  cout << "A+c =\n" << A + c << endl;
  cout << "c+A =\n" << c + A << endl;
  A2 = A;
  A2 += B;
  cout << "A+=B; A =\n" << A2 << endl;
  A2 = A;
  A2 += c;
  cout << "A+=c; A =\n" << A2 << endl << endl;

  // subtraction
  cout << "A-B =\n" << A - B << endl;
  cout << "A-c =\n" << A - c << endl;
  cout << "c-A =\n" << c - A << endl;
  A2 = A;
  A2 -= B;
  cout << "A-=B; A =\n" << A2 << endl;
  A2 = A;
  A2 -= c;
  cout << "A-=c; A =\n" << A2 << endl;
  cout << "-A =\n" << -A << endl << endl;

  // multiplication
  cout << "A*C =\n" << A*C << endl;
  A2 = A;
  A2 *= C;
  cout << "A*=C; A =\n" << A2 << endl;
  cout << "A*c =\n" << A*c << endl;
  cout << "c*A =\n" << c*A << endl;
  A2 = A;
  A2 *= c;
  cout << "A*=c; A =\n" << A2 << endl;
  cout << "A*v = " << A*v << endl;
  cout << "elem_mult(A,B) =\n" << elem_mult(A, B) << endl;
  elem_mult_out(A, B, A2);
  cout << "elem_mult_out(A,B,out); out =\n" << A2 << endl;
  Mat<T> B2 = B;
  elem_mult_inplace(A, B2);
  cout << "elem_mult_inplace(A,B); B =\n" << B2 << endl;
  cout << "elem_mult_sum(A,B) = " << elem_mult_sum(A, B) << endl << endl;

  // division
  cout << "A/c =\n" << A / c << endl;
  A2 = A;
  A2 /= c;
  cout << "A/=c; A =\n" << A2 << endl;
  A2 = A;
  A2 /= B;
  cout << "A/=B; A =\n" << A2 << endl;
  cout << "elem_div(A,B) =\n" << elem_div(A, B) << endl;
  elem_div_out(A, B, A2);
  cout << "elem_div_out(A,B,out); out =\n" << A2 << endl;
  cout << "elem_div_sum(A,B) = " << elem_div_sum(A, B) << endl << endl;
}


int main()
{
  cout << "=============================" << endl
       << "   Testing Mat<bin> (bmat)" << endl
       << "=============================" << endl;
  bmat bM1 = randb(3, 4);
  bmat bM2 = randb(3, 4);
  bmat bM3 = randb(4, 3);
  bvec bv1 = randb(3);
  bvec bv2 = randb(4);
  bin bx = randb();
  common_operators(bM1, bM2, bM3, bv1, bv2, bx);

  cout << "=============================" << endl
       << "   Testing Mat<int> (imat)" << endl
       << "=============================" << endl;
  imat iM1 = randi(3, 4, 1, 9);
  imat iM2 = randi(3, 4, 1, 9);
  imat iM3 = randi(4, 3, 1, 9);
  ivec iv1 = randi(3, 1, 9);
  ivec iv2 = randi(4, 1, 9);
  int ix = randi(1, 9);
  common_operators(iM1, iM2, iM3, iv1, iv2, ix);

  cout << "===============================" << endl
       << "   Testing Mat<double> (mat)" << endl
       << "===============================" << endl;
  mat dM1 = randn(3, 4);
  mat dM2 = randn(3, 4);
  mat dM3 = randn(4, 3);
  vec dv1 = randn(3);
  vec dv2 = randn(4);
  double dx = randn();
  common_operators(dM1, dM2, dM3, dv1, dv2, dx);

  cout << "==========================================" << endl
       << "   Testing Mat<complex<double> > (cmat)" << endl
       << "==========================================" << endl;
  cmat cM1 = randn_c(3, 4);
  cmat cM2 = randn_c(3, 4);
  cmat cM3 = randn_c(4, 3);
  cvec cv1 = randn_c(3);
  cvec cv2 = randn_c(4);
  complex<double> cx = randn_c();
  common_operators(cM1, cM2, cM3, cv1, cv2, cx);


  cout << "========================================" << endl;
  cout << "   Testing initialisation with string" << endl;
  cout << "========================================" << endl;

  cout << "bmat M = \" 1 1 0; 0 1; 1 1 ,1,   1; ; 0 1\"" << endl;
  bmat bM = " 1 1 0; 0 1; 1 1 ,1,   1; ; 0 1";
  cout << "M =\n" << bM << endl;

  cout << "smat M = \"0xFF, -021 ,   100; 0,-0x01; 0xA, 10 012;  \"" << endl;
  smat sM = "0xFF, -021 ,   100; 0,-0x01; 0xA, 10 012;  ";
  cout << "M =\n" << sM << endl;

  cout << "imat M = \"0xFAC0, -021, 100000; 0,-0x01; 0xA, 10 012; \"" << endl;
  imat iM = "0xFAC0, -021, 100000; 0,-0x01; 0xA, 10 012; ";
  cout << "M =\n" << iM << endl;

  cout << "mat M = \".77 1e9 35; 0x7 3.5 -1000 \"";
  mat M = ".77 1.89e5 35; 7 3.5 -1000 ";
  cout << "M =\n" << M << endl;

  cout << "cmat M = \" 1.5+3i, (.33,1) ;  (333,-1) 2-0.2E-3i\"" << endl;
  cmat cM = " 1.5+3i, (.33,1) ;  (333,-1) 2-0.2E-3i";
  cout << "M =\n" << cM << endl << endl;

  return 0;
}
