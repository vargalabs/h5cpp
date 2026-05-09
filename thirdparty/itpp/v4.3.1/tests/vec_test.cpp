/*!
 * \file
 * \brief Vector class test program
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
void common_operators(const Vec<T> &a, const Vec<T> &b, T c)
{
  cout.setf(ios::fixed);
  cout.precision(4);

  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "c = " << c << endl;

  cout << "a+b = " << a + b << endl;
  cout << "a+c = " << a + c << endl;
  cout << "c+a = " << c + a << endl;
  Vec<T> a2 = a;
  a2 += b;
  cout << "a+=b; a = " << a2 << endl;
  a2 = a;
  a2 += c;
  cout << "a+=c; a = " << a2 << endl;

  cout << "a-b = " << a - b << endl;
  cout << "a-c = " << a - c << endl;
  cout << "c-a = " << c - a << endl;
  a2 = a;
  a2 -= b;
  cout << "a-=b; a = " << a2 << endl;
  a2 = a;
  a2 -= c;
  cout << "a-=c; a = " << a2 << endl;
  cout << "-a = " << -a << endl;

  cout << "a*b = " << a*b << endl;
  cout << "dot(a,b) = " << dot(a, b) << endl;
  cout << "outer_product(a,b) = " << outer_product(a, b) << endl;
  cout << "a*c = " << a*c << endl;
  cout << "c*a = " << c*a << endl;
  a2 = a;
  a2 *= c;
  cout << "a*=c; a = " << a2 << endl;
  cout << "elem_mult(a,b) = " << elem_mult(a, b) << endl;
  Vec<T> x1;
  elem_mult_out(a, b, x1);
  cout << "elem_mult_out(a,b,x); x = " << x1 << endl;
  Vec<T> b2 = b;
  elem_mult_inplace(a, b2);
  cout << "elem_mult_inplace(a,b); b = " << b2 << endl;
  cout << "elem_mult_sum(a,b) = " << elem_mult_sum(a, b) << endl;

  cout << "a/c = " << a / c << endl;
  cout << "c/a = " << c / a << endl;
  a2 = a;
  a2 /= c;
  cout << "a/=c; a = " << a2 << endl;
  a2 = a;
  a2 /= b;
  cout << "a/=b; a = " << a2 << endl;
  cout << "elem_div(a,b) = " << elem_div(a, b) << endl;
  Vec<T> x2;
  elem_div_out(a, b, x2);
  cout << "elem_div_out(a,b,x); x = " << x2 << endl;
  cout << "elem_div_sum(a,b) = " << elem_div_sum(a, b) << endl;

  cout << "concat(a,b) = " << concat(a, b) << endl;
  cout << "concat(a,c) = " << concat(a, c) << endl;
  cout << "concat(c,a) = " << concat(c, a) << endl;
  cout << "concat(a,b,a) = " << concat(a, b, a) << endl;
  cout << "concat(a,b,a,b) = " << concat(a, b, a, b) << endl;
  cout << "concat(a,b,a,b,a) = " << concat(a, b, a, b, a) << endl;

  cout << "a.T() = " << a.T() << endl;
  cout << "a.H() = " << a.H() << endl;

  cout << "a.size() = " << a.size() << endl;
  a2 = a;
  a2.set_size(a2.size() + 3, true);
  cout << "a.set_size(a.size()+3, true); a = " << a2 << endl;
  a2.set_size(a2.size() - 6, true);
  cout << "a.set_size(a.size()-6, true); a = " << a2 << endl;

  cout << "a(5) = " << a(5) << endl;
  cout << "a.get(5) = " << a.get(5) << endl;
  cout << "a(0,5) = " << a(0, 5) << endl;
  cout << "a.get(0,5) = " << a.get(0, 5) << endl;
  cout << "a(6,-1) = " << a(6, -1) << endl;
  ivec idx_list = "0 5 6 7";
  cout << "idx_list = " << idx_list << endl;
  cout << "a(idx_list) = " << a(idx_list) << endl;
  cout << "a.get(idx_list) = " << a.get(idx_list) << endl;
  bvec bin_list = "1 0 0 0 0 1 1 1 0 0";
  cout << "bin_list = " << bin_list << endl;
  cout << "a(bin_list) = " << a(bin_list) << endl;
  cout << "a.get(bin_list) = " << a.get(bin_list) << endl;
  cout << "a.right(3) = " << a.right(3) << endl;
  cout << "a.left(4) = " << a.left(4) << endl;
  cout << "a.mid(3,2) = " << a.mid(3, 2) << endl;
  a2 = a;
  cout << "a.split(0) = " << a2.split(0) << ";   a = " << a2 << endl;
  a2 = a;
  cout << "a.split(a.size()) = " << a2.split(a2.size()) << ";   a = " << a2
       << endl;
  a2 = a;
  cout << "a.split(4) = " << a2.split(4) << ";   a = " << a2 << endl;
  a2(5) = a(6);
  cout << "a(5) = a(6); a = " << a2 << endl;

  a2 = a;
  a2.shift_left(c, 2);
  cout << "a.shift_left(c,2) = " << a2 << endl;
  a2 = a;
  a2.shift_right(c);
  cout << "a.shift_right(c) = " << a2 << endl;
  a2 = a;
  a2.shift_left(b.mid(0, 2));
  cout << "a.shift_left(b.mid(0,2)) = " << a2 << endl;
  a2 = a;
  a2.shift_right(b.right(5));
  cout << "a.shift_right(b.right(5)) = " << a2 << endl;

  a2 = a;
  a2.set_subvector(0, b);
  cout << "a.set_subvector(0, b) = " << a2 << endl;
  a2 = a;
  a2.set_subvector(4, b(3, 5));
  cout << "a.set_subvector(4, b(3,5)) = " << a2 << endl;
  a2 = a;
  a2.replace_mid(4, b(3, 5));
  cout << "a.replace_mid(4, b(3,5)) = " << a2 << endl;
  a2 = a;
  a2.del(6);
  cout << "a.del(6) = " << a2 << endl;
  a2 = a;
  a2.del(3, 9);
  cout << "a.del(3,9) = " << a2 << endl;
  a2 = a;
  a2.ins(0, c);
  cout << "a.ins(0,c) = " << a2 << endl;
  a2 = a;
  a2.ins(2, c);
  cout << "a.ins(2,c) = " << a2 << endl;
  a2 = a;
  a2.ins(10, c);
  cout << "a.ins(10,c) = " << a2 << endl;
  a2 = a;
  a2.ins(3, b(0, 2));
  cout << "a.ins(3, b(0,2)) = " << a2 << endl;

  a2 = a;
  a2.zeros();
  cout << "a.zeros(); a = " << a2 << endl;
  a2 = a;
  a2.ones();
  cout << "a.ones(); a = " << a2 << endl;
  a2 = a;
  a2 = c;
  cout << "a = c; a = " << a2 << endl;
  a2 = a;
  a2 = b(0, 4);
  cout << "a = b(0,4); a = " << a2 << endl;
  a2 = a;
  a2 = b.T();
  cout << "a = b.T(); a = " << a2 << endl;
  a2 = a;
  a2 = b.T().T();
  cout << "a = b.T().T(); a = " << a2 << endl << endl;
}


template <typename T>
void logical_operators(const Vec<T> &a, const Vec<T> &b, T c)
{
  cout << "(a == c) = " << (a == c) << endl;
  cout << "(a != c) = " << (a != c) << endl;
  cout << "(a <= c) = " << (a <= c) << endl;
  cout << "(a >= c) = " << (a >= c) << endl;
  cout << "(a < c) = " << (a < c) << endl;
  cout << "(a > c) = " << (a > c) << endl;
  cout << "(a == b) = " << (a == b) << endl;
  cout << "(a != b) = " << (a != b) << endl;
  Vec<T> a2 = a;
  cout << "a2 = a; (a2 == a) = " << (a2 == a) << endl;
  cout << "a2 = a; (a2 != a) = " << (a2 != a) << endl << endl;
}


int main()
{
  cout << "=============================" << endl
       << "   Testing Vec<bin> (bvec)" << endl
       << "=============================" << endl;
  bvec bv1 = randb(10);
  bvec bv2 = randb(10);
  bin bx = randb();
  common_operators(bv1, bv2, bx);
  logical_operators(bv1, bv2, bx);

  cout << "=============================" << endl
       << "   Testing Vec<int> (ivec)" << endl
       << "=============================" << endl;
  ivec iv1 = randi(10, 1, 9);
  ivec iv2 = randi(10, 1, 9);
  int ix = randi(1, 9);
  common_operators(iv1, iv2, ix);
  logical_operators(iv1, iv2, ix);

  cout << "===============================" << endl
       << "   Testing Vec<double> (vec)" << endl
       << "===============================" << endl;
  vec dv1 = randu(10);
  vec dv2 = randu(10);
  double dx = randu();
  common_operators(dv1, dv2, dx);
  logical_operators(dv1, dv2, dx);

  cout << "===============================================" << endl
       << "   Testing Vec<std::complex<double> > (cvec)" << endl
       << "===============================================" << endl;
  cvec cv1 = randn_c(10);
  cvec cv2 = randn_c(10);
  complex<double> cx = randn_c();
  common_operators(cv1, cv2, cx);
  cout << "(a == c) = " << (cv1 == cx) << endl;
  cout << "(a != c) = " << (cv1 != cx) << endl << endl;


  // Test vectror initialisation with string
  vec v = "23.3 1232.7 0.111 1.525 0.333";
  cout << "Testing double vector initialisation with: "
    "\"23.3 1232.7 0.111 1.525 0.333\":" << endl << "v = " << v << endl;

  v = "-10.000 :.5:-4.5  1.33e+1, -.9, 1e0:1.5:1E+1";
  cout << "Testing double vector initialisation with: "
    "\"-10.000 :.5:-4.5  1.33e+1, -.9, 1e0:1.5:1E+1\":" << endl << "v = "
       << v << endl;

  ivec iv = "0xA :-0x1: -010";
  cout << "Testing int vector initialisation with: \"0xA :-0x1: -010\":"
       << endl << "iv = " << iv << endl;

  iv = "-5:3:9, 7, 1:10";
  cout << "Testing int vector initialisation with: \"-5:3:9, 7, 1:10\":"
       << endl << "iv = " << iv << endl;

  svec sv = "3 0xF -10, 0133 0177, 0x0 ";
  cout << "Testing short int vector initialisation with: \"3 0xF -10, 0133 0177, 0x0 \":"
       << endl << "sv = " << sv << endl;

  cvec cv = " (0.3, 0.4)  .2-.01i, 1e-3+0.25i";
  cout << "Testing complex vector initialisation with: \" (0.3, 0.4)  .2-.01i, 1e-3+0.25i\":"
       << endl << "cv = " << cv << endl;

  bvec bv = "1 1 0,1  1  ,  0 ,1  ";
  cout << "Testing bit vector initialisation with: \"1 1 0,1  1  ,  0 ,1  \":"
       << endl << "bv = " << bv << endl << endl;

  // Test of rem:
  v = "1.0 2.0 3.4 -4.5 6.7";
  double y = 0.76;
  cout << "v = " << v << endl;
  cout << "y = " << y << endl;
  cout << "rem(v, y) = " << rem(v, y) << endl;
  cout << "rem(10, v) = " << rem(10, v) << endl;
  mat M = "1.0 2.3; 4.5 -6.7";
  cout << "M = " << M << endl;
  cout << "rem(M, y) = " << rem(M, y) << endl;
  cout << "rem(10, M) = " << rem(10, M) << endl << endl;

  // Test of all and any:
  bvec b1 = "0 0 0 0 0 0 0 1 0 0";
  bvec b2 = "0 0 0 0 0 0 0 0 0 0";
  bvec b3 = "1 1 1 1 1 1 1 1 1 1 1 1 1";
  bvec b4 = "1 1 1 1 1 1 1 1 1 1 1 0 1";

  cout << "any(b1) = " << any(b1) << endl;
  cout << "any(b2) = " << any(b2) << endl;
  cout << "all(b3) = " << all(b3) << endl;
  cout << "all(b4) = " << all(b4) << endl;

  return 0;
}
