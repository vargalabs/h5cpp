/*!
 * \file
 * \brief Matrix class test program
 * \author Tony Ottosson, Adam Piatyszek and Bogdan Cristea
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

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
const double tol = 1e-5;

static
void common_operators(const bmat &A, const bmat &B, const bmat &C,
                      const bvec &u, const bvec &v, bin c)
{
  // indexing
  ASSERT_TRUE(bin(0) == A(1,2));
  ASSERT_TRUE(bin(0) == A(2, 3));
  ASSERT_TRUE(bin(1) == A(6));
  bmat ref_m = "1 1 1;"
               "0 0 0;"
               "1 0 0";
  ASSERT_TRUE(ref_m == A(0, 2, 1, 3));
  bvec ref_v = "1 0 0 0";
  ASSERT_TRUE(ref_v == A.get_row(1));
  ref_m = "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A.get_rows(1, 2));
  ref_v = "1 0 0";
  ASSERT_TRUE(ref_v == A.get_col(2));
  ref_m = "1 1;"
          "0 0;"
          "0 0";
  ASSERT_TRUE(ref_m == A.get_cols(2, 3));

  // setting, copying, swapping
  bmat Mv(v);
  ref_m = "1; 1; 0; 1";
  ASSERT_TRUE(ref_m == Mv);
  bmat D(A);
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  ref_m = "1 1 1 1 0;"
          "1 0 0 0 0;"
          "1 1 0 0 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0";
  D.set_size(6, 5, true);
  ASSERT_TRUE(ref_m == D);
  D.set_size(3, 2, true);
  ref_m = "1 1;"
          "1 0;"
          "1 1";
  ASSERT_TRUE(ref_m == D);
  D.zeros();
  ref_m = "0 0;"
          "0 0;"
          "0 0";
  ASSERT_TRUE(ref_m == D);
  D.ones();
  ref_m = "1 1;"
          "1 1;"
          "1 1";
  ASSERT_TRUE(ref_m == D);
  D = A;
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D(2, 2) = c;
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D(9) = c;
  ref_m = "1 1 1 0;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.set(0, 1, c);
  ref_m = "1 0 1 0;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.set_row(1, v);
  ref_m = "1 0 1 0;"
          "1 1 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.set_col(2, u);
  ref_m = "1 0 1 0;"
          "1 1 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.set_rows(0, B.get_rows(1, 2));
  ref_m = "1 0 0 1;"
          "1 1 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.set_cols(2, B.get_cols(0, 1));
  ref_m = "1 0 0 1;"
          "1 1 1 0;"
          "1 1 1 1";
  ASSERT_TRUE(ref_m == D);
  D.copy_row(1, 2);
  ref_m = "1 0 0 1;"
          "1 1 1 1;"
          "1 1 1 1";
  ASSERT_TRUE(ref_m == D);
  D.copy_col(2, 3);
  ref_m = "1 0 1 1;"
          "1 1 1 1;"
          "1 1 1 1";
  ASSERT_TRUE(ref_m == D);
  D.swap_rows(0, 2);
  ref_m = "1 1 1 1;"
          "1 1 1 1;"
          "1 0 1 1";
  ASSERT_TRUE(ref_m == D);
  D.swap_cols(0, 3);
  ref_m = "1 1 1 1;"
          "1 1 1 1;"
          "1 0 1 1";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(1, 2, A(0, 1, 0, 1));
  ref_m = "1 1 1 1;"
          "1 1 1 1;"
          "1 0 1 0";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(0, 0, A(0, 1, 0, 1));
  ref_m = "1 1 1 1;"
          "1 0 1 1;"
          "1 0 1 0";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(1, 2, 2, 3, c);
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 0 0 0";
  ASSERT_TRUE(ref_m == D);

  // transposition
  ref_m = "1 1 1;"
          "1 0 1;"
          "1 0 0;"
          "1 0 0";
  ASSERT_TRUE(ref_m == A.T());
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A.T().T());
  ref_m = "1 1 1;"
          "1 0 1;"
          "1 0 0;"
          "1 0 0";
  ASSERT_TRUE(ref_m == A.H());

  // concatenation
  D = concat_horizontal(A, B);
  ref_m = "1 1 1 1 0 1 1 0;"
          "1 0 0 0 1 0 0 1;"
          "1 1 0 0 1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D = concat_vertical(A, B);
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0;"
          "0 1 1 0;"
          "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);

  // deleting rows, cols
  D.del_row(2);
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "0 1 1 0;"
          "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.del_rows(0, 2);
  ref_m = "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == D);
  D.del_col(3);
  ref_m = "1 0 0;"
          "1 1 0";
  ASSERT_TRUE(ref_m == D);
  D.del_cols(0, 1);
  ref_m = "0; 0";
  ASSERT_TRUE(ref_m == D);

  // inserting, appending rows cols
  bmat A2 = A;
  A2.ins_row(1, v);
  ref_m = "1 1 1 1;"
          "1 1 0 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2.ins_col(0, v);
  ref_m = "1 1 1 1 1;"
          "1 1 1 0 1;"
          "0 1 0 0 0;"
          "1 1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2.append_col(A2.get_col(3));
  ref_m = "1 1 1 1 1 1;"
          "1 1 1 0 1 0;"
          "0 1 0 0 0 0;"
          "1 1 1 0 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2.append_row(A2.get_row(0));
  ref_m = "1 1 1 1 1 1;"
          "1 1 1 0 1 0;"
          "0 1 0 0 0 0;"
          "1 1 1 0 0 0;"
          "1 1 1 1 1 1";
  ASSERT_TRUE(ref_m == A2);

  // addition
  ref_m = "1 0 0 1;"
          "0 0 0 1;"
          "0 0 0 0";
  ASSERT_TRUE(ref_m == (A + B));
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == (A + c));
  ASSERT_TRUE(ref_m == (c + A));
  A2 = A;
  A2 += B;
  ref_m = "1 0 0 1;"
          "0 0 0 1;"
          "0 0 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 += c;
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);

  // subtraction
  ref_m = "1 0 0 1;"
          "0 0 0 1;"
          "0 0 0 0";
  ASSERT_TRUE(ref_m == (A - B));
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == (A - c));
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == (c - A));
  A2 = A;
  A2 -= B;
  ref_m = "1 0 0 1;"
          "0 0 0 1;"
          "0 0 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 -= c;
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  ASSERT_TRUE(ref_m == -A);

  // multiplication
  ref_m = "0 0 1;"
          "0 1 1;"
          "0 1 0";
  ASSERT_TRUE(ref_m == A*C);
  A2 = A;
  A2 *= C;
  ref_m = "0 0 1;"
          "0 1 1;"
          "0 1 0";
  ASSERT_TRUE(ref_m == A2);
  ref_m = "0 0 0 0;"
          "0 0 0 0;"
          "0 0 0 0";
  ASSERT_TRUE(ref_m == A*c);
  ASSERT_TRUE(ref_m == c*A);
  A2 = A;
  A2 *= c;
  ASSERT_TRUE(ref_m == A2);
  ref_v = "1 1 0;";
  ASSERT_TRUE(ref_v == A*v);
  ref_m = "0 1 1 0;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == elem_mult(A, B));
  elem_mult_out(A, B, A2);
  ref_m = "0 1 1 0;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m ==A2);
  bmat B2 = B;
  elem_mult_inplace(A, B2);
  ref_m = "0 1 1 0;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == B2);
  ASSERT_TRUE(bin(1) == elem_mult_sum(A, B));

  // division
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A / c);
  A2 = A;
  A2 /= c;
  ref_m = "1 1 1 1;"
          "1 0 0 0;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 /= B;
  ref_m = "1 1 1 1;"
          "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  ref_m = "1 1 1 1;"
          "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == elem_div(A, B));
  elem_div_out(A, B, A2);
  ref_m = "1 1 1 1;"
          "1 0 0 1;"
          "1 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  ASSERT_TRUE(bin(0) == elem_div_sum(A, B));
}

static
void common_operators(const imat &A, const imat &B, const imat &C,
                      const ivec &u, const ivec &v, int c)
{
  // indexing
  ASSERT_TRUE(4 == A(1,2));
  ASSERT_TRUE(4 == A(2, 3));
  ASSERT_TRUE(6 == A(6));
  imat ref_m = "2 6 6;"
               "3 4 9;"
               "8 7 4";
  ASSERT_TRUE(ref_m == A(0, 2, 1, 3));
  ivec ref_v = "7 3 4 9";
  ASSERT_TRUE(ref_v == A.get_row(1));
  ref_m = "7 3 4 9;"
          "2 8 7 4";
  ASSERT_TRUE(ref_m == A.get_rows(1, 2));
  ref_v = "6 4 7";
  ASSERT_TRUE(ref_v == A.get_col(2));
  ref_m = "6 6;"
          "4 9;"
          "7 4";
  ASSERT_TRUE(ref_m == A.get_cols(2, 3));

  // setting, copying, swapping
  imat Mv(v);
  ref_m = "8; 5; 8; 3";
  ASSERT_TRUE(ref_m == Mv);
  imat D(A);
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 8 7 4";
  ASSERT_TRUE(ref_m == D);
  ref_m = "4 2 6 6 0;"
          "7 3 4 9 0;"
          "2 8 7 4 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0";
  D.set_size(6, 5, true);
  ASSERT_TRUE(ref_m == D);
  D.set_size(3, 2, true);
  ref_m = "4 2;"
          "7 3;"
          "2 8";
  ASSERT_TRUE(ref_m == D);
  D.zeros();
  ref_m = "0 0;"
          "0 0;"
          "0 0";
  ASSERT_TRUE(ref_m == D);
  D.ones();
  ref_m = "1 1;"
          "1 1;"
          "1 1";
  ASSERT_TRUE(ref_m == D);
  D = A;
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 8 7 4";
  ASSERT_TRUE(ref_m == D);
  D(2, 2) = c;
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 8 8 4";
  ASSERT_TRUE(ref_m == D);
  D(9) = c;
  ref_m = "4 2 6 8;"
          "7 3 4 9;"
          "2 8 8 4";
  ASSERT_TRUE(ref_m == D);
  D.set(0, 1, c);
  ref_m = "4 8 6 8;"
          "7 3 4 9;"
          "2 8 8 4";
  ASSERT_TRUE(ref_m == D);
  D.set_row(1, v);
  ref_m = "4 8 6 8;"
          "8 5 8 3;"
          "2 8 8 4";
  ASSERT_TRUE(ref_m == D);
  D.set_col(2, u);
  ref_m = "4 8 8 8;"
          "8 5 9 3;"
          "2 8 5 4";
  ASSERT_TRUE(ref_m == D);
  D.set_rows(0, B.get_rows(1, 2));
  ref_m = "5 9 8 8;"
          "6 7 8 2;"
          "2 8 5 4";
  ASSERT_TRUE(ref_m == D);
  D.set_cols(2, B.get_cols(0, 1));
  ref_m = "5 9 2 2;"
          "6 7 5 9;"
          "2 8 6 7";
  ASSERT_TRUE(ref_m == D);
  D.copy_row(1, 2);
  ref_m = "5 9 2 2;"
          "2 8 6 7;"
          "2 8 6 7";
  ASSERT_TRUE(ref_m == D);
  D.copy_col(2, 3);
  ref_m = "5 9 2 2;"
          "2 8 7 7;"
          "2 8 7 7";
  ASSERT_TRUE(ref_m == D);
  D.swap_rows(0, 2);
  ref_m = "2 8 7 7;"
          "2 8 7 7;"
          "5 9 2 2";
  ASSERT_TRUE(ref_m == D);
  D.swap_cols(0, 3);
  ref_m = "7 8 7 2;"
          "7 8 7 2;"
          "2 9 2 5";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(1, 2, A(0, 1, 0, 1));
  ref_m = "7 8 7 2;"
          "7 8 4 2;"
          "2 9 7 3";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(0, 0, A(0, 1, 0, 1));
  ref_m = "4 2 7 2;"
          "7 3 4 2;"
          "2 9 7 3";
  ASSERT_TRUE(ref_m == D);
  D.set_submatrix(1, 2, 2, 3, c);
  ref_m = "4 2 7 2;"
          "7 3 8 8;"
          "2 9 8 8";
  ASSERT_TRUE(ref_m == D);

  // transposition
  ref_m = "4 7 2;"
          "2 3 8;"
          "6 4 7;"
          "6 9 4";
  ASSERT_TRUE(ref_m == A.T());
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 8 7 4";
  ASSERT_TRUE(ref_m == A.T().T());
  ref_m = "4 7 2;"
          "2 3 8;"
          "6 4 7;"
          "6 9 4";
  ASSERT_TRUE(ref_m == A.H());

  // concatenation
  D = concat_horizontal(A, B);
  ref_m = "4 2 6 6 2 2 1 1;"
          "7 3 4 9 5 9 8 8;"
          "2 8 7 4 6 7 8 2";
  ASSERT_TRUE(ref_m == D);
  D = concat_vertical(A, B);
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 8 7 4;"
          "2 2 1 1;"
          "5 9 8 8;"
          "6 7 8 2";
  ASSERT_TRUE(ref_m == D);

  // deleting rows, cols
  D.del_row(2);
  ref_m = "4 2 6 6;"
          "7 3 4 9;"
          "2 2 1 1;"
          "5 9 8 8;"
          "6 7 8 2";
  ASSERT_TRUE(ref_m == D);
  D.del_rows(0, 2);
  ref_m = "5 9 8 8;"
          "6 7 8 2";
  ASSERT_TRUE(ref_m == D);
  D.del_col(3);
  ref_m = "5 9 8;"
          "6 7 8";
  ASSERT_TRUE(ref_m == D);
  D.del_cols(0, 1);
  ref_m = "8; 8";
  ASSERT_TRUE(ref_m == D);

  // inserting, appending rows cols
  imat A2 = A;
  A2.ins_row(1, v);
  ref_m = "4 2 6 6;"
          "8 5 8 3;"
          "7 3 4 9;"
          "2 8 7 4";
  ASSERT_TRUE(ref_m == A2);
  A2.ins_col(0, v);
  ref_m = "8 4 2 6 6;"
          "5 8 5 8 3;"
          "8 7 3 4 9;"
          "3 2 8 7 4";
  ASSERT_TRUE(ref_m == A2);
  A2.append_col(A2.get_col(3));
  ref_m = "8 4 2 6 6 6;"
          "5 8 5 8 3 8;"
          "8 7 3 4 9 4;"
          "3 2 8 7 4 7";
  ASSERT_TRUE(ref_m == A2);
  A2.append_row(A2.get_row(0));
  ref_m = "8 4 2 6 6 6;"
          "5 8 5 8 3 8;"
          "8 7 3 4 9 4;"
          "3 2 8 7 4 7;"
          "8 4 2 6 6 6";
  ASSERT_TRUE(ref_m == A2);

  // addition
  ref_m = "6 4 7 7;"
          "12 12 12 17;"
          "8 15 15 6";
  ASSERT_TRUE(ref_m == (A + B));
  ref_m = "12 10 14 14;"
          "15 11 12 17;"
          "10 16 15 12";
  ASSERT_TRUE(ref_m == (A + c));
  ASSERT_TRUE(ref_m == (c + A));
  A2 = A;
  A2 += B;
  ref_m = "6 4 7 7;"
          "12 12 12 17;"
          "8 15 15 6";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 += c;
  ref_m = "12 10 14 14;"
          "15 11 12 17;"
          "10 16 15 12";
  ASSERT_TRUE(ref_m == A2);

  // subtraction
  ref_m = "2 0 5 5;"
          "2 -6 -4 1;"
          "-4 1 -1 2";
  ASSERT_TRUE(ref_m == (A - B));
  ref_m = "-4 -6 -2 -2;"
          "-1 -5 -4 1;"
          "-6 0 -1 -4";
  ASSERT_TRUE(ref_m == (A - c));
  ref_m = "4 6 2 2;"
          "1 5 4 -1;"
          "6 0 1 4";
  ASSERT_TRUE(ref_m == (c - A));
  A2 = A;
  A2 -= B;
  ref_m = "2 0 5 5;"
          "2 -6 -4 1;"
          "-4 1 -1 2";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 -= c;
  ref_m = "-4 -6 -2 -2;"
          "-1 -5 -4 1;"
          "-6 0 -1 -4";
  ASSERT_TRUE(ref_m == A2);
  ref_m = "-4 -2 -6 -6;"
          "-7 -3 -4 -9;"
          "-2 -8 -7 -4";
  ASSERT_TRUE(ref_m == -A);

  // multiplication
  ref_m = "124 68 72;"
          "152 81 91;"
          "156 91 78";
  ASSERT_TRUE(ref_m == A*C);
  A2 = A;
  A2 *= C;
  ref_m = "124 68 72;"
          "152 81 91;"
          "156 91 78";
  ASSERT_TRUE(ref_m == A2);
  ref_m = "32 16 48 48;"
          "56 24 32 72;"
          "16 64 56 32";
  ASSERT_TRUE(ref_m == A*c);
  ASSERT_TRUE(ref_m == c*A);
  A2 = A;
  A2 *= c;
  ASSERT_TRUE(ref_m == A2);
  ref_v = "108 130 124";
  ASSERT_TRUE(ref_v == A*v);
  ref_m = "8 4 6 6;"
          "35 27 32 72;"
          "12 56 56 8";
  ASSERT_TRUE(ref_m == elem_mult(A, B));
  elem_mult_out(A, B, A2);
  ref_m = "8 4 6 6;"
          "35 27 32 72;"
          "12 56 56 8";
  ASSERT_TRUE(ref_m ==A2);
  imat B2 = B;
  elem_mult_inplace(A, B2);
  ref_m = "8 4 6 6;"
          "35 27 32 72;"
          "12 56 56 8";
  ASSERT_TRUE(ref_m == B2);
  ASSERT_TRUE(322 == elem_mult_sum(A, B));

  // division
  ref_m = "0 0 0 0;"
          "0 0 0 1;"
          "0 1 0 0";
  ASSERT_TRUE(ref_m == A / c);
  A2 = A;
  A2 /= c;
  ref_m = "0 0 0 0;"
          "0 0 0 1;"
          "0 1 0 0";
  ASSERT_TRUE(ref_m == A2);
  A2 = A;
  A2 /= B;
  ref_m = "2 1 6 6;"
          "1 0 0 1;"
          "0 1 0 2";
  ASSERT_TRUE(ref_m == A2);
  ref_m = "2 1 6 6;"
          "1 0 0 1;"
          "0 1 0 2";
  ASSERT_TRUE(ref_m == elem_div(A, B));
  elem_div_out(A, B, A2);
  ref_m = "2 1 6 6;"
          "1 0 0 1;"
          "0 1 0 2";
  ASSERT_TRUE(ref_m == A2);
  ASSERT_TRUE(20 == elem_div_sum(A, B));
}

static
void assert_vec(const vec &ref, const vec &act)
{
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref[n], act[n], tol);
  }
}
static
void assert_mat(const mat &ref, const mat &act)
{
  ASSERT_EQ(ref.rows(), act.rows());
  ASSERT_EQ(ref.cols(), act.cols());
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol);
    }
  }
}
static
void common_operators(const mat &A, const mat &B, const mat &C,
                      const vec &u, const vec &v, double c)
{
  // indexing
  ASSERT_NEAR(0.236765, A(1,2), tol);
  ASSERT_NEAR(-1.36471, A(2, 3), tol);
  ASSERT_NEAR(-0.833452, A(6), tol);
  mat ref_m = "1.21572 -0.833452 -0.0596825;"
              "-1.48074 0.236765 -0.38619;"
              "-0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, A(0, 2, 1, 3));
  vec ref_v = "-1.85156 -1.48074 0.236765 -0.38619";
  assert_vec(ref_v, A.get_row(1));
  ref_m = "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, A.get_rows(1, 2));
  ref_v = "-0.833452 0.236765 -0.864985";
  assert_vec(ref_v, A.get_col(2));
  ref_m = "-0.833452 -0.0596825;"
          "0.236765 -0.38619;"
          "-0.864985 -1.36471";
  assert_mat(ref_m, A.get_cols(2, 3));

  // setting, copying, swapping
  mat Mv(v);
  ref_m = "0.33227; 1.0172; 2.01805; 0.396482";
  assert_mat(ref_m, Mv);
  mat D(A);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, D);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825 0;"
          "-1.85156 -1.48074 0.236765 -0.38619 0;"
          "-1.91608 -0.231393 -0.864985 -1.36471 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0;"
          "0 0 0 0 0";
  D.set_size(6, 5, true);
  assert_mat(ref_m, D);
  D.set_size(3, 2, true);
  ref_m = "0.23575 1.21572;"
          "-1.85156 -1.48074;"
          "-1.91608 -0.231393";
  assert_mat(ref_m, D);
  D.zeros();
  ref_m = "0 0;"
          "0 0;"
          "0 0";
  assert_mat(ref_m, D);
  D.ones();
  ref_m = "1 1;"
          "1 1;"
          "1 1";
  assert_mat(ref_m, D);
  D = A;
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, D);
  D(2, 2) = c;
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 1.84251 -1.36471";
  assert_mat(ref_m, D);
  D(9) = c;
  ref_m = "0.23575 1.21572 -0.833452 1.84251;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 1.84251 -1.36471";
  assert_mat(ref_m, D);
  D.set(0, 1, c);
  ref_m = "0.23575 1.84251 -0.833452 1.84251;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 1.84251 -1.36471";
  assert_mat(ref_m, D);
  D.set_row(1, v);
  ref_m = "0.23575 1.84251 -0.833452 1.84251;"
          "0.33227 1.0172 2.01805 0.396482;"
          "-1.91608 -0.231393 1.84251 -1.36471";
  assert_mat(ref_m, D);
  D.set_col(2, u);
  ref_m = "0.23575 1.84251 -0.275123 1.84251;"
          "0.33227 1.0172 0.10947 0.396482;"
          "-1.91608 -0.231393 0.251788 -1.36471";
  assert_mat(ref_m, D);
  D.set_rows(0, B.get_rows(1, 2));
  ref_m = "0.542602 1.63002 -0.794633 1.1841;"
          "0.57172 0.643432 0.502798 -0.237625;"
          "-1.91608 -0.231393 0.251788 -1.36471";
  assert_mat(ref_m, D);
  D.set_cols(2, B.get_cols(0, 1));
  ref_m = "0.542602 1.63002 0.459558 -0.665486;"
          "0.57172 0.643432 0.542602 1.63002;"
          "-1.91608 -0.231393 0.57172 0.643432";
  assert_mat(ref_m, D);
  D.copy_row(1, 2);
  ref_m = "0.542602 1.63002 0.459558 -0.665486;"
          "-1.91608 -0.231393 0.57172 0.643432;"
          "-1.91608 -0.231393 0.57172 0.643432";
  assert_mat(ref_m, D);
  D.copy_col(2, 3);
  ref_m = "0.542602 1.63002 -0.665486 -0.665486;"
          "-1.91608 -0.231393 0.643432 0.643432;"
          "-1.91608 -0.231393 0.643432 0.643432";
  assert_mat(ref_m, D);
  D.swap_rows(0, 2);
  ref_m = "-1.91608 -0.231393 0.643432 0.643432;"
          "-1.91608 -0.231393 0.643432 0.643432;"
          "0.542602 1.63002 -0.665486 -0.665486";
  assert_mat(ref_m, D);
  D.swap_cols(0, 3);
  ref_m = "0.643432 -0.231393 0.643432 -1.91608;"
          "0.643432 -0.231393 0.643432 -1.91608;"
          "-0.665486 1.63002 -0.665486 0.542602";
  assert_mat(ref_m, D);
  D.set_submatrix(1, 2, A(0, 1, 0, 1));
  ref_m = "0.643432 -0.231393 0.643432 -1.91608;"
          "0.643432 -0.231393 0.23575 1.21572;"
          "-0.665486 1.63002 -1.85156 -1.48074";
  assert_mat(ref_m, D);
  D.set_submatrix(0, 0, A(0, 1, 0, 1));
  ref_m = "0.23575 1.21572 0.643432 -1.91608;"
          "-1.85156 -1.48074 0.23575 1.21572;"
          "-0.665486 1.63002 -1.85156 -1.48074";
  assert_mat(ref_m, D);
  D.set_submatrix(1, 2, 2, 3, c);
  ref_m = "0.23575 1.21572 0.643432 -1.91608;"
          "-1.85156 -1.48074 1.84251 1.84251;"
          "-0.665486 1.63002 1.84251 1.84251";
  assert_mat(ref_m, D);

  // transposition
  ref_m = "0.23575 -1.85156 -1.91608;"
          "1.21572 -1.48074 -0.231393;"
          "-0.833452 0.236765 -0.864985;"
          "-0.0596825 -0.38619 -1.36471";
  assert_mat(ref_m, A.T());
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, A.T().T());
  ref_m = "0.23575 -1.85156 -1.91608;"
          "1.21572 -1.48074 -0.231393;"
          "-0.833452 0.236765 -0.864985;"
          "-0.0596825 -0.38619 -1.36471";
  assert_mat(ref_m, A.H());

  // concatenation
  D = concat_horizontal(A, B);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825 0.459558 -0.665486 0.751416 -1.30033;"
          "-1.85156 -1.48074 0.236765 -0.38619 0.542602 1.63002 -0.794633 1.1841;"
          "-1.91608 -0.231393 -0.864985 -1.36471 0.57172 0.643432 0.502798 -0.237625";
  assert_mat(ref_m, D);
  D = concat_vertical(A, B);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471;"
          "0.459558 -0.665486 0.751416 -1.30033;"
          "0.542602 1.63002 -0.794633 1.1841;"
          "0.57172 0.643432 0.502798 -0.237625";
  assert_mat(ref_m, D);

  // deleting rows, cols
  D.del_row(2);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "0.459558 -0.665486 0.751416 -1.30033;"
          "0.542602 1.63002 -0.794633 1.1841;"
          "0.57172 0.643432 0.502798 -0.237625";
  assert_mat(ref_m, D);
  D.del_rows(0, 2);
  ref_m = "0.542602 1.63002 -0.794633 1.1841;"
          "0.57172 0.643432 0.502798 -0.237625";
  assert_mat(ref_m, D);
  D.del_col(3);
  ref_m = "0.542602 1.63002 -0.794633;"
          "0.57172 0.643432 0.502798";
  assert_mat(ref_m, D);
  D.del_cols(0, 1);
  ref_m = "-0.794633; 0.502798";
  assert_mat(ref_m, D);

  // inserting, appending rows cols
  mat A2 = A;
  A2.ins_row(1, v);
  ref_m = "0.23575 1.21572 -0.833452 -0.0596825;"
          "0.33227 1.0172 2.01805 0.396482;"
          "-1.85156 -1.48074 0.236765 -0.38619;"
          "-1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, A2);
  A2.ins_col(0, v);
  ref_m = "0.33227 0.23575 1.21572 -0.833452 -0.0596825;"
          "1.0172 0.33227 1.0172 2.01805 0.396482;"
          "2.01805 -1.85156 -1.48074 0.236765 -0.38619;"
          "0.396482 -1.91608 -0.231393 -0.864985 -1.36471";
  assert_mat(ref_m, A2);
  A2.append_col(A2.get_col(3));
  ref_m = "0.33227 0.23575 1.21572 -0.833452 -0.0596825 -0.833452;"
          "1.0172 0.33227 1.0172 2.01805 0.396482 2.01805;"
          "2.01805 -1.85156 -1.48074 0.236765 -0.38619 0.236765;"
          "0.396482 -1.91608 -0.231393 -0.864985 -1.36471 -0.864985";
  assert_mat(ref_m, A2);
  A2.append_row(A2.get_row(0));
  ref_m = "0.33227 0.23575 1.21572 -0.833452 -0.0596825 -0.833452;"
          "1.0172 0.33227 1.0172 2.01805 0.396482 2.01805;"
          "2.01805 -1.85156 -1.48074 0.236765 -0.38619 0.236765;"
          "0.396482 -1.91608 -0.231393 -0.864985 -1.36471 -0.864985;"
          "0.33227 0.23575 1.21572 -0.833452 -0.0596825 -0.833452";
  assert_mat(ref_m, A2);

  // addition
  ref_m = "0.695308 0.550239 -0.0820359 -1.36001;"
          "-1.30896 0.149281 -0.557868 0.797913;"
          "-1.34436 0.412039 -0.362187 -1.60234";
  assert_mat(ref_m, (A + B));
  ref_m = "2.07826 3.05823 1.00906 1.78283;"
          "-0.00904911 0.361768 2.07928 1.45632;"
          "-0.0735664 1.61112 0.977525 0.4778";
  assert_mat(ref_m, (A + c));
  assert_mat(ref_m, (c + A));
  A2 = A;
  A2 += B;
  ref_m = "0.695308 0.550239 -0.0820359 -1.36001;"
          "-1.30896 0.149281 -0.557868 0.797913;"
          "-1.34436 0.412039 -0.362187 -1.60234";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 += c;
  ref_m = "2.07826 3.05823 1.00906 1.78283;"
          "-0.00904911 0.361768 2.07928 1.45632;"
          "-0.0735664 1.61112 0.977525 0.4778";
  assert_mat(ref_m, A2);

  // subtraction
  ref_m = "-0.223808 1.88121 -1.58487 1.24064;"
          "-2.39416 -3.11076 1.0314 -1.57029;"
          "-2.4878 -0.874826 -1.36778 -1.12709";
  assert_mat(ref_m, (A - B));
  ref_m = "-1.60676 -0.626786 -2.67596 -1.90219;"
          "-3.69407 -3.32325 -1.60574 -2.2287;"
          "-3.75859 -2.0739 -2.7075 -3.20722";
  assert_mat(ref_m, (A - c));
  ref_m = "1.60676 0.626786 2.67596 1.90219;"
          "3.69407 3.32325 1.60574 2.2287;"
          "3.75859 2.0739 2.7075 3.20722";
  assert_mat(ref_m, (c - A));
  A2 = A;
  A2 -= B;
  ref_m = "-0.223808 1.88121 -1.58487 1.24064;"
          "-2.39416 -3.11076 1.0314 -1.57029;"
          "-2.4878 -0.874826 -1.36778 -1.12709";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 -= c;
  ref_m = "-1.60676 -0.626786 -2.67596 -1.90219;"
          "-3.69407 -3.32325 -1.60574 -2.2287;"
          "-3.75859 -2.0739 -2.7075 -3.20722";
  assert_mat(ref_m, A2);
  ref_m = "-0.23575 -1.21572 0.833452 0.0596825;"
          "1.85156 1.48074 -0.236765 0.38619;"
          "1.91608 0.231393 0.864985 1.36471";
  assert_mat(ref_m, -A);

  // multiplication
  ref_m = "-2.54149 -1.07673 -0.293434;"
          "2.36426 1.1283 2.94484;"
          "-2.00345 -1.34274 2.18407";
  assert_mat(ref_m, A*C);
  A2 = A;
  A2 *= C;
  ref_m = "-2.54149 -1.07673 -0.293434;"
          "2.36426 1.1283 2.94484;"
          "-2.00345 -1.34274 2.18407";
  assert_mat(ref_m, A2);
  ref_m = "0.434372 2.23998 -1.53564 -0.109966;"
          "-3.41152 -2.72828 0.436242 -0.711559;"
          "-3.53039 -0.426344 -1.59374 -2.51449";
  assert_mat(ref_m, A*c);
  assert_mat(ref_m, c*A);
  A2 = A;
  A2 *= c;
  ref_m = "0.434372 2.23998 -1.53564 -0.109966;"
          "-3.41152 -2.72828 0.436242 -0.711559;"
          "-3.53039 -0.426344 -1.59374 -2.51449";
  assert_mat(ref_m, A2);
  ref_v = "-0.390636 -1.79675 -3.15869";
  assert_vec(ref_v, A*v);
  ref_m = "0.108341 -0.809047 -0.626268 0.0776068;"
          "-1.00466 -2.41364 -0.188141 -0.457289;"
          "-1.09546 -0.148886 -0.434913 0.324289";
  assert_mat(ref_m, elem_mult(A, B));
  elem_mult_out(A, B, A2);
  ref_m = "0.108341 -0.809047 -0.626268 0.0776068;"
          "-1.00466 -2.41364 -0.188141 -0.457289;"
          "-1.09546 -0.148886 -0.434913 0.324289";
  assert_mat(ref_m, A2);
  mat B2 = B;
  elem_mult_inplace(A, B2);
  ref_m = "0.108341 -0.809047 -0.626268 0.0776068;"
          "-1.00466 -2.41364 -0.188141 -0.457289;"
          "-1.09546 -0.148886 -0.434913 0.324289";
  assert_mat(ref_m, B2);
  ASSERT_NEAR(-6.66807, elem_mult_sum(A, B), tol);

  // division
  ref_m = "0.127951 0.65982 -0.452346 -0.032392;"
          "-1.00491 -0.803655 0.128501 -0.2096;"
          "-1.03993 -0.125586 -0.46946 -0.74068";
  assert_mat(ref_m, A / c);
  A2 = A;
  A2 /= c;
  ref_m = "0.127951 0.65982 -0.452346 -0.032392;"
          "-1.00491 -0.803655 0.128501 -0.2096;"
          "-1.03993 -0.125586 -0.46946 -0.74068";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 /= B;
  ref_m = "0.512993 -1.82682 -1.10918 0.0458981;"
          "-3.41237 -0.908418 -0.297955 -0.326145;"
          "-3.35142 -0.359623 -1.72034 5.74313";
  assert_mat(ref_m, A2);
  ref_m = "0.512993 -1.82682 -1.10918 0.0458981;"
          "-3.41237 -0.908418 -0.297955 -0.326145;"
          "-3.35142 -0.359623 -1.72034 5.74313";
  assert_mat(ref_m, elem_div(A, B));
  elem_div_out(A, B, A2);
  ref_m = "0.512993 -1.82682 -1.10918 0.0458981;"
          "-3.41237 -0.908418 -0.297955 -0.326145;"
          "-3.35142 -0.359623 -1.72034 5.74313";
  assert_mat(ref_m, A2);
  ASSERT_NEAR(-7.01026, elem_div_sum(A, B), tol);
}

static
void assert_vec(const cvec &ref, const cvec &act)
{
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref[n].real(), act[n].real(), tol);
    ASSERT_NEAR(ref[n].imag(), act[n].imag(), tol);
  }
}
static
void assert_mat_p(const cmat &ref, const cmat &act, int line)
{
  ASSERT_EQ(ref.rows(), act.rows()) << "line:" << line;
  ASSERT_EQ(ref.cols(), act.cols()) << "line:" << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k).real(), act(n,k).real(), tol) << "line:" << line;
      ASSERT_NEAR(ref(n,k).imag(), act(n,k).imag(), tol) << "line:" << line;
    }
  }
}
#define assert_mat(ref, act) assert_mat_p(ref, act, __LINE__)
static
void common_operators(const cmat &A, const cmat &B, const cmat &C,
                      const cvec &u, const cvec &v, complex<double> c)
{
  // indexing
  ASSERT_NEAR(-0.595499, A(1,2).real(), tol);
  ASSERT_NEAR(-0.56388, A(1,2).imag(), tol);
  ASSERT_NEAR(-0.817562, A(2, 3).real(), tol);
  ASSERT_NEAR(-1.10836, A(2, 3).imag(), tol);
  ASSERT_NEAR(0.722271, A(6).real(), tol);
  ASSERT_NEAR(-0.153507, A(6).imag(), tol);
  cmat ref_m = "-1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
              "0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
              "0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A(0, 2, 1, 3));
  cvec ref_v = "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i";
  assert_vec(ref_v, A.get_row(1));
  ref_m = "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A.get_rows(1, 2));
  ref_v = "0.722271-0.153507i -0.595499-0.56388i 1.0276-0.175089i";
  assert_vec(ref_v, A.get_col(2));
  ref_m = "0.722271-0.153507i -0.0211019-0.365373i;"
          "-0.595499-0.56388i -0.687348-0.0869053i;"
          "1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A.get_cols(2, 3));

  // setting, copying, swapping
  cmat Mv(v);
  ref_m = "0.174329+0.0466934i; -0.51944+0.527715i; -0.410049-0.311715i; -0.253286-0.813752i";
  assert_mat(ref_m, Mv);
  cmat D(A);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i 0+0i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i 0+0i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i 0+0i;"
          "0+0i 0+0i 0+0i 0+0i 0+0i;"
          "0+0i 0+0i 0+0i 0+0i 0+0i;"
          "0+0i 0+0i 0+0i 0+0i 0+0i";
  D.set_size(6, 5, true);
  assert_mat(ref_m, D);
  D.set_size(3, 2, true);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i;"
          "0.0690522+0.282504i 0.308802+0.109446i;"
          "0.570089-0.350927i 0.248223-1.13728i";
  assert_mat(ref_m, D);
  D.zeros();
  ref_m = "0+0i 0+0i;"
          "0+0i 0+0i;"
          "0+0i 0+0i";
  assert_mat(ref_m, D);
  D.ones();
  ref_m = "1+0i 1+0i;"
          "1+0i 1+0i;"
          "1+0i 1+0i";
  assert_mat(ref_m, D);
  D = A;
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D(2, 2) = c;
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.117174+0.468096i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D(9) = c;
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.117174+0.468096i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.117174+0.468096i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D.set(0, 1, c);
  ref_m = "-1.0992+0.0482272i -0.117174+0.468096i 0.722271-0.153507i -0.117174+0.468096i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.117174+0.468096i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D.set_row(1, v);
  ref_m = "-1.0992+0.0482272i -0.117174+0.468096i 0.722271-0.153507i -0.117174+0.468096i;"
          "0.174329+0.0466934i -0.51944+0.527715i -0.410049-0.311715i -0.253286-0.813752i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.117174+0.468096i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D.set_col(2, u);
  ref_m = "-1.0992+0.0482272i -0.117174+0.468096i 0.455057-0.444672i -0.117174+0.468096i;"
          "0.174329+0.0466934i -0.51944+0.527715i -0.8556-0.192441i -0.253286-0.813752i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.110127+1.22645i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D.set_rows(0, B.get_rows(1, 2));
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i 1.21948+1.39547i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i 0.0988444-0.333968i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.110127+1.22645i -0.817562-1.10836i";
  assert_mat(ref_m, D);
  D.set_cols(2, B.get_cols(0, 1));
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i -0.940531+1.37009i -0.287338+0.801736i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.221685-0.581717i -0.647685-0.916925i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.537457-0.362192i -0.221416-0.892514i";
  assert_mat(ref_m, D);
  D.copy_row(1, 2);
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i -0.940531+1.37009i -0.287338+0.801736i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.537457-0.362192i -0.221416-0.892514i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.537457-0.362192i -0.221416-0.892514i";
  assert_mat(ref_m, D);
  D.copy_col(2, 3);
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i -0.287338+0.801736i -0.287338+0.801736i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.221416-0.892514i -0.221416-0.892514i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.221416-0.892514i -0.221416-0.892514i";
  assert_mat(ref_m, D);
  D.swap_rows(0, 2);
  ref_m = "0.570089-0.350927i 0.248223-1.13728i -0.221416-0.892514i -0.221416-0.892514i;"
          "0.570089-0.350927i 0.248223-1.13728i -0.221416-0.892514i -0.221416-0.892514i;"
          "-0.221685-0.581717i -0.647685-0.916925i -0.287338+0.801736i -0.287338+0.801736i";
  assert_mat(ref_m, D);
  D.swap_cols(0, 3);
  ref_m = "-0.221416-0.892514i 0.248223-1.13728i -0.221416-0.892514i 0.570089-0.350927i;"
          "-0.221416-0.892514i 0.248223-1.13728i -0.221416-0.892514i 0.570089-0.350927i;"
          "-0.287338+0.801736i -0.647685-0.916925i -0.287338+0.801736i -0.221685-0.581717i";
  assert_mat(ref_m, D);
  D.set_submatrix(1, 2, A(0, 1, 0, 1));
  ref_m = "-0.221416-0.892514i 0.248223-1.13728i -0.221416-0.892514i 0.570089-0.350927i;"
          "-0.221416-0.892514i 0.248223-1.13728i -1.0992+0.0482272i -1.44302-0.163398i;"
          "-0.287338+0.801736i -0.647685-0.916925i 0.0690522+0.282504i 0.308802+0.109446i";
  assert_mat(ref_m, D);
  D.set_submatrix(0, 0, A(0, 1, 0, 1));
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i -0.221416-0.892514i 0.570089-0.350927i;"
          "0.0690522+0.282504i 0.308802+0.109446i -1.0992+0.0482272i -1.44302-0.163398i;"
          "-0.287338+0.801736i -0.647685-0.916925i 0.0690522+0.282504i 0.308802+0.109446i";
  assert_mat(ref_m, D);
  D.set_submatrix(1, 2, 2, 3, c);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i -0.221416-0.892514i 0.570089-0.350927i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.117174+0.468096i -0.117174+0.468096i;"
          "-0.287338+0.801736i -0.647685-0.916925i -0.117174+0.468096i -0.117174+0.468096i";
  assert_mat(ref_m, D);

  // transposition
  ref_m = "-1.0992+0.0482272i 0.0690522+0.282504i 0.570089-0.350927i;"
          "-1.44302-0.163398i 0.308802+0.109446i 0.248223-1.13728i;"
          "0.722271-0.153507i -0.595499-0.56388i 1.0276-0.175089i;"
          "-0.0211019-0.365373i -0.687348-0.0869053i -0.817562-1.10836i";
  assert_mat(ref_m, A.T());
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A.T().T());
  ref_m = "-1.0992-0.0482272i 0.0690522-0.282504i 0.570089+0.350927i;"
          "-1.44302+0.163398i 0.308802-0.109446i 0.248223+1.13728i;"
          "0.722271+0.153507i -0.595499+0.56388i 1.0276+0.175089i;"
          "-0.0211019+0.365373i -0.687348+0.0869053i -0.817562+1.10836i";
  assert_mat(ref_m, A.H());

  // concatenation
  D = concat_horizontal(A, B);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i -0.940531+1.37009i -0.287338+0.801736i -0.332921+0.669844i -0.331418+0.133073i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i -0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i 1.21948+1.39547i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i -0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i 0.0988444-0.333968i";
  assert_mat(ref_m, D);
  D = concat_vertical(A, B);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i;"
          "-0.940531+1.37009i -0.287338+0.801736i -0.332921+0.669844i -0.331418+0.133073i;"
          "-0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i 1.21948+1.39547i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i 0.0988444-0.333968i";
  assert_mat(ref_m, D);

  // deleting rows, cols
  D.del_row(2);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "-0.940531+1.37009i -0.287338+0.801736i -0.332921+0.669844i -0.331418+0.133073i;"
          "-0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i 1.21948+1.39547i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i 0.0988444-0.333968i";
  assert_mat(ref_m, D);
  D.del_rows(0, 2);
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i 1.21948+1.39547i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i 0.0988444-0.333968i";
  assert_mat(ref_m, D);
  D.del_col(3);
  ref_m = "-0.221685-0.581717i -0.647685-0.916925i 0.459866+0.264737i;"
          "-0.537457-0.362192i -0.221416-0.892514i -0.337512-0.476078i";
  assert_mat(ref_m, D);
  D.del_cols(0, 1);
  ref_m = "0.459866+0.264737i; -0.337512-0.476078i";
  assert_mat(ref_m, D);

  // inserting, appending rows cols
  cmat A2 = A;
  A2.ins_row(1, v);
  ref_m = "-1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "0.174329+0.0466934i -0.51944+0.527715i -0.410049-0.311715i -0.253286-0.813752i;"
          "0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A2);
  A2.ins_col(0, v);
  ref_m = "0.174329+0.0466934i -1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i;"
          "-0.51944+0.527715i 0.174329+0.0466934i -0.51944+0.527715i -0.410049-0.311715i -0.253286-0.813752i;"
          "-0.410049-0.311715i 0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i;"
          "-0.253286-0.813752i 0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i";
  assert_mat(ref_m, A2);
  A2.append_col(A2.get_col(3));
  ref_m = "0.174329+0.0466934i -1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i 0.722271-0.153507i;"
          "-0.51944+0.527715i 0.174329+0.0466934i -0.51944+0.527715i -0.410049-0.311715i -0.253286-0.813752i -0.410049-0.311715i;"
          "-0.410049-0.311715i 0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i -0.595499-0.56388i;"
          "-0.253286-0.813752i 0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i 1.0276-0.175089i";
  assert_mat(ref_m, A2);
  A2.append_row(A2.get_row(0));
  ref_m = "0.174329+0.0466934i -1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i 0.722271-0.153507i;"
          "-0.51944+0.527715i 0.174329+0.0466934i -0.51944+0.527715i -0.410049-0.311715i -0.253286-0.813752i -0.410049-0.311715i;"
          "-0.410049-0.311715i 0.0690522+0.282504i 0.308802+0.109446i -0.595499-0.56388i -0.687348-0.0869053i -0.595499-0.56388i;"
          "-0.253286-0.813752i 0.570089-0.350927i 0.248223-1.13728i 1.0276-0.175089i -0.817562-1.10836i 1.0276-0.175089i;"
          "0.174329+0.0466934i -1.0992+0.0482272i -1.44302-0.163398i 0.722271-0.153507i -0.0211019-0.365373i 0.722271-0.153507i";
  assert_mat(ref_m, A2);

  // addition
  ref_m = "-2.03973+1.41831i -1.73036+0.638337i 0.38935+0.516337i -0.35252-0.232301i;"
          "-0.152633-0.299213i -0.338883-0.807479i -0.135633-0.299143i 0.532135+1.30856i;"
          "0.0326323-0.71312i 0.0268068-2.02979i 0.690083-0.651167i -0.718717-1.44232i";
  assert_mat(ref_m, (A + B));
  ref_m = "-1.21637+0.516323i -1.5602+0.304697i 0.605097+0.314588i -0.138276+0.102722i;"
          "-0.0481216+0.750599i 0.191628+0.577542i -0.712673-0.0957844i -0.804522+0.38119i;"
          "0.452915+0.117168i 0.13105-0.669184i 0.910422+0.293007i -0.934736-0.64026i";
  assert_mat(ref_m, (A + c));
  assert_mat(ref_m, (c + A));
  A2 = A;
  A2 += B;
  ref_m = "-2.03973+1.41831i -1.73036+0.638337i 0.38935+0.516337i -0.35252-0.232301i;"
          "-0.152633-0.299213i -0.338883-0.807479i -0.135633-0.299143i 0.532135+1.30856i;"
          "0.0326323-0.71312i 0.0268068-2.02979i 0.690083-0.651167i -0.718717-1.44232i";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 += c;
  ref_m = "-1.21637+0.516323i -1.5602+0.304697i 0.605097+0.314588i -0.138276+0.102722i;"
          "-0.0481216+0.750599i 0.191628+0.577542i -0.712673-0.0957844i -0.804522+0.38119i;"
          "0.452915+0.117168i 0.13105-0.669184i 0.910422+0.293007i -0.934736-0.64026i";
  assert_mat(ref_m, A2);

  // subtraction
  ref_m = "-0.158667-1.32186i -1.15569-0.965134i 1.05519-0.823352i 0.310316-0.498446i;"
          "0.290738+0.864221i 0.956487+1.02637i -1.05537-0.828617i -1.90683-1.48237i;"
          "1.10755+0.0112649i 0.46964-0.244766i 1.36511+0.30099i -0.916406-0.774388i";
  assert_mat(ref_m, (A - B));
  ref_m = "-0.982024-0.419868i -1.32585-0.631494i 0.839445-0.621603i 0.0960718-0.833469i;"
          "0.186226-0.185592i 0.425976-0.35865i -0.478326-1.03198i -0.570174-0.555001i;"
          "0.687263-0.819023i 0.365397-1.60538i 1.14477-0.643184i -0.700388-1.57645i";
  assert_mat(ref_m, (A - c));
  ref_m = "0.982024+0.419868i 1.32585+0.631494i -0.839445+0.621603i -0.0960718+0.833469i;"
          "-0.186226+0.185592i -0.425976+0.35865i 0.478326+1.03198i 0.570174+0.555001i;"
          "-0.687263+0.819023i -0.365397+1.60538i -1.14477+0.643184i 0.700388+1.57645i";
  assert_mat(ref_m, (c - A));
  A2 = A;
  A2 -= B;
  ref_m = "-0.158667-1.32186i -1.15569-0.965134i 1.05519-0.823352i 0.310316-0.498446i;"
          "0.290738+0.864221i 0.956487+1.02637i -1.05537-0.828617i -1.90683-1.48237i;"
          "1.10755+0.0112649i 0.46964-0.244766i 1.36511+0.30099i -0.916406-0.774388i";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 -= c;
  ref_m = "-0.982024-0.419868i -1.32585-0.631494i 0.839445-0.621603i 0.0960718-0.833469i;"
          "0.186226-0.185592i 0.425976-0.35865i -0.478326-1.03198i -0.570174-0.555001i;"
          "0.687263-0.819023i 0.365397-1.60538i 1.14477-0.643184i -0.700388-1.57645i";
  assert_mat(ref_m, A2);
  ref_m = "1.0992-0.0482272i 1.44302+0.163398i -0.722271+0.153507i 0.0211019+0.365373i;"
          "-0.0690522-0.282504i -0.308802-0.109446i 0.595499+0.56388i 0.687348+0.0869053i;"
          "-0.570089+0.350927i -0.248223+1.13728i -1.0276+0.175089i 0.817562+1.10836i";
  assert_mat(ref_m, -A);

  // multiplication
  ref_m = "-0.910805-1.18033i -1.44173-0.454669i -1.57228+0.766807i;"
          "1.10653+0.339386i -0.365877+0.320608i 1.006+0.9153i;"
          "2.11969+0.360853i -1.86799-1.10692i -1.79365+0.0201111i";
  assert_mat(ref_m, A*C);
  A2 = A;
  A2 *= C;
  ref_m = "-0.910805-1.18033i -1.44173-0.454669i -1.57228+0.766807i;"
          "1.10653+0.339386i -0.365877+0.320608i 1.006+0.9153i;"
          "2.11969+0.360853i -1.86799-1.10692i -1.79365+0.0201111i";
  assert_mat(ref_m, A2);
  ref_m = "0.106222-0.520181i 0.245571-0.656327i -0.0127751+0.356079i 0.173502+0.0329345i;"
          "-0.14033-0.000778999i -0.0874147+0.131725i 0.333727-0.212679i 0.121219-0.311562i;"
          "0.0974682+0.307976i 0.50327+0.249452i -0.038449+0.501529i 0.614613-0.252827i";
  assert_mat(ref_m, A*c);
  assert_mat(ref_m, c*A);
  A2 = A;
  A2 *= c;
  ref_m = "0.106222-0.520181i 0.245571-0.656327i -0.0127751+0.356079i 0.173502+0.0329345i;"
          "-0.14033-0.000778999i -0.0874147+0.131725i 0.333727-0.212679i 0.121219-0.311562i;"
          "0.0974682+0.307976i 0.50327+0.249452i -0.038449+0.501529i 0.614613-0.252827i";
  assert_mat(ref_m, A2);
  ref_v = "0.00592265-0.772031i -0.0475239+1.15677i -0.5838+1.38468i";
  assert_vec(ref_v, A*v);
  ref_m = "0.967754-1.55136i 0.545637-1.10997i -0.137633+0.534915i 0.0556148+0.118283i;"
          "0.149029-0.102796i -0.0996524-0.354035i -0.12457-0.41696i -0.716936-1.06515i;"
          "-0.433501-0.0178736i -1.07+0.0302697i -0.430182-0.430121i -0.450966+0.163484i";
  assert_mat(ref_m, elem_mult(A, B));
  elem_mult_out(A, B, A2);
  ref_m = "0.967754-1.55136i 0.545637-1.10997i -0.137633+0.534915i 0.0556148+0.118283i;"
          "0.149029-0.102796i -0.0996524-0.354035i -0.12457-0.41696i -0.716936-1.06515i;"
          "-0.433501-0.0178736i -1.07+0.0302697i -0.430182-0.430121i -0.450966+0.163484i";
  assert_mat(ref_m, A2);
  cmat B2 = B;
  elem_mult_inplace(A, B2);
  ref_m = "0.967754-1.55136i 0.545637-1.10997i -0.137633+0.534915i 0.0556148+0.118283i;"
          "0.149029-0.102796i -0.0996524-0.354035i -0.12457-0.41696i -0.716936-1.06515i;"
          "-0.433501-0.0178736i -1.07+0.0302697i -0.430182-0.430121i -0.450966+0.163484i";
  assert_mat(ref_m, B2);
  complex<double> ref_s = elem_mult_sum(A, B);
  ASSERT_NEAR(-1.7454, ref_s.real(), tol);
  ASSERT_NEAR(-4.20131, ref_s.imag(), tol);

  // division
  ref_m = "0.650103+2.1855i 0.397686+2.98321i -0.672071-1.37477i -0.723908+0.226289i;"
          "0.533181-0.280983i 0.0646261-0.675875i -0.833921+1.48092i 0.171184+1.42554i;"
          "-0.992372-0.969479i -2.41124+0.0732986i -0.869106-1.97771i -1.81676+2.20134i";
  assert_mat(ref_m, A / c);
  A2 = A;
  A2 /= c;
  ref_m = "0.650103+2.1855i 0.397686+2.98321i -0.672071-1.37477i -0.723908+0.226289i;"
          "0.533181-0.280983i 0.0646261-0.675875i -0.833921+1.48092i 0.171184+1.42554i;"
          "-0.992372-0.969479i -2.41124+0.0732986i -0.869106-1.97771i -1.81676+2.20134i";
  assert_mat(ref_m, A2);
  A2 = A;
  A2 /= B;
  ref_m = "0.398266+0.528884i 0.391033+1.65973i -0.613527-0.773337i -0.326374+0.971408i;"
          "-0.463553-0.0579506i -0.238334+0.168429i -1.50279-0.361051i -0.279369+0.248421i;"
          "-0.426849+0.940595i 1.13537+0.559782i -0.773626+1.61i 2.38526-3.154i";
  assert_mat(ref_m, A2);
  ref_m = "0.398266+0.528884i 0.391033+1.65973i -0.613527-0.773337i -0.326374+0.971408i;"
          "-0.463553-0.0579506i -0.238334+0.168429i -1.50279-0.361051i -0.279369+0.248421i;"
          "-0.426849+0.940595i 1.13537+0.559782i -0.773626+1.61i 2.38526-3.154i";
  assert_mat(ref_m, elem_div(A, B));
  elem_div_out(A, B, A2);
  ref_m = "0.398266+0.528884i 0.391033+1.65973i -0.613527-0.773337i -0.326374+0.971408i;"
          "-0.463553-0.0579506i -0.238334+0.168429i -1.50279-0.361051i -0.279369+0.248421i;"
          "-0.426849+0.940595i 1.13537+0.559782i -0.773626+1.61i 2.38526-3.154i";
  assert_mat(ref_m, A2);
  ref_s = elem_div_sum(A, B);
  ASSERT_NEAR(-0.314489, ref_s.real(), tol);
  ASSERT_NEAR(2.34092, ref_s.imag(), tol);
}

TEST(Mat, All)
{
  RNG_reset(0);

  // Testing Mat<bin> (bmat)
  bmat bM1 = randb(3, 4);
  bmat bM2 = randb(3, 4);
  bmat bM3 = randb(4, 3);
  bvec bv1 = randb(3);
  bvec bv2 = randb(4);
  bin bx = randb();
  common_operators(bM1, bM2, bM3, bv1, bv2, bx);

  // Testing Mat<int> (imat)
  imat iM1 = randi(3, 4, 1, 9);
  imat iM2 = randi(3, 4, 1, 9);
  imat iM3 = randi(4, 3, 1, 9);
  ivec iv1 = randi(3, 1, 9);
  ivec iv2 = randi(4, 1, 9);
  int ix = randi(1, 9);
  common_operators(iM1, iM2, iM3, iv1, iv2, ix);

  // Testing Mat<double> (mat)
  mat dM1 = randn(3, 4);
  mat dM2 = randn(3, 4);
  mat dM3 = randn(4, 3);
  vec dv1 = randn(3);
  vec dv2 = randn(4);
  double dx = randn();
  common_operators(dM1, dM2, dM3, dv1, dv2, dx);

  // Testing Mat<complex<double> > (cmat)
  cmat cM1 = randn_c(3, 4);
  cmat cM2 = randn_c(3, 4);
  cmat cM3 = randn_c(4, 3);
  cvec cv1 = randn_c(3);
  cvec cv2 = randn_c(4);
  complex<double> cx = randn_c();
  common_operators(cM1, cM2, cM3, cv1, cv2, cx);

  //initializations with strings are tested implicitly
}
