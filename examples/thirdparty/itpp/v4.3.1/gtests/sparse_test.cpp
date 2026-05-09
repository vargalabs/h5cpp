/*!
 * \file
 * \brief Sparse vectors and matrices test program
 * \author 
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

static const double tol = 1e-4;

static
void assert_vec_p(const vec &ref, const vec &act, int line)
{
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol) << line;
  }
}
#define assert_vec(ref, act) assert_vec_p(ref, act, __LINE__)

static
void assert_mat_p(const mat &ref, const mat &act, int line)
{
  ASSERT_EQ(ref.rows(), act.rows()) << line;
  ASSERT_EQ(ref.cols(), act.cols()) << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol) << line;
    }
  }
}
#define assert_mat(ref, act) assert_mat_p(ref, act, __LINE__)

TEST(Sparse, All)
{
  RNG_reset(0);

  //********* Testing Sparse_Vec *********

  Sparse_Vec<double> v1(5), v2(5);
  Sparse_Vec<double> v_out(5);
  ivec index_vec;
  vec v;

  v1.set(1, 3);
  v2.set(1, 1);
  v2.set(2, 2);
  vec ref = "0 3 0 0 0";
  assert_vec(ref, v1.full());
  ref = "0 1 2 0 0";
  assert_vec(ref, v2.full());
  v_out = v1 + v2;
  ref = "0 4 2 0 0";
  assert_vec(ref, v_out.full());
  ASSERT_NEAR(3, v1*v2, tol);
  ASSERT_NEAR(0.2, v1.density(), tol);
  v1.zeros();
  ref = "0 0 0 0 0";
  assert_vec(ref, v1.full());

  index_vec = "0 2 2 3";
  v = "1 2 4 3";
  v1.set(index_vec, v);
  ref = "1 0 6 3 0";
  assert_vec(ref, v1.full());
  v1.set_new(index_vec, v);
  // Unnoticed error in set_new() if the same index is used several times
  ref = "1 0 4 3 0";
  assert_vec(ref, v1.full());

  v1.zeros();
  v1 -= v2;
  ref = "0 -1 -2 0 0";
  assert_vec(ref, v1.full());
  v1 /= 4;
  ref = "0 -0.25 -0.5 0 0";
  assert_vec(ref, v1.full());
  v1 *= 2;
  ref = "0 -0.5 -1 0 0";
  assert_vec(ref, v1.full());
  v2 /= 2;
  v1 += v2;
  ref = "0 0 0 0 0";
  assert_vec(ref, v1.full());
  ASSERT_EQ(0, v1.nnz());

  index_vec = "0 2 2 1";
  v = "-1 5 -2 4";
  v1.add(index_vec, v);
  ref = "-1 4 3 0 0";
  assert_vec(ref, v1.full());

  v1.clear_elem(2);
  ref = "-1 4 0 0 0";
  assert_vec(ref, v1.full());


  Sparse_Vec<double> v3;
  v3.set_size(2);
  v3.set(0, 3);
  ref = "3 0";
  assert_vec(ref, v3.full());
  v3.set_size(5);
  ref = "0 0 0 0 0";
  assert_vec(ref, v3.full());

  Sparse_Vec<double> v4(5), v5(5);
  v4.set(1, 1);
  v5.set(1, 1);
  v4.set(2, 2);
  v5.set(2, 4);
  ref = "0 1 2 0 0";
  assert_vec(ref, v4.full());
  ref = "0 1 4 0 0";
  assert_vec(ref, v5.full());

  //********* Testing Sparse_Mat *********
  Sparse_Mat<double> m1(3, 3), m2(3, 3);

  m1.set(1, 1, 3);
  m2.set(1, 2, 1);
  m2.set(2, 1, 2);
  mat ref_m = "0 0 0; 0 3 0; 0 0 0";
  assert_mat(ref_m, full(m1));
  ref_m = "0 0 0; 0 0 1; 0 2 0";
  assert_mat(ref_m, full(m2));
  ref_m = "0 0 0; 0 3 1; 0 2 0";
  assert_mat(ref_m, full(m1 + m2));
  ref_m = "0 0 0; 0 0 3; 0 0 0";
  assert_mat(ref_m, full(m1*m2));
  ASSERT_NEAR(0.111111, m1.density(), tol);
  ref_m = "0 0 0; 0 0 2; 0 1 0";
  assert_mat(ref_m, full(transpose(m2)));

  m1.zeros();
  ref_m = "0 0 0; 0 0 0; 0 0 0";
  assert_mat(ref_m, full(m1));
  m1 -= m2;
  ref_m = "0 0 0; 0 0 -1; 0 -2 0";
  assert_mat(ref_m, full(m1));
  m1 /= 4;
  ref_m = "0 0 0; 0 0 -0.25; 0 -0.5 0";
  assert_mat(ref_m, full(m1));
  m1 *= 2;
  ref_m = "0 0 0; 0 0 -0.5; 0 -1 0";
  assert_mat(ref_m, full(m1));

  m1.add_elem(0, 2, 4);
  ref_m = "0 0 4; 0 0 -0.5; 0 -1 0";
  assert_mat(ref_m, full(m1));
  m1.clear_elem(1, 2);
  ref_m = "0 0 4; 0 0 0; 0 -1 0";
  assert_mat(ref_m, full(m1));

  Sparse_Mat<double> m3(2, 3), m4(3, 2);

  m3.set(0, 0, 3);
  m4.set(0, 1, 1);
  m4.set(2, 0, 2);
  ref_m = "3 0 0; 0 0 0";
  assert_mat(ref_m, full(m3));
  ref_m = "0 1; 0 0; 2 0";
  assert_mat(ref_m, full(m4));
  ref_m = "0 3; 0 0";
  assert_mat(ref_m, full(m3*m4));
  ref_m = "0 0 2; 1 0 0";
  assert_mat(ref_m, full(transpose(m4)));
  ref_m = "4 0; 0 1";
  assert_mat(ref_m, full(trans_mult(m4, m4)));
  ref_m = "1 0 0; 0 0 0; 0 0 4";
  assert_mat(ref_m, full(mult_trans(m4, m4)));
  ref_m = "1 0 0; 0 0 0; 0 0 4";
  assert_mat(ref_m, full(mult_trans(m4, m4)));

  v = "1 2 3";
  ref = "3 0";
  assert_vec(ref, m3*v);

  Sparse_Mat<double> A(3, 5), B(5, 3), C;
  vec x(3), y(5), z1, z2;

  A = randn(3, 5);
  B = randn(5, 3);
  x = randn(3);
  y = randn(5);

  C = A * B;
  z1 = A * y;
  z2 = x * A;
  
  ref_m = "0.0511847 1.39674 2.48276;"
          "0.916172 -1.00692 -1.16731;"
          "-1.21703 -4.20494 -4.17853";
  assert_mat(ref_m, full(C));
  ref = "-1.31897 -1.49449 5.33006";
  assert_vec(ref, z1);
  ref = "0.123808 -0.590395 0.45136 -1.52973 0.00504626";
  assert_vec(ref, z2);
}
