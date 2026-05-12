/*!
 * \file
 * \brief Test program of various functions on vectors and matrices
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/itstat.h>
#include "gtest/gtest.h"

using namespace itpp;

static const double tol = 1e-4;

static
void assert_vec(const vec &ref, const vec &act)
{
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n)
  {
    ASSERT_NEAR(ref[n], act[n], tol);
  }
}
static
void assert_mat(const mat &ref, const mat &act)
{
  ASSERT_EQ(ref.rows(), act.rows());
  ASSERT_EQ(ref.cols(), act.cols());
  for (int n = 0; n < ref.rows(); ++n)
  {
    for (int k = 0; k < ref.cols(); ++k)
    {
      ASSERT_NEAR(ref(n,k), act(n,k), tol);
    }
  }
}

TEST(MatFunc, All)
{
  // Test of matfunc routines
  RNG_reset(0);
  vec act_v, ref_v;
  mat act_m, ref_m;

  vec a = randn(5);
  ASSERT_NEAR(-1.2785, sum(a), tol);
  act_v = cumsum(a);
  ref_v = "-0.538679 -1.81359 -0.78536 -0.913054 -1.2785";
  assert_vec(ref_v, act_v);
  ASSERT_NEAR(0.0329534, prod(a), tol);
  ASSERT_NEAR(3.1227, sum_sqr(a), tol);

  mat A = randn(5, 5);
  act_v = sum(A);
  ref_v = "-0.154075 -5.06832 -0.386055 1.2232 2.46118";
  assert_vec(ref_v, act_v);
  act_v = sum(A, 1);
  ref_v = "-0.154075 -5.06832 -0.386055 1.2232 2.46118";
  assert_vec(ref_v, act_v);
  act_v = sum(A, 2);
  ref_v = "-2.77077 -2.02482 0.621612 0.255076 1.99483";
  assert_vec(ref_v, act_v);

  act_m = cumsum(A);
  ref_m = "-0.878434 -0.531968 -1.27628 0.382527 -0.466618;"
          "0.411677 -4.26467 -0.888925 -0.227868 0.174195;"
          "-0.263705 -4.37742 -1.27645 0.121915 1.62168;"
          "-0.964984 -5.31757 -1.26185 1.41868 2.20682;"
          "-0.154075 -5.06832 -0.386055 1.2232 2.46118";
  assert_mat(ref_m, act_m);
  act_m = cumsum(A, 1);
  ref_m = "-0.878434 -0.531968 -1.27628 0.382527 -0.466618;"
          "0.411677 -4.26467 -0.888925 -0.227868 0.174195;"
          "-0.263705 -4.37742 -1.27645 0.121915 1.62168;"
          "-0.964984 -5.31757 -1.26185 1.41868 2.20682;"
          "-0.154075 -5.06832 -0.386055 1.2232 2.46118";
  assert_mat(ref_m, act_m);
  act_m = cumsum(A, 2);
  ref_m = "-0.878434 -1.4104 -2.68668 -2.30415 -2.77077;"
          "1.29011 -2.44259 -2.05523 -2.66563 -2.02482;"
          "-0.675383 -0.788138 -1.17566 -0.825877 0.621612;"
          "-0.701279 -1.64143 -1.62682 -0.330059 0.255076;"
          "0.810909 1.06016 1.93595 1.74047 1.99483";
  assert_mat(ref_m, act_m);

  act_v = prod(A);
  ref_v = "-0.435261 0.0524651 0.00244989 0.0207031 -0.0644198";
  assert_vec(ref_v, act_v);
  act_v = prod(A, 1);
  ref_v = "-0.435261 0.0524651 0.00244989 0.0207031 -0.0644198";
  assert_vec(ref_v, act_v);
  act_v = prod(A, 2);
  ref_v = "0.106454 0.729624 -0.0149416 0.0073047 -0.00880157";
  assert_vec(ref_v, act_v);

  act_v = sum_sqr(A);
  ref_v = "4.04154 15.1747 2.69632 2.36107 3.13068";
  assert_vec(ref_v, act_v);
  act_v = sum_sqr(A, 1);
  ref_v = "4.04154 15.1747 2.69632 2.36107 3.13068";
  assert_vec(ref_v, act_v);
  act_v = sum_sqr(A, 2);
  ref_v = "3.04758 16.5307 2.8366 3.39987 1.58962";
  assert_vec(ref_v, act_v);

  act_m = repmat(a, 1, 3);
  ref_m = "-0.538679 -0.538679 -0.538679;"
          "-1.27491 -1.27491 -1.27491;"
          "1.02823 1.02823 1.02823;"
          "-0.127694 -0.127694 -0.127694;"
          "-0.36545 -0.36545 -0.36545;";
  assert_mat(ref_m, act_m);
  act_m = repmat(a, 3, 1, true);
  ref_m = "-0.538679 -1.27491 1.02823 -0.127694 -0.36545;"
          "-0.538679 -1.27491 1.02823 -0.127694 -0.36545;"
          "-0.538679 -1.27491 1.02823 -0.127694 -0.36545";
  assert_mat(ref_m, act_m);
  act_m = repmat(A, 2, 2);
  ref_m = "-0.878434 -0.531968 -1.27628 0.382527 -0.466618 -0.878434 -0.531968 -1.27628 0.382527 -0.466618;"
          "1.29011 -3.7327 0.387353 -0.610395 0.640813 1.29011 -3.7327 0.387353 -0.610395 0.640813;"
          "-0.675383 -0.112755 -0.387522 0.349783 1.44749 -0.675383 -0.112755 -0.387522 0.349783 1.44749;"
          "-0.701279 -0.940147 0.0146015 1.29677 0.585136 -0.701279 -0.940147 0.0146015 1.29677 0.585136;"
          "0.810909 0.249247 0.875791 -0.19548 0.254363 0.810909 0.249247 0.875791 -0.19548 0.254363;"
          "-0.878434 -0.531968 -1.27628 0.382527 -0.466618 -0.878434 -0.531968 -1.27628 0.382527 -0.466618;"
          "1.29011 -3.7327 0.387353 -0.610395 0.640813 1.29011 -3.7327 0.387353 -0.610395 0.640813;"
          "-0.675383 -0.112755 -0.387522 0.349783 1.44749 -0.675383 -0.112755 -0.387522 0.349783 1.44749;"
          "-0.701279 -0.940147 0.0146015 1.29677 0.585136 -0.701279 -0.940147 0.0146015 1.29677 0.585136;"
          "0.810909 0.249247 0.875791 -0.19548 0.254363 0.810909 0.249247 0.875791 -0.19548 0.254363";
  assert_mat(ref_m, act_m);

  // Kronecker test
  mat X = to_mat(randi(2, 2, 1, 4));
  mat Y = randn(3, 3);
  act_m = kron(X, Y);
  ref_m = "-2.93292 1.78242 -2.87218 -2.93292 1.78242 -2.87218;"
          "-0.301315 -1.01854 -0.296556 -0.301315 -1.01854 -0.296556;"
          "-1.28323 -3.50521 0.336507 -1.28323 -3.50521 0.336507;"
          "-2.93292 1.78242 -2.87218 -4.39937 2.67363 -4.30826;"
          "-0.301315 -1.01854 -0.296556 -0.451972 -1.52781 -0.444834;"
          "-1.28323 -3.50521 0.336507 -1.92484 -5.25782 0.504761";
  assert_mat(ref_m, act_m);

  // sqrtm of a real matrix
  A = randn(3, 3);
  cmat A_sqrtm = sqrtm(A);
  ASSERT_NEAR(0, norm(A_sqrtm * A_sqrtm - to_cmat(A)), tol);

  // sqrtm of a complex matrix
  cmat B = randn_c(3, 3);
  cmat B_sqrtm = sqrtm(B);
  ASSERT_NEAR(0, norm(B_sqrtm * B_sqrtm - B), tol);

  // Rank test
  A = randn(3, 3);
  ASSERT_EQ(3, itpp::rank(A));
  A.set_row(1, 3.0 * A.get_row(0));
  ASSERT_EQ(2, itpp::rank(A));
  B = randn_c(3, 3);
  ASSERT_EQ(3, itpp::rank(B));
  B.set_col(1, B.get_col(0));
  ASSERT_EQ(2, itpp::rank(B));
}
