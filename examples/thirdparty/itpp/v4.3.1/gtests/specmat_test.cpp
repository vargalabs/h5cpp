/*!
 * \file
 * \brief Test program of special vectors and matrices
 * \author Adam Piatyszek
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

using namespace std;
using namespace itpp;

static const double tol = 1e-4;
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
static
void assert_cmat_p(const cmat &ref, const cmat &act, int line)
{
  ASSERT_EQ(ref.rows(), act.rows()) << line;
  ASSERT_EQ(ref.cols(), act.cols()) << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k).real(), act(n,k).real(), tol) << line;
      ASSERT_NEAR(ref(n,k).imag(), act(n,k).imag(), tol) << line;
    }
  }
}
#define assert_cmat(ref, act) assert_cmat_p(ref, act, __LINE__)

TEST(SpecMat, All)
{
  RNG_reset(0);

  // Test of specmat routines

  bvec b1 = randb(7);
  bvec b2 = randb(13);
  bmat ref_b = "1 0 1 0 0 0 1 1 1 0 1 1 0;"
               "1 1 0 1 0 0 0 1 1 1 0 1 1;"
               "1 1 1 0 1 0 0 0 1 1 1 0 1;"
               "1 1 1 1 0 1 0 0 0 1 1 1 0;"
               "0 1 1 1 1 0 1 0 0 0 1 1 1;"
               "1 0 1 1 1 1 0 1 0 0 0 1 1;"
               "1 1 0 1 1 1 1 0 1 0 0 0 1";
  ASSERT_TRUE(ref_b == toeplitz(b1, b2));
  ref_b = "1 1 1 1 0 1 1;"
          "1 1 1 1 1 0 1;"
          "1 1 1 1 1 1 0;"
          "1 1 1 1 1 1 1;"
          "0 1 1 1 1 1 1;"
          "1 0 1 1 1 1 1;"
          "1 1 0 1 1 1 1";
  ASSERT_TRUE(ref_b == toeplitz(b1));

  vec v1 = randn(5);
  vec v2 = randn(7);
  mat ref = "0.0146015 -0.19548 -0.466618 0.640813 1.44749 0.585136 0.254363;"
            "0.875791 0.0146015 -0.19548 -0.466618 0.640813 1.44749 0.585136;"
            "0.382527 0.875791 0.0146015 -0.19548 -0.466618 0.640813 1.44749;"
            "-0.610395 0.382527 0.875791 0.0146015 -0.19548 -0.466618 0.640813;"
            "0.349783 -0.610395 0.382527 0.875791 0.0146015 -0.19548 -0.466618";
  assert_mat(ref, toeplitz(v1, v2));
  ref = "0.0146015 0.875791 0.382527 -0.610395 0.349783;"
        "0.875791 0.0146015 0.875791 0.382527 -0.610395;"
        "0.382527 0.875791 0.0146015 0.875791 0.382527;"
        "-0.610395 0.382527 0.875791 0.0146015 0.875791;"
        "0.349783 -0.610395 0.382527 0.875791 0.0146015";
  assert_mat(ref, toeplitz(v1));

  cvec c1 = randn_c(4);
  cvec c2 = randn_c(3);
  cmat ref_c = "0.480134+0.0415308i -1.01547-0.104848i 0.118973+1.59312i;"
               "0.388173+0.90716i 0.480134+0.0415308i -1.01547-0.104848i;"
               "-1.03694-0.106531i 0.388173+0.90716i 0.480134+0.0415308i;"
               "-0.45369+0.630181i -1.03694-0.106531i 0.388173+0.90716i";
  assert_cmat(ref_c, toeplitz(c1, c2));
  ref_c = "0.480134+0.0415308i 0.388173+0.90716i -1.03694-0.106531i -0.45369+0.630181i;"
          "0.388173+0.90716i 0.480134+0.0415308i 0.388173+0.90716i -1.03694-0.106531i;"
          "-1.03694-0.106531i 0.388173+0.90716i 0.480134+0.0415308i 0.388173+0.90716i;"
          "-0.45369+0.630181i -1.03694-0.106531i 0.388173+0.90716i 0.480134+0.0415308i";
  assert_cmat(ref_c, toeplitz(c1, c1));
  ref_c = "0.480134+0.0415308i 0.388173+0.90716i -1.03694-0.106531i -0.45369+0.630181i;"
          "0.388173-0.90716i 0.480134+0.0415308i 0.388173+0.90716i -1.03694-0.106531i;"
          "-1.03694+0.106531i 0.388173-0.90716i 0.480134+0.0415308i 0.388173+0.90716i;"
          "-0.45369-0.630181i -1.03694+0.106531i 0.388173-0.90716i 0.480134+0.0415308i";
  assert_cmat(ref_c, toeplitz(c1));
}
