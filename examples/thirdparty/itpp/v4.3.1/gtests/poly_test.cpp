/*!
 * \file
 * \brief Polynomial routines test program
 * \author Tony Ottosson, Adam Piatyszek and Kumar Appaiah
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

#include <itpp/itsignal.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
const double tol = 1e-4;
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
void assert_cvec_p(const cvec &ref, const cvec &act, int line)
{
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n).real(), act(n).real(), tol) << line;
    ASSERT_NEAR(ref(n).imag(), act(n).imag(), tol) << line;
  }
}
#define assert_cvec(ref, act) assert_cvec_p(ref, act, __LINE__)


TEST(Poly, All)
{
  // Test of polynomial routines
  RNG_reset(0);

  {
    vec r = randn(3);
    vec p = poly(r);
    cvec r2 = roots(p);
    vec ref = "1 0.78536 -1.17803 -0.706159";
    assert_vec(ref, p);
    cvec ref_c = "1.02823+0i -1.27491+0i -0.538679+0i";
    assert_cvec(ref_c, r2);

    r = randn(7);
    p = poly(r);
    r2 = roots(p);
    ref = "1 0.647219 -1.88626 -1.74412 0.41764 0.836933 0.256015 0.0203118";
    assert_vec(ref, p);
    ref_c = "1.29011+0i 0.810909+0i -0.878434+0i -0.701279+0i -0.675383+0i -0.36545+0i -0.127694+0i";
    assert_cvec(ref_c, r2);

    vec x = randn(9);
    vec y = polyval(p, x);
    ref = "0.00138573 -7328.4084433054286 0.00124015 -0.00711672 0.134252 -0.770298 0.217109 0.000511238 0.0242297";
    assert_vec(ref, y);
  }

  {
    // Complex polynomials
    cvec r = randn_c(3);
    cvec p = poly(r);
    cvec r2 = roots(p);
    cvec ref_c = "1+0i -1.10462-0.379596i -0.0905352+0.4853i 0.301401-0.079591i";
    assert_cvec(ref_c, p);
    ref_c = "-0.431614+0.247334i 0.916952-0.138225i 0.619278+0.270487i";
    assert_cvec(ref_c, r2);

    r = randn_c(7);
    p = poly(r);
    r2 = roots(p);
    ref_c = "1+0i -2.34578+0.115557i 2.61793+0.471743i -1.41308-1.24427i 0.50249+0.79891i -0.119019-0.3526i "
    "0.0871739+0.177806i 0.0285549-0.0501647i";
    assert_cvec(ref_c, p);
    ref_c = "0.90716-1.03694i 1.02353+0.413753i 0.630181-0.360109i -0.106531-0.45369i -0.329949+0.453123i "
    "0.179862+0.480134i 0.0415308+0.388173i";
    assert_cvec(ref_c, r2);

    cvec x = randn_c(9);
    cvec y = polyval(p, x);
    ref_c = "62.9076+74.3315i -0.00908609-0.0502987i -12.7646-15.3766i -0.218263-0.631535i -1.23503-9.01664i "
    "-3.51764+7.19977i 0.273831+0.59537i -0.023377-0.00661798i 0.0748359+0.191528i";
    assert_cvec(ref_c, y);
  }

  {
    // Chebyshev polynomial
    vec x = randn(8);
    vec ref = "2061173.0875914963 151.508 0.186741 0.187703 -0.767701 12825.119165891341 1968.4127832291322 923658.76088428847";
    assert_vec(ref, cheb(10, x));
    ref = "-4184915494.8777189 2637.3130287898421 0.48261 -0.481323 -0.864083 -2054031.5518480411 -123506.33117417867 "
    "1255400590.6012175";
    assert_vec(ref, cheb(15, x));
  }
}
