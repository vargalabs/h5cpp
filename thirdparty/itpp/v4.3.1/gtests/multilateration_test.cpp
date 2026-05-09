/*!
 * \file
 * \brief Unit tests for multilateration class
 * \author Bogdan Cristea
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

#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

#ifdef _MSC_VER
#include <float.h>
#define isnan _isnan
#endif

#define BASE_LINE_M 10.0
#define SPHERE_RADIUS_M 450

//uncomment the line below in order to run full tests
//this is needed in order to make sure that the tests pass on non-i386 architectures
//#define FULL_TESTS

static
bool point_eq(const vec &expect, const vec &actual, double eps)
{
  if((3 != length(expect)) || (3 != length(actual))) {
    it_warning("invalid input");
    return false;
  }
  bool out = ((fabs(expect(0) - actual(0)) < eps) && (fabs(expect(1) - actual(1)) < eps) && (fabs(expect(2) - actual(2)) < eps)
              && !isnan(actual(0)) && !isnan(actual(1)) && !isnan(actual(2)));
  if(false == out) {
    cout << "expect " << expect << endl;
    cout << "actual " << actual << endl;
  }
  return out;
}

static
double get_dist(const vec &p0, const vec &p1)
{
  double out = 0.0;
  for(int n = 0; n < 3; ++n) {
    out += (p1[n] - p0[n]) * (p1[n] - p0[n]);
  }
  return sqrt(out);
}

//used for both spheric and hybrid multilateration
static
bool generate_meas(vec &meas, const bvec &method, const mat &bs_pos, const vec &ms_pos)
{
  unsigned int method_len = length(method);
  unsigned int nb_bs = bs_pos.cols();
  bool column = true;
  if(3 == nb_bs) {
    nb_bs = bs_pos.rows();
    column = false;
  }
  if((nb_bs < method_len) || (3 != length(ms_pos))) {
    return false;
  }
  meas.set_size(method_len);
  vec pos(3);
  vec pos_ref(3);
  pos_ref = column ? bs_pos.get_col(0) : bs_pos.get_row(0);
  for(unsigned int k = 0; k < method_len; ++k) {
    if(bin(1) == method[k]) {  /* hyperbolic */
      pos = column ? bs_pos.get_col(k + 1) : bs_pos.get_row(k + 1);
      meas[k] = get_dist(pos, ms_pos) - get_dist(pos_ref, ms_pos);
    }
    else { /* spherical */
      pos = column ? bs_pos.get_col(k) : bs_pos.get_row(k);
      meas[k] = get_dist(pos, ms_pos);
    }
  }
  return true;
}

//used only for hyperbolic multilateration
static
bool generate_meas(mat &meas, const bvec &method, const mat &bs_pos, const vec &ms_pos)
{
  unsigned int k;
  unsigned int i;
  unsigned int method_len = length(method);
  unsigned int nb_bs = bs_pos.cols();
  bool column = true;
  if(3 == nb_bs) {
    nb_bs = bs_pos.rows();
    column = false;
  }
  if((nb_bs < method_len) || (3 != length(ms_pos))) {
    return false;
  }
  meas.set_size(nb_bs, nb_bs);
  vec pos_i(3);
  vec pos_k(3);
  for(k = 0; k < nb_bs; ++k) {
    pos_k = column ? bs_pos.get_col(k) : bs_pos.get_row(k);
    for(i = 0; i < nb_bs; ++i) {
      pos_i = column ? bs_pos.get_col(i) : bs_pos.get_row(i);
      meas(i, k) = get_dist(pos_i, ms_pos) - get_dist(pos_k, ms_pos);
    }
  }
  return true;
}

static
bool get_bs(mat &bs_pos, unsigned int nb_bs, double l)
{
  int n, i, j, k;
  static const double SQRT2M1 = 0.7071067811865475;

  if(8 < nb_bs) {
    it_warning("at most 8 BSs");
    return false;
  }
  n = 0;
  i = j = k = 0;
  bs_pos.set_size(3, nb_bs);
  while((unsigned int)n < nb_bs) {
    bs_pos(0, n) = (0.5 - i) * l;
    bs_pos(1, n) = (0.5 - j) * l;
    bs_pos(2, n++) = (1 - 2 * k) * l * SQRT2M1;
    switch(n) { /*ensure that the first 4 BSs are not on the same plane*/
    case 3:
      i = j = 0;
      k = 1;
      break;
    case 5:
      i = k = 1;
      j = 0;
      break;
    case 4:
      k = 0;
      i = j = 1;
      break;
    default:
      i = (i + 1) % 2;
      if((0 != n) && (0 == n % 2)) {
        j = (j + 1) % 2;
      }
    }
  }
  return true;
}

static
vec get_ms(double radius)
{
  static const double SQRT3M1 = 0.5773502691896258;
  vec out;
  out.set_size(3);
  for(int n = 0; n < 3; ++n) {
    out[n] = SQRT3M1 * radius * (2.0 * randu() - 1.0);
  }
  return out;
}

TEST(Multilateration, get_pos)
{
  RNG_reset(0);

  const unsigned int nb_points = 10;
  const double eps = 1e-3;
  mat bs_pos;
  unsigned int method_len = 0;//NB: the number of BSs is greater by one for hyprid and hyperbolic
  bvec method;
  vec ms_pos;
  vec actual_ms_pos;
  unsigned int i;
  Multilateration multi;
  vec meas;//measurements vector for spherical and hybrid multilateration
  mat meas_hyper;//measurements matrix for hyperbolic multilateration

  //test inputs
  bs_pos.set_size(3, 4);
  method.set_size(4);
  method.zeros();
  multi.setup(method, bs_pos);
  ASSERT_TRUE(Multilateration::MULTI_SPHERICAL == multi.get_type());
  method(0) = 1;
  bs_pos.set_size(3, 5);//nb of BSs is the number of measures plus one
  multi.setup(method, bs_pos);
  ASSERT_TRUE(Multilateration::MULTI_HYBRID == multi.get_type());
  method.ones();
  multi.setup(method, bs_pos);
  ASSERT_TRUE(Multilateration::MULTI_HYPERBOLIC == multi.get_type());

  for(method_len = 4; method_len < 8; ++method_len) {
    //spherical multilateration
    ASSERT_TRUE(get_bs(bs_pos, method_len, BASE_LINE_M));
    method.set_size(method_len);
    method.zeros();
    multi.setup(method, bs_pos);
    ms_pos = get_ms(SPHERE_RADIUS_M);
    ASSERT_TRUE(generate_meas(meas, method, bs_pos, ms_pos));
    ASSERT_TRUE(multi.get_pos(actual_ms_pos, meas));
    ASSERT_TRUE(point_eq(ms_pos, actual_ms_pos, eps));
    actual_ms_pos.zeros();

    //hybrid multilateration
    ASSERT_TRUE(get_bs(bs_pos, method_len + 1, BASE_LINE_M));
    for(i = 0; i < (method_len - 1); ++i) {
      method[i] = 1;
      multi.setup(method, bs_pos);
      ms_pos = get_ms(SPHERE_RADIUS_M);
      ASSERT_TRUE(generate_meas(meas, method, bs_pos, ms_pos));
      ASSERT_TRUE(multi.get_pos(actual_ms_pos, meas));
      ASSERT_TRUE(point_eq(ms_pos, actual_ms_pos, eps));
      actual_ms_pos.zeros();
    }

    //hyperbolic multilateration
    method[i] = 1;
    multi.setup(method, bs_pos);
    ms_pos = get_ms(SPHERE_RADIUS_M);
    ASSERT_TRUE(generate_meas(meas_hyper, method, bs_pos, ms_pos));
    ASSERT_TRUE(multi.get_pos(actual_ms_pos, meas_hyper));
    ASSERT_TRUE(point_eq(ms_pos, actual_ms_pos, eps));
    actual_ms_pos.zeros();
  }
#ifdef FULL_TESTS
  //test case when the last measure is always from TDOA
  for(method_len = 5; method_len < 8; ++method_len) {
    ASSERT_TRUE(get_bs(bs_pos, method_len + 1, BASE_LINE_M));
    method.set_size(method_len);
    method.ones();
    for(i = 0; i < (method_len - 1); ++i) {
      method[i] = 0;
      multi.setup(method, bs_pos);
      ms_pos = get_ms(SPHERE_RADIUS_M);
      EXPECT_TRUE(generate_meas(meas, method, bs_pos, ms_pos));
      EXPECT_TRUE(multi.get_pos(actual_ms_pos, meas));
      EXPECT_TRUE(point_eq(ms_pos, actual_ms_pos, eps));
      actual_ms_pos.zeros();
    }
  }
#endif
}

TEST(Multilateration, get_crlb)
{
  vec ms_pos(3);
  mat bs_pos;
  bvec method;
  unsigned int method_len = 4;
  double sigma2 = 0.0;
  Multilateration multi;
  unsigned int i;

  method.set_size(method_len);
  method.zeros();
  ASSERT_TRUE(get_bs(bs_pos, method_len, BASE_LINE_M));
  multi.setup(method, bs_pos);

  ms_pos = get_ms(SPHERE_RADIUS_M);

  double crlb = multi.get_crlb(ms_pos, sigma2);
  ASSERT_EQ(0.0, crlb);

  sigma2 = 1e-6;
  for(method_len = 4; method_len < 8; ++method_len) {
    ASSERT_TRUE(get_bs(bs_pos, method_len + 1, BASE_LINE_M));
    method.set_size(method_len);
    method.zeros();
    for(i = 0; i < method_len; ++i) {
      method(i) = 1;
      multi.setup(method, bs_pos);
      ms_pos = get_ms(SPHERE_RADIUS_M);
      crlb = multi.get_crlb(ms_pos, sigma2);
      ASSERT_NEAR(0.0, crlb, (i < 3) ? 1e-1 : 12);
    }
  }

  ms_pos.ones();
  sigma2 = 0.1;
  for(method_len = 4; method_len < 8; ++method_len) {
    ASSERT_TRUE(get_bs(bs_pos, method_len, BASE_LINE_M));
    method.set_size(method_len);
    method.zeros();
    multi.setup(method, bs_pos);
    crlb = multi.get_crlb(ms_pos, sigma2);
    ASSERT_NEAR(0.5, crlb, 0.2);
    ASSERT_TRUE(get_bs(bs_pos, method_len + 1, BASE_LINE_M));
    for(i = 0; i < method_len; ++i) {
      method(i) = 1;
      multi.setup(method, bs_pos);
      crlb = multi.get_crlb(ms_pos, sigma2);
      ASSERT_NEAR(0.5, crlb, 0.6);
    }
  }
}
