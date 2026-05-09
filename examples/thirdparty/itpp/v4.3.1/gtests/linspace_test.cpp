/*!
 * \file
 * \brief linspace functions test program
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

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace itpp;

static
void assert_vec(const vec &ref, const vec &act)
{
  static const double TOL = 1e-12;
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n)
  {
    ASSERT_NEAR(ref(n), act(n), TOL);
  }
}

TEST(Linspace, All)
{
    double from = 0;
    double to = 2;
    int points = 1;
    vec v = linspace(from, to, points);
    ASSERT_EQ(points, v.length());
    ASSERT_EQ(to, v(0));

    from = 0;
    to = 10;
    points = 30;
    v = linspace(from, to, points);
    ASSERT_EQ(points, v.length());
    vec actual_v(v.length());
    int n;
    for (n = 0; n < v.length()-1; ++n)
    {
        actual_v(n) = n*(to-from)/(points-1);
    }
    actual_v(n) = to;
    assert_vec(actual_v, v);

    from = 0;
    to = 10;
    double step = 1.0;
    v = linspace_fixed_step(from, to, step);
    int actual_len = floor_i((to-from)/step)+1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    assert_vec(actual_v, v);

    from = 0;
    to = 10;
    step = 1.1;
    v = linspace_fixed_step(from, to, step);
    actual_len = floor_i((to-from)/step)+1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    assert_vec(actual_v, v);

    from = 0;
    to = 1;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    actual_len = 1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    assert_vec(actual_v, v);

    from = 0;
    to = -1;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    actual_len = 1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    assert_vec(actual_v, v);

    from = 0;
    to = -1;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    actual_len = 0;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    ASSERT_TRUE(v == actual_v);

    from = -1;
    to = 0;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    actual_len = 1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    actual_v(0) = from;
    ASSERT_TRUE(v == actual_v);

    from = -1;
    to = 0;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    actual_len = 0;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    ASSERT_TRUE(v == actual_v);

    from = -5;
    to = 5;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    actual_len = floor_i((to-from)/step)+1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    ASSERT_TRUE(v == actual_v);

    from = 5;
    to = -5;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    actual_len = floor_i((to-from)/step)+1;
    ASSERT_EQ(actual_len, v.length());
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_len; ++n)
    {
        actual_v(n) = from+n*step;
    }
    ASSERT_TRUE(v == actual_v);

    int ifrom = 0;
    int ito = 10;
    ivec iv = linspace_fixed_step(ifrom, ito);
    actual_len = floor_i(ito-ifrom)+1;
    ASSERT_EQ(actual_len, iv.length());
    ivec actual_iv(actual_len);
    for (n = 0; n < actual_iv.length(); ++n)
    {
        actual_iv(n) = ifrom+n;
    }
    ASSERT_TRUE(iv == actual_iv);

    short int sfrom = 0;
    short int sto = 10;
    svec sv = linspace_fixed_step(sfrom, sto);
    actual_len = floor_i(sto-sfrom)+1;
    ASSERT_EQ(actual_len, sv.length());
    svec actual_sv(actual_len);
    for (n = 0; n < actual_sv.length(); ++n)
    {
        actual_sv(n) = sfrom+n;
    }
    ASSERT_TRUE(sv == actual_sv);
}
