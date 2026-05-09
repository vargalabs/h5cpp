/*!
 * \file
 * \brief Fast math test program
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

static
void assert_mat_p(const mat &ref, const mat &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.rows(), act.rows()) << line;
  ASSERT_EQ(ref.cols(), ref.cols()) << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol) << line;
    }
  }
}
#define assert_mat(ref, act) assert_mat_p(ref, act, __LINE__)

TEST(FastMath, All)
{
  // Test of fastmath

  mat m0("1 2 3;4 5 6;7 8 9"), mv0("2;3;1");
  vec v0("2 3 1");

  // sub_v_vT_m: the slow and fast way
  mat act = m0 - mv0*transpose(mv0)*m0;
  mat ref = "-41 -52 -63; -59 -76 -93; -14 -19 -24";
  assert_mat(ref, act);
  sub_v_vT_m(m0, v0);
  ref = "-41 -52 -63; -59 -76 -93; -14 -19 -24";
  assert_mat(ref, m0);

  m0 = "1 2 3;4 5 6;7 8 9";

  // sub_m_v_vT: the slow and fast way
  act = m0 - m0*mv0*transpose(mv0);
  ref = "-21 -31 -8; -54 -82 -23; -87 -133 -38";
  assert_mat(ref, act);
  sub_m_v_vT(m0, v0);
  ref = "-21 -31 -8; -54 -82 -23; -87 -133 -38";
  assert_mat(ref, m0);
}
