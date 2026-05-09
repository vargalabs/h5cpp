/*!
 * \file
 * \brief Deterministic sources test program
 * \author Tobias Ringstrom, Tony Ottosson and Adam Piatyszek
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

#include <itpp/itsignal.h>
#include <itpp/itstat.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

#define LOOP_SIZE 100000

static
const double tol = 1e-4;
static
double m;
static
double var;

#define REALRUN(name,s)     \
  for (int i = 0; i < LOOP_SIZE; i++)   \
    real_result(i) = s();    \
  m = mean(real_result), var = variance(real_result);

TEST(Source, All)
{
  Sine_Source      s10(20.0 / LOOP_SIZE);
  Square_Source    s11(20.0 / LOOP_SIZE);
  Triangle_Source  s12(20.0 / LOOP_SIZE);
  Sawtooth_Source  s13(20.0 / LOOP_SIZE);
  Impulse_Source   s14(20.0 / LOOP_SIZE);
  Pattern_Source   s15(vec("1 3"));

  RNG_reset(12345);

  vec real_result(LOOP_SIZE);
  REALRUN("Sine",        s10);
  ASSERT_NEAR(5.38742e-14, m, tol);
  ASSERT_NEAR(0.500005, var, tol);
  REALRUN("Square",      s11);
  ASSERT_NEAR(2e-05, m, tol);
  ASSERT_NEAR(1.00001, var, tol);
  REALRUN("Triangle",    s12);
  ASSERT_NEAR(-1.52056e-17, m, tol);
  ASSERT_NEAR(0.333337, var, tol);
  REALRUN("Sawtooth",    s13);
  ASSERT_NEAR(0.0002, m, tol);
  ASSERT_NEAR(0.333337, var, tol);
  REALRUN("Impulse",     s14);
  ASSERT_NEAR(0.00019, m, tol);
  ASSERT_NEAR(0.000189966, var, tol);
  REALRUN("Pattern",     s15);
  ASSERT_NEAR(2, m, tol);
  ASSERT_NEAR(1.00001, var, tol);
}
