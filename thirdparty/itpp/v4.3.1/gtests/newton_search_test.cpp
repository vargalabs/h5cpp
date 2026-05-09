/*!
 * \file
 * \brief Newton search test program
 * \author Tony Ottosson
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

#include <itpp/itoptim.h>
#include "gtest/gtest.h"

using namespace itpp;

double rosenbrock(const vec &x)
{
  double f1 = x(1) - sqr(x(0)), f2 = 1 - x(0);
  double F = 50 * sqr(f1) + 0.5 * sqr(f2) + 0.0;

  return F;
}

vec rosenbrock_gradient(const vec &x)
{
  double f1 = x(1) - sqr(x(0)), f2 = 1 - x(0);
  vec g(2);
  g(0) = -200.0 * x(0) * f1 - f2;
  g(1) = 100.0 * f1;

  return g;
}

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

TEST(NewtonSearch, All)
{

  // Test of Numerical optimization

  {
    // Line_Search
    vec x0 = "-1.2 1";
    double F0 = rosenbrock(x0);
    vec g0 = rosenbrock_gradient(x0);
    vec h = -g0;

    Line_Search line;
    line.set_functions(rosenbrock, rosenbrock_gradient);

    ASSERT_NEAR(12.1, F0, tol);
    vec ref = "-107.8 -44";
    assert_vec(ref, g0);
    ref = "107.8 44";
    assert_vec(ref, h);

    double Fn;
    vec xn, gn;

    line.enable_trace();

    ASSERT_TRUE(line.search(x0, F0, g0, h, xn, Fn, gn));

    // Soft
    ref = "-0.876168 1.13218";
    assert_vec(ref, xn);
    ASSERT_NEAR(8.40327, Fn, tol);
    ref = "61.9977 36.4507";
    assert_vec(ref, gn);
    ASSERT_EQ(4, line.get_no_function_evaluations());

    vec alpha_values, F_values, dF_values;
    line.get_trace(alpha_values, F_values, dF_values);

    // trace
    ref = "0 1 0.1 0.01 0.00300401";
    assert_vec(ref, alpha_values);
    ref = "12.1 6405495599.3599997 373080.93204800034 102.17722267280001 8.40327";
    assert_vec(ref, F_values);
    ref = "-13556.84 25963610893.439999 17461465.02672001 9898.07 8287.1771139224893";
    assert_vec(ref, dF_values);

    line.set_method(Exact);
    ASSERT_TRUE(line.search(x0, F0, g0, h, xn, Fn, gn));

    // Exact
    ref = "-1.03008 1.06935";
    assert_vec(ref, xn);
    ASSERT_NEAR(2.06405, Fn, tol);
    ref = "-0.322427 0.828893";
    assert_vec(ref, gn);
    ASSERT_EQ(7, line.get_no_function_evaluations());

    line.get_trace(alpha_values, F_values, dF_values);

    // trace
    ref = "0 1 0.1 0.01 0.00300401 0.00165196 0.00148676 0.00157625";
    assert_vec(ref, alpha_values);
    ref = "12.1 6405495599.3599997 373080.93204800034 102.17722267280001 8.40327 2.08432 2.09243 2.06405";
    assert_vec(ref, F_values);
    ref = "-13556.84 25963610893.439999 17461465.02672001 9898.07 8287.1771139224893 531.883 -638.60641049642675 1.71365";
    assert_vec(ref, dF_values);
  }

  {
    // Newton_Search
    vec x0 = "-1.2 1";

    Newton_Search newton;
    newton.set_functions(rosenbrock, rosenbrock_gradient);
    newton.enable_trace();

    vec xn, gn;

    ASSERT_TRUE(newton.search(x0, xn));

    vec ref = "1 1";
    assert_vec(ref, xn);
    ASSERT_NEAR(5.24508e-14, newton.get_function_value(), tol);
    ASSERT_NEAR(1.41687e-06, newton.get_stop_1(), tol);
    ASSERT_NEAR(9.32294e-06, newton.get_stop_2(), tol);
    ASSERT_EQ(38, newton.get_no_function_evaluations());
    ASSERT_EQ(34, newton.get_no_iterations());

    Array<vec> xv;
    vec Fv, ngv, dv;
    newton.get_trace(xv, Fv, ngv, dv);

    // trace
    ref = "-0.916319 1.11579";
    assert_vec(ref, xv(0));
    ref = "1 1";
    assert_vec(ref, xv(xv.length()-1));
    ref = "5.64904 1.5618 1.41629 1.2203 1.20012 1.18688 1.16266 1.11426 1.03385 0.841214 0.682837 0.634992 "
    "0.528517 0.395891 0.349639 0.224737 0.202608 0.149765 0.109608 0.0815102 0.0465651 0.0379422 0.0229203 "
    "0.0131056 0.00545953 0.0026574 0.000729924 0.000189087 2.7111e-05 4.66978e-06 5.98346e-07 4.12168e-08 "
    "2.55577e-11 5.24508e-14";
    assert_vec(ref, Fv);
    ref = "48.6917 6.83035 9.16562 1.45963 2.87505 3.70968 4.45032 5.01409 4.59566 3.08462 1.58691 2.36257 "
    "3.0891 1.40777 3.52707 0.767488 3.04116 0.907026 1.22519 1.64112 0.351909 2.13814 0.643864 1.16604 "
    "0.458514 0.87564 0.575896 0.384089 0.0473799 0.0132776 0.0183764 0.00561985 0.000112443 1.41687e-06";
    assert_vec(ref, ngv);
    ref = "0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 "
    "0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 0.1225 "
    "0.1225 0.1225";
    assert_vec(ref, dv);

    newton.disable_trace();

    // fminunc
    xn = fminunc(rosenbrock, rosenbrock_gradient, x0);
    ref = "-1.2 1";
    assert_vec(ref, x0);
    ref = "1 1";
    assert_vec(ref, xn);
  }
}
