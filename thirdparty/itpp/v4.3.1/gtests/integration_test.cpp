/*!
 * \file
 * \brief Transforms test program
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood, Adam Piatyszek, Andy Panov and Bogdan Cristea
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

//Integration is tested with x*ln(x) function
double test_integrand(const double x)
{
  return x * log(x);
}

struct Test_Integrand {
  double operator()(double x) {
    return x * log(x);
  }
};

//set test tolerance (measure of relative and absolute error)
const double max_rel_error = 1e-6;
const double max_abs_error = 1e-6;

//results tester
inline void test_result(double q, double ref)
{
  if(abs(q - ref) < max_abs_error) return; //handle numbers with absolute value close to zero (relative error can be huge for them)
  double rel_error = abs(q - ref) / abs(q);
  ASSERT_LE(rel_error, max_rel_error);
}

//analytically computed integral
double integral_value(double a, double b)
{
  return 0.25 * b * b * (2 * log(b) - 1.0) - 0.25 * a * a * (2 * log(a) - 1.0);
}
//----------------------------------------------
//Gtest test cases
//----------------------------------------------
static const double low_lim = 1.5;
static const double hi_lim = 3.5;
static const double ref_integral_value = integral_value(low_lim, hi_lim);

TEST(Integration, Simpson)
{
  {
    SCOPED_TRACE("Testing function object");
    double q = itpp::quad(Test_Integrand(), low_lim, hi_lim);
    test_result(q, ref_integral_value);
  }
  {
    SCOPED_TRACE("Testing double f(double)");
    double q = itpp::quad(test_integrand, low_lim, hi_lim);
    test_result(q, ref_integral_value);
  }
}

TEST(Integration, Lobatto)
{
  {
    SCOPED_TRACE("Testing function object");
    double q = itpp::quadl(Test_Integrand(), low_lim, hi_lim);
    test_result(q, ref_integral_value);
  }
  {
    SCOPED_TRACE("Testing double f(double)");
    double q = itpp::quadl(test_integrand, low_lim, hi_lim);
    test_result(q, ref_integral_value);
  }
}



