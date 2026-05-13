/*!
 * \file
 * \brief Timer classes test program
 * \author Thomas Eriksson, Tony Ottosson, Tobias Ringstrom and Adam Piatyszek
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

TEST(Timer, All)
{
  const double period = 2.0;
  const double relative_error = 0.05;
  CPU_Timer t1;
  Real_Timer t2;

  t1.start();
  while (t1.get_time() < period) ;
  t1.stop();

  tic();
  t2.start();
  while (t2.get_time() < period) ;
  t2.stop();
  double t3 = toc();

  ASSERT_TRUE(fabs(t1.get_time() - period) <= relative_error * period);

  ASSERT_TRUE(fabs(t2.get_time() - period) <= relative_error * period);

  ASSERT_TRUE(t3 >= t2.get_time());
}
