/*!
 * \file
 * \brief Galois Field algebra classes test program
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

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST (Galois, All)
{
  GF a(8), b(8), c(8);

  a = 4;
  b = 2;

  c = a + b;

  ostringstream ss(ostringstream::out);
  string ref[] = {"a=alpha^4, b=alpha^2", "c=alpha^1"};

  ss << "a=" << a << ", b=" << b;
  ASSERT_TRUE(ss.str() == ref[0]);
  ss.str("");
  ss << "c=" << c;
  ASSERT_TRUE(ss.str() == ref[1]);
}

