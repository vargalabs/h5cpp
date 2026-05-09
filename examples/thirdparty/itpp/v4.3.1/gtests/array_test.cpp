/*!
 * \file
 * \brief Array class test program
 * \author Adam Piatyszek
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

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace std;
using namespace itpp;

TEST(Array, Operations)
{
  ostringstream ss(ostringstream::out);//replaces cout
  const string ref[] = {"Testing simple Arrays of integers:\
A1 = {1 1 42 1 1 1 1 1 1 1}\
A2 = {5 5 5 5 5 5 5 5 5 5 5 5 5 5 5}",\
"Testing Array initialisation with: \"{[1 1] [1 0 1 0] [0 0 1]}\":\
A3 = {[1 1] [1 0 1 0] [0 0 1]}",\
"Testing Array initialisation with: \"{{[5 3; 020 4] [1 0; 3 9]} {[0 -3; 1 0xa]}}\":\
A4 = {{[[5 3]\n\
 [16 4]] [[1 0]\n\
 [3 9]]} {[[0 -3]\n\
 [1 10]]}}",\
"Testing Array::operator():\
A4(1) = {[[0 -3]\n\
 [1 10]]}\
A4(0)(1) = [[1 0]\n\
 [3 9]]",\
"Testing Array::left(), Array::right() and Array::mid():\
A1.left(4) = {1 1 42 1}\
A1.right(5) = {1 1 1 1 1}\
A1.mid(2, 3) = {42 1 1}",\
"Testing A4.swap(0, 1):\
{{[[0 -3]\n\
 [1 10]]} {[[5 3]\n\
 [16 4]] [[1 0]\n\
 [3 9]]}}"};

  Array<int> A1(10), A2(15);

  A1 = 1;
  A1(2) = 42;
  A2 = 5;

  ss << "Testing simple Arrays of integers:" << "A1 = " << A1 << "A2 = " << A2;
  ASSERT_TRUE(ss.str() == ref[0]);
  ss.str("");

  // Test of Array initialisation by string:
  Array<bvec> A3 = "{[1 1] [1 0 1 0] [0 0 1]}";
  ss << "Testing Array initialisation with: \"{[1 1] [1 0 1 0] [0 0 1]}\":"
       << "A3 = " << A3;
  ASSERT_TRUE(ss.str() == ref[1]);
  ss.str("");

  Array<Array<imat> > A4 = "{{[5 3; 020 4] [1 0; 3 9]} {[0 -3; 1 0xa]}}";
  ss << "Testing Array initialisation with: \"{{[5 3; 020 4] [1 0; 3 9]} {[0 -3; 1 0xa]}}\":"
       << "A4 = " << A4;
  ASSERT_TRUE(ss.str() == ref[2]);
  ss.str("");

  // Test of operator()
  ss << "Testing Array::operator():"
       << "A4(1) = " << A4(1) << "A4(0)(1) = " << A4(0)(1);
  ASSERT_TRUE(ss.str() == ref[3]);
  ss.str("");

  // Test of left(), right() and mid() methods:
  ss << "Testing Array::left(), Array::right() and Array::mid():"
       << "A1.left(4) = " << A1.left(4) << "A1.right(5) = " << A1.right(5)
       << "A1.mid(2, 3) = " << A1.mid(2, 3);
  ASSERT_TRUE(ss.str() == ref[4]);
  ss.str("");

  // Test of swap function
  A4.swap(0, 1);
  ss << "Testing A4.swap(0, 1):" << A4;
  ASSERT_TRUE(ss.str() == ref[5]);
}
