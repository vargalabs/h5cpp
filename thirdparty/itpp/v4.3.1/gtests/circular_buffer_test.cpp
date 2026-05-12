/*!
 * \file
 * \brief Circular buffer test program
 * \author 
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/circular_buffer.h>
#include <itpp/base/random.h>
#include <itpp/base/matfunc.h>
#include "gtest/gtest.h"

using namespace std;
using namespace itpp;

TEST(CircularBuffer, All)
{
//Set a fixed seed to the random number generator
  RNG_reset(12345);

  int  nrof_elements;
  ivec index_vec;
  vec  a = randn(3);
  vec  b = randn(8);
  vec  out_vec;
  double   out;
  Array <double> out_array;

  Circular_Buffer<double> cb1(10);

//Put the elements of a to the buffer
  cb1.put(a);

//Vector output: Peek at the two oldest elements of the buffer (the two first extracted)
  cb1.peek(out_vec, 2);
  ASSERT_TRUE(a(0,1) == out_vec);

//Vector output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_vec);
  ASSERT_TRUE(reverse(a) == out_vec);

//Array output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_array);
  nrof_elements = out_array.size();
  for (int n = 0; n < nrof_elements; ++n) {
    ASSERT_EQ(a[n], out_array(nrof_elements-n-1));
  }

//Put the elements of \c b to the buffer
  cb1.put(b);

//Extract the oldest element of the buffer
  cb1.get(out_vec, 1);
  ASSERT_EQ(a[1], out_vec[0]);

//Extract all elements of the buffer
  cb1.get(out_vec);
  ASSERT_EQ(a[2], out_vec[0]);
  ASSERT_TRUE(b == out_vec(1,-1));

  cb1.get(out_vec);
  ASSERT_EQ(0, out_vec.size());
  cb1.get(out_array, 0);
  ASSERT_EQ(0, out_array.size());

  for (int i = 0;i < a.length();i++) { cb1.put(a(i)); }
  for (int i = 0;i < b.length();i++) { cb1.put(b(i)); }
  ASSERT_EQ(cb1.size(), cb1.nrof_elements());

  nrof_elements = cb1.nrof_elements();
  for (int i = 0;i < nrof_elements;i++) {
    cb1.peek(i, out);
    if (2 > i) {
      ASSERT_EQ(a[i+1], out);
    }
    else {
      ASSERT_EQ(b[i-2], out);
    }
  }

//Test of error handling
//cout << "get() = " << cb1.get() << endl;

  cb1.put(b);
  ASSERT_EQ(10, cb1.size());
  ASSERT_EQ(10, cb1.nrof_elements());

  cb1.peek(out_vec, -1);
  for (int i = 0; i < out_vec.length(); ++i) {
    if (2 > i) {
      ASSERT_EQ(b[6+i], out_vec[i]);
    }
    else {
      ASSERT_EQ(b[i-2], out_vec[i]);
    }
  }

  index_vec = "5 3 7 1";
  for (int i = 0;i < index_vec.size();i++) {
    cb1.peek(index_vec(i), out);
    if (3 > i) {
      ASSERT_EQ(b[index_vec[i]-2], out);
    }
    else {
      ASSERT_EQ(b[7], out);
    }
  }

  cb1.peek(index_vec, out_vec);
  for (int i = 0;i < index_vec.size();i++) {
    if (3 > i) {
      ASSERT_EQ(b[index_vec[i]-2], out_vec[i]);
    }
    else {
      ASSERT_EQ(b[7], out_vec[i]);
    }
  }

  cb1.set_size(15, true);
  ASSERT_EQ(15, cb1.size());
  ASSERT_EQ(10, cb1.nrof_elements());

  cb1.peek(out_vec, -1);
  for (int i = 0; i < out_vec.length(); ++i) {
    if (2 > i) {
      ASSERT_EQ(b[6+i], out_vec[i]);
    }
    else {
      ASSERT_EQ(b[i-2], out_vec[i]);
    }
  }

  cb1.set_size(5, true);
  cb1.peek(out_vec, -1);
  for (int i = 0; i < out_vec.length(); ++i) {
    ASSERT_EQ(b[i+3], out_vec[i]);
  }

  // "Copy constructor: "
  Circular_Buffer<double> cb2(cb1);
  cb2.peek(out_vec, -1);
  for (int i = 0; i < out_vec.length(); ++i) {
    ASSERT_EQ(b[i+3], out_vec[i]);
  }

  // "Copy operator: "
  Circular_Buffer<double> cb3 = cb1;
  cb3.peek(out_array, -1);
  for (int i = 0; i < out_vec.length(); ++i) {
    ASSERT_EQ(b[i+3], out_vec[i]);
  }
}
