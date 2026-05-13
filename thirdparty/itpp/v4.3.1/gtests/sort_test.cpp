/*!
 * \file
 * \brief Sort class test program
 * \author Mark Dobossy
 *
 * $Date$
 * $Revision$
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

TEST(Sort, All)
{
  vec A;
  ivec B;
  ivec I;
  int failedSort;
  Uniform_RNG URNG;
  I_Uniform_RNG IURNG(0, 10);

  // Sorting Tests
  // Vectors of length 15, using uniform random doubles (0,1)
  A = URNG(15);
  sort(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  sort(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  sort(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  sort(A, INSERTSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  I = sort_index(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  I = sort_index(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  I = sort_index(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(15);
  I = sort_index(A, INSERTSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  // Vectors of length 20000, using uniform random doubles (0,1)
  A = URNG(20000);
  sort(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(20000);
  sort(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(20000);
  sort(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(20000);
  I = sort_index(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(20000);
  I = sort_index(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  A = URNG(20000);
  I = sort_index(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  // Vectors of length 15, using uniform random integers (0,10)
  B = IURNG(15);
  sort(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  sort(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  sort(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  sort(B, INSERTSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  I = sort_index(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  I = sort_index(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  I = sort_index(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(15);
  I = sort_index(B, INSERTSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  // Vectors of length 20000, using uniform random integers (0,10)
  B = IURNG(20000);
  sort(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(20000);
  sort(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(20000);
  sort(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(20000);
  I = sort_index(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(20000);
  I = sort_index(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);

  B = IURNG(20000);
  I = sort_index(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      failedSort++;
    }
  }
  ASSERT_TRUE(0 == failedSort);
}
