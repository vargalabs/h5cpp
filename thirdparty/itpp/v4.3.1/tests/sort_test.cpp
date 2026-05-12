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
#include <iomanip>

using namespace itpp;
using namespace std;

int main()
{
  vec A;
  ivec B;
  ivec I;
  int failedSort;
  Uniform_RNG URNG;
  I_Uniform_RNG IURNG(0, 10);

  cout.setf(ios::fixed);
  cout.precision(3);

  cout << "-----------------" << endl;
  cout << "  Sorting Tests" << endl;
  cout << "-----------------" << endl << endl;
  cout << "Vectors of length 15, using uniform random doubles (0,1):" << endl;
  A = URNG(15);
  cout << "Introsort:" << endl;
  cout << "    A (before sort): " << A << endl;
  sort(A);
  cout << "    A (after sort):  " << A << endl;

  A = URNG(15);
  cout << "Quicksort:" << endl;
  cout << "    A (before sort): " << A << endl;
  sort(A, QUICKSORT);
  cout << "    A (after sort):  " << A << endl;

  A = URNG(15);
  cout << "Heapsort:" << endl;
  cout << "    A (before sort): " << A << endl;
  sort(A, HEAPSORT);
  cout << "    A (after sort):  " << A << endl;

  A = URNG(15);
  cout << "Insertion Sort:" << endl;
  cout << "    A (before sort): " << A << endl;
  sort(A, INSERTSORT);
  cout << "    A (after sort):  " << A << endl;

  A = URNG(15);
  cout << "Introsort (index):" << endl;
  cout << "    A:    " << A << endl;
  I = sort_index(A);
  cout << "    I:    " << I << endl;
  cout << "    A(I): " << A(I) << endl;

  A = URNG(15);
  cout << "Quicksort (index):" << endl;
  cout << "    A:    " << A << endl;
  I = sort_index(A, QUICKSORT);
  cout << "    I:    " << I << endl;
  cout << "    A(I): " << A(I) << endl;

  A = URNG(15);
  cout << "Heapsort (index):" << endl;
  cout << "    A:    " << A << endl;
  I = sort_index(A, HEAPSORT);
  cout << "    I:    " << I << endl;
  cout << "    A(I): " << A(I) << endl;

  A = URNG(15);
  cout << "Insertion sort (index):" << endl;
  cout << "    A:    " << A << endl;
  I = sort_index(A, INSERTSORT);
  cout << "    I:    " << I << endl;
  cout << "    A(I): " << A(I) << endl;

  cout << endl << endl;


  cout << "Vectors of length 20000, using uniform random doubles (0,1):" << endl;
  A = URNG(20000);
  sort(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Introsort:                Successful!" << endl;
  else
    cout << "Introsort:                Failed!" << endl;

  A = URNG(20000);
  sort(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Quicksort:                Successful!" << endl;
  else
    cout << "Quicksort:                Failed!" << endl;


  A = URNG(20000);
  sort(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(i - 1) > A(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Heapsort:                 Successful!" << endl;
  else
    cout << "Heapsort:                 Failed!" << endl;

  A = URNG(20000);
  I = sort_index(A);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Introsort (index):        Successful!" << endl;
  else
    cout << "Introsort (index):        Failed!" << endl;


  A = URNG(20000);
  I = sort_index(A, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Quicksort (index):        Successful!" << endl;
  else
    cout << "Quicksort (index):        Failed!" << endl;

  A = URNG(20000);
  I = sort_index(A, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < A.length(); i++) {
    if (A(I(i - 1)) > A(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Heapsort (index):         Successful!" << endl;
  else
    cout << "Heapsort (index):         Failed!" << endl;
  cout << endl << endl;

  cout << "Vectors of length 15, using uniform random integers (0,10):" << endl;
  B = IURNG(15);
  cout << "Introsort:" << endl;
  cout << "    B (before sort): " << B << endl;
  sort(B);
  cout << "    B (after sort):  " << B << endl;

  B = IURNG(15);
  cout << "Quicksort:" << endl;
  cout << "    B (before sort): " << B << endl;
  sort(B, QUICKSORT);
  cout << "    B (after sort):  " << B << endl;

  B = IURNG(15);
  cout << "Heapsort:" << endl;
  cout << "    B (before sort): " << B << endl;
  sort(B, HEAPSORT);
  cout << "    B (after sort):  " << B << endl;

  B = IURNG(15);
  cout << "Insertion sort:" << endl;
  cout << "    B (before sort): " << B << endl;
  sort(B, INSERTSORT);
  cout << "    B (after sort):  " << B << endl;

  B = IURNG(15);
  cout << "Introsort (index):" << endl;
  cout << "    B:    " << B << endl;
  I = sort_index(B);
  cout << "    I:    " << I << endl;
  cout << "    B(I): " << B(I) << endl;

  B = IURNG(15);
  cout << "Quicksort (index):" << endl;
  cout << "    B:    " << B << endl;
  I = sort_index(B, QUICKSORT);
  cout << "    I:    " << I << endl;
  cout << "    B(I): " << B(I) << endl;

  B = IURNG(15);
  cout << "Heapsort (index):" << endl;
  cout << "    B:    " << B << endl;
  I = sort_index(B, HEAPSORT);
  cout << "    I:    " << I << endl;
  cout << "    B(I): " << B(I) << endl;

  B = IURNG(15);
  cout << "Insertion sort (index):" << endl;
  cout << "    B:    " << B << endl;
  I = sort_index(B, INSERTSORT);
  cout << "    I:    " << I << endl;
  cout << "    B(I): " << B(I) << endl;

  cout << endl << endl;

  cout << "Vectors of length 20000, using uniform random integers (0,10):" << endl;
  B = IURNG(20000);
  sort(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Introsort:                Successful!" << endl;
  else
    cout << "Introsort:                Failed!" << endl;

  B = IURNG(20000);
  sort(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Quicksort:                Successful!" << endl;
  else
    cout << "Quicksort:                Failed!" << endl;


  B = IURNG(20000);
  sort(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(i - 1) > B(i)) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Heapsort:                 Successful!" << endl;
  else
    cout << "Heapsort:                 Failed!" << endl;

  B = IURNG(20000);
  I = sort_index(B);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Introsort (index):        Successful!" << endl;
  else
    cout << "Introsort (index):        Failed!" << endl;

  B = IURNG(20000);
  I = sort_index(B, QUICKSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Quicksort (index):        Successful!" << endl;
  else
    cout << "Quicksort (index):        Failed!" << endl;

  B = IURNG(20000);
  I = sort_index(B, HEAPSORT);
  failedSort = 0;
  for (int i = 1; i < B.length(); i++) {
    if (B(I(i - 1)) > B(I(i))) {
      cout << "error in sort!" << endl;
      failedSort++;
    }
  }
  if (failedSort == 0)
    cout << "Heapsort (index):         Successful!" << endl;
  else
    cout << "Heapsort (index):         Failed!" << endl;
}
