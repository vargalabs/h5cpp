/*!
 * \file
 * \brief Test program for a class for algebra on GF(2) (binary) matrices
 * \author Erik G. Larsson
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

// To run extensive tests uncomment the following definition
//#define EXTENSIVE_TESTS

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
GF2mat random_matrix(int m, int n)
{
  GF2mat Z(m, n);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      Z.set(i, j, randi(0, 1));
    }
  }
  return Z;
}

TEST (GF2Mat, All)
{
  RNG_reset(0);
  int i,j;

  // gf2mat_test: Test of matrix operations in gfmat.h/gfmat.cpp

  GF2mat A(3, 3);
  A.set(0, 0, 1);
  A.set(1, 2, 1);
  A.set(2, 1, 1);

  ASSERT_DOUBLE_EQ(3.0/9, A.density());
  ASSERT_EQ(3, A.rows());
  ASSERT_EQ(3, A.cols());
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if ((0 == i && 0 == j) || (1 == i && 2 == j) || (2 == i && 1 == j)) {
        ASSERT_EQ(bin(1), A(i,j));
      }
      else {
        ASSERT_EQ(bin(0), A(i,j));
      }
    }
  }

  GF2mat B = A*A;
  ASSERT_DOUBLE_EQ(3.0/9, B.density());
  ASSERT_EQ(3, B.rows());
  ASSERT_EQ(3, B.cols());
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if ((0 == i && 0 == j) || (1 == i && 1 == j) || (2 == i && 2 == j)) {
        ASSERT_EQ(bin(1), B(i,j));
      }
      else {
        ASSERT_EQ(bin(0), B(i,j));
      }
    }
  }

  B = A*A.transpose();
  ASSERT_DOUBLE_EQ(3.0/9, B.density());
  ASSERT_EQ(3, B.rows());
  ASSERT_EQ(3, B.cols());
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if ((0 == i && 0 == j) || (1 == i && 1 == j) || (2 == i && 2 == j)) {
        ASSERT_EQ(bin(1), B(i,j));
      }
      else {
        ASSERT_EQ(bin(0), B(i,j));
      }
    }
  }

  B = A;
  bvec v = B.get_row(1);
  bvec v_ref = "0 0 1";
  ASSERT_TRUE(v == v_ref);
  v = B.get_col(2);
  v_ref = "0 1 0";
  ASSERT_TRUE(v == v_ref);

  v.set_size(3);
  v(0) = 1;
  v(1) = 1;
  v(2) = 0;
  v = A*v;
  v_ref = "1 0 1";
  ASSERT_TRUE(v == v_ref);

  ASSERT_EQ(3, A.row_rank());
  B = A.inverse();
  ASSERT_DOUBLE_EQ(3.0/9, B.density());
  ASSERT_EQ(3, B.rows());
  ASSERT_EQ(3, B.cols());
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if ((0 == i && 0 == j) || (1 == i && 2 == j) || (2 == i && 1 == j)) {
        ASSERT_EQ(bin(1), B(i,j));
      }
      else {
        ASSERT_EQ(bin(0), B(i,j));
      }
    }
  }

  GF2mat C, D;
  ivec p;
  A.T_fact(C, D, p);
  ASSERT_DOUBLE_EQ(3.0/9, C.density());
  ASSERT_EQ(3, C.rows());
  ASSERT_EQ(3, C.cols());
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      if ((0 == i && 0 == j) || (1 == i && 1 == j) || (2 == i && 2 == j)) {
        ASSERT_EQ(bin(1), C(i,j));
      }
      else {
        ASSERT_EQ(bin(0), C(i,j));
      }
    }
  }
  ASSERT_EQ(C, D);
  ivec p_ref = "0 2 1";
  ASSERT_EQ(p, p_ref);

// Test Alist functionality
  string file = "gf2mat_test.alist";
  GF2mat_sparse_alist alist;
  alist.from_sparse(A.sparsify());
  alist.write(file);
  GF2mat_sparse_alist alist2(file);
  ASSERT_EQ(GF2mat(alist2.to_sparse()), A) << "Alist test failed";

#ifdef EXTENSIVE_TESTS

// ========== EXTENSIVE RANDOM TESTS ==========

// The following code can be used to verify the behavior of the GF2
// class for large matrix dimensions.  Note that with debugging
// options enabled, this code takes a while to run.  To run these
// extensive tests, simply remove the comments around this code segment.

// Test of file I/O
  GF2mat Z = random_matrix(301, 179);

  it_file f1("gf2mat_test.it");
  f1 << Name("Z") << Z;
  f1.close();

  it_ifile f2("gf2mat_test.it");
  GF2mat Z_temp;
  f2 >> Name("Z") >> Z_temp;
  f2.close();
  ASSERT_EQ(Z, Z_temp);

// Binary vector
  bvec b = randb(Z.cols());
  ASSERT_EQ(GF2mat(Z*b, 1), Z*GF2mat(b, 1));

// Multiplication test
  GF2mat W = random_matrix(139, Z.rows());
  GF2mat temp1 = W * Z;
  GF2mat temp2 = GF2mat(W.sparsify() * Z.sparsify());
  ASSERT_EQ(temp1, temp2);
  Z = Z.transpose();
  ASSERT_EQ(W*Z.transpose(), mult_trans(W, Z));

// Transpose
  ASSERT_EQ(GF2mat(b, 0), GF2mat(b, 1).transpose());
  ASSERT_EQ(GF2mat(b, 1), GF2mat(b, 0).transpose());
  GF2mat Y = random_matrix(Z.cols(), 73);
  ASSERT_EQ((Z*Y).transpose(), Y.transpose()*Z.transpose());

// Concatenation
  int m = Z.rows();
  int n = Z.cols();
  ASSERT_EQ(Z, Z.get_submatrix(0, 0, m - 1, 27).concatenate_horizontal(Z.get_submatrix(0, 28, m - 1, n - 1)));
  ASSERT_EQ(Z, Z.get_submatrix(0, 0, 13, n - 1).concatenate_vertical(Z.get_submatrix(14, 0, m - 1, n - 1)));

// Assignment operator
  GF2mat P = Z;
  ASSERT_EQ(P, Z);
  ASSERT_TRUE((P + Z).is_zero());

// Sparse-dense conversions
  GF2mat_sparse As(Z.rows(), Z.cols());
  for (i = 0; i < Z.rows(); i++) {
    for (j = 0; j < Z.cols(); j++) {
      if (Z.get(i, j) == 1) {
        As.set(i, j, 1);
      }
    }
  }
  ASSERT_EQ(GF2mat(As), Z);
  GF2mat_sparse Cs = Z.sparsify();
  ASSERT_EQ(Cs.full(), As.full());

  Z = random_matrix(100, 75);

  // Get rows and columns
  v_ref = "0 0 1 1 1 1 1 1 1 0 0 0 0 1 1 0 1 0 0 1 0 1 0 1 1 1 0 0 1 0 1 1 1 1 0 0 0 0 0 1 0 0 1 0 0 1 0 "
  "1 0 1 0 1 1 1 1 0 1 0 1 1 0 0 1 1 1 1 0 0 0 1 0 0 1 0 1";
  ASSERT_EQ(v_ref, Z.get_row(1));
  v_ref = "0 1 0 1 1 1 1 1 0 1 0 0 1 1 1 0 1 0 0 1 1 1 1 1 0 0 1 1 1 0 1 0 0 1 0 0 0 1 0 0 0 1 0 0 1 1 0 "
  "1 0 1 0 1 1 0 1 0 0 0 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 0 0 0 1 1 0 0 1 1 1 1 0 0 0 1 0 0 1 1 1 1 1 1 0 0 1 0";
  ASSERT_EQ(v_ref, Z.get_col(2));

  // Print a submatrix on the screen
  B = Z.get_submatrix(1, 1, 6, 4);
  ASSERT_DOUBLE_EQ(0.54166666666666663, B.density());
  ASSERT_EQ(6, B.rows());
  ASSERT_EQ(4, B.cols());

// Test of T-factorization
  int dim = 250;

  for (int trial = 0; trial < 100; trial++) {
    GF2mat X = random_matrix(rand() % dim + 1, rand() % dim + 1);
    GF2mat T, U;
    ivec perm;
    X.T_fact(T, U, perm);
    GF2mat W = T * X;
    W.permute_cols(perm, 0);
    ASSERT_EQ(U, W);
  }

// Test of inversion
  for (int trial = 0; trial < 100; trial++) {
    GF2mat X = random_matrix(dim, dim);
    while (X.row_rank() != dim) {
      X.set(rand() % dim, rand() % dim, rand() % 2);
    }
    ASSERT_EQ(X*X.inverse(), gf2dense_eye(dim));
    ASSERT_EQ(X.inverse()*X, gf2dense_eye(dim));
  }

// Test of the T-factorization bitflip update
  for (int trial = 0; trial < 100; trial++) {
    GF2mat X = random_matrix(rand() % dim + 1, rand() % dim + 1);
    GF2mat T, U;
    ivec perm;
    int rank = X.T_fact(T, U, perm);

    GF2mat Tnew = T;
    GF2mat Unew = U;
    ivec permnew = perm;
    for (int trial2 = 0; trial2 < 10; trial2++) {
      int i = rand() % X.rows();
      int j = rand() % X.cols();
      X.addto_element(i, j, 1);
      X.T_fact_update_bitflip(Tnew, Unew, permnew, rank, i, j);

      GF2mat W = Tnew * X;
      W.permute_cols(permnew, 0);
      ASSERT_EQ(Unew, W);
    }
  }


// Test of the T-factorization add-column update
  for (int trial = 0; trial < 100; trial++) {
    bvec c = randb(dim);
    GF2mat X(c, 1);

    GF2mat T, U;
    ivec perm;
    X.T_fact(T, U, perm);

    for (int trial2 = 0; trial2 < 100; trial2++) {
      bvec c = randb(dim);
      //      cerr << X << endl;
      //      cerr << GF2mat(c,1) << endl;
      GF2mat Xtemp = X.concatenate_horizontal(GF2mat(c, 1));
      int success = Xtemp.T_fact_update_addcol(T, U, perm, c);
      if (success == 1) {
        X = Xtemp;
      }
      //      cerr << "rank was: " << X.row_rank() << endl;

      GF2mat W = T * X;
      W.permute_cols(perm, 0);
      ASSERT_EQ(U, W);
    }
  }

#endif // #ifdef(EXTENSIVE_TESTS)

}
