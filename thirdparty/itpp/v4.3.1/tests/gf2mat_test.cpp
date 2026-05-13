/*!
 * \file
 * \brief Test program for a class for algebra on GF(2) (binary) matrices
 * \author Erik G. Larsson
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

// To run extensive tests uncomment the following definition
// #define EXTENSIVE_TESTS

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


GF2mat random_matrix(int m, int n)
{
  GF2mat Z(m, n);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      Z.set(i, j, (rand() % 2 == 0 ? 1 : 0));
    }
  }
  return Z;
}

int main()
{

// ========== SIMPLE DEMOS/TESTS ====================

  cout << "===========================================================" << endl;
  cout << "gf2mat_test: Test of matrix operations in gfmat.h/gfmat.cpp" << endl;
  cout << "===========================================================" << endl;

  GF2mat A(3, 3);
  A.set(0, 0, 1);
  A.set(1, 2, 1);
  A.set(2, 1, 1);
  cout << "A=" << A << endl;
  cout << "A*A=" << A*A << endl;
  cout << "A*A'=" << A*A.transpose() << endl;

  GF2mat B = A;
  cout << "B=" << B << endl;
  cout << "B.get_row(1)=" << B.get_row(1) << endl;
  cout << "B.get_col(2)=" << B.get_col(2) << endl;

  bvec v(3);
  v(0) = 1;
  v(1) = 1;
  v(2) = 0;
  cout << "v=" << v << endl;
  cout << "A*v=" << A*v << endl;

  cout << "rank(A)=" << A.row_rank() << endl;
  cout << "inv(A)=" << A.inverse() << endl;

  GF2mat C, D;
  ivec p;
  A.T_fact(C, D, p);
  cout << "C=" << C << endl;
  cout << "D=" << D << endl;
  cout << "p=" << p << endl;

// Test Alist functionality
  string file = "gf2mat_test.alist";
  GF2mat_sparse_alist alist;
  alist.from_sparse(A.sparsify());
  alist.write(file);
  GF2mat_sparse_alist alist2(file);
  it_assert(GF2mat(alist2.to_sparse()) == A, "Alist test failed");

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
  cout << "Saved successfully." << endl;

  it_ifile f2("gf2mat_test.it");
  GF2mat Z_temp;
  f2 >> Name("Z") >> Z_temp;
  f2.close();
  it_assert(Z == Z_temp, "error");
  cout << "Read successfully." << endl;

// Binary vector
  bvec b = randb(Z.cols());
  it_assert(GF2mat(Z*b, 1) == Z*GF2mat(b, 1), "error");
  cout << "Passed matrix-vector multiplication test" << endl;

// Multiplication test
  GF2mat W = random_matrix(139, Z.rows());
  GF2mat temp1 = W * Z;
  cout << "computed product." << endl;
  GF2mat temp2 = GF2mat(W.sparsify() * Z.sparsify());
  it_assert(temp1 == temp2, "error");
  Z = Z.transpose();
  it_assert(W*Z.transpose() == mult_trans(W, Z), "error");
  cout << "Passed matrix-matrix multiplication test" << endl;

// Transpose
  it_assert(GF2mat(b, 0) == GF2mat(b, 1).transpose(), "error");
  it_assert(GF2mat(b, 1) == GF2mat(b, 0).transpose(), "error");
  GF2mat Y = random_matrix(Z.cols(), 73);
  it_assert((Z*Y).transpose() == Y.transpose()*Z.transpose(), "error");
  cout << "Passed transpose test." << endl;

// Concatenation
  int m = Z.rows();
  int n = Z.cols();
  it_assert(Z == Z.get_submatrix(0, 0, m - 1, 27).concatenate_horizontal(Z.get_submatrix(0, 28, m - 1, n - 1)), "error");
  it_assert(Z == Z.get_submatrix(0, 0, 13, n - 1).concatenate_vertical(Z.get_submatrix(14, 0, m - 1, n - 1)), "error");
  cout << "Passed concatenation test." << endl;

// Assignment operator
  GF2mat P = Z;
  it_assert(P == Z, "error");
  it_assert((P + Z).is_zero(), "error");
  cout << "Passed assignment operator test." << endl;

// Sparse-dense conversions
  GF2mat_sparse As(Z.rows(), Z.cols());
  for (int i = 0; i < Z.rows(); i++) {
    for (int j = 0; j < Z.cols(); j++) {
      if (Z.get(i, j) == 1) {
        As.set(i, j, 1);
      }
    }
  }
  it_assert(GF2mat(As) == Z, "error");
  GF2mat_sparse Cs = Z.sparsify();
  it_assert(Cs.full() == As.full(), "error");
  cout << "Passed sparse test." << endl;

  Z = random_matrix(100, 75);

// Get rows and columns
  cout << "Z.get_row(1)=" << Z.get_row(1) << endl;
  cout << "Z.get_col(2)=" << Z.get_col(2) << endl;

// Print a submatrix on the screen
  cout << "Z.get_submatrix(1,1,6,4): " << Z.get_submatrix(1, 1, 6, 4) << endl;


// Test of T-factorization
  int dim = 250;

  for (int trial = 0; trial < 100; trial++) {
    cout << "Testing T-factorization, realization: " << trial << endl;
    GF2mat X = random_matrix(rand() % dim + 1, rand() % dim + 1);
    GF2mat T, U;
    ivec perm;
    X.T_fact(T, U, perm);
    GF2mat W = T * X;
    W.permute_cols(perm, 0);
    it_assert(U == W, "error");
  }

// Test of inversion
  for (int trial = 0; trial < 100; trial++) {
    cout << "Testing inversion, realization: " << trial << endl;
    GF2mat X = random_matrix(dim, dim);
    while (X.row_rank() != dim) {
      X.set(rand() % dim, rand() % dim, rand() % 2);
    }
    it_assert(X*X.inverse() == gf2dense_eye(dim), "error");
    it_assert(X.inverse()*X == gf2dense_eye(dim), "error");
  }

// Test of the T-factorization bitflip update
  for (int trial = 0; trial < 100; trial++) {
    cout << "Testing the T-factorization bitflip update, realization: " << trial;
    GF2mat X = random_matrix(rand() % dim + 1, rand() % dim + 1);
    cout << " dimension: " << X.rows() << "*" << X.cols();
    GF2mat T, U;
    ivec perm;
    int rank = X.T_fact(T, U, perm);
    cout << " rank:" << rank << endl;

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
      it_assert(Unew == W, "error");
    }
  }


// Test of the T-factorization add-column update
  for (int trial = 0; trial < 100; trial++) {
    cout << "Testing the T-factorization add-column update, realization number: " << trial << endl;
    bvec c = randb(dim);
    GF2mat X(c, 1);

    GF2mat T, U;
    ivec perm;
    X.T_fact(T, U, perm);

    for (int trial2 = 0; trial2 < 100; trial2++) {
      bvec c = randb(dim);
      //      cout << X << endl;
      //      cout << GF2mat(c,1) << endl;
      GF2mat Xtemp = X.concatenate_horizontal(GF2mat(c, 1));
      int success = Xtemp.T_fact_update_addcol(T, U, perm, c);
      if (success == 1) {
        X = Xtemp;
      }
      //      cout << "rank was: " << X.row_rank() << endl;

      GF2mat W = T * X;
      W.permute_cols(perm, 0);
      it_assert(U == W, "error");
    }
  }

  cout << "All tests successfully passed." << endl;

#endif // #ifdef(EXTENSIVE_TESTS)

}
