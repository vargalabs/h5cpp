/*!
 * \file
 * \brief Implementation of a class for algebra on GF(2) (binary) matrices
 * \author Erik G. Larsson and Adam Piatyszek
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

#include <itpp/base/gf2mat.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/converters.h>
#include <iostream>

namespace itpp
{

// ====== IMPLEMENTATION OF THE ALIST CLASS ==========

GF2mat_sparse_alist::GF2mat_sparse_alist(const std::string &fname)
    : data_ok(false)
{
  read(fname);
}

void GF2mat_sparse_alist::read(const std::string &fname)
{
  std::ifstream file;
  std::string line;
  std::stringstream ss;

  file.open(fname.c_str());
  it_assert(file.is_open(),
            "GF2mat_sparse_alist::read(): Could not open file \""
            << fname << "\" for reading");

  // parse N and M:
  getline(file, line);
  ss << line;
  ss >> N >> M;
  it_assert(!ss.fail(),
            "GF2mat_sparse_alist::read(): Wrong alist data (N or M)");
  it_assert((N > 0) && (M > 0),
            "GF2mat_sparse_alist::read(): Wrong alist data");
  ss.seekg(0, std::ios::end);
  ss.clear();

  // parse max_num_n and max_num_m
  getline(file, line);
  ss << line;
  ss >> max_num_n >> max_num_m;
  it_assert(!ss.fail(),
            "GF2mat_sparse_alist::read(): Wrong alist data (max_num_{n,m})");
  it_assert((max_num_n >= 0) && (max_num_n <= N) &&
            (max_num_m >= 0) && (max_num_m <= M),
            "GF2mat_sparse_alist::read(): Wrong alist data");
  ss.seekg(0, std::ios::end);
  ss.clear();

  // parse weight of each column n
  num_nlist.set_size(N);
  num_nlist.clear();
  getline(file, line);
  ss << line;
  for (int i = 0; i < N; i++) {
    ss >> num_nlist(i);
    it_assert(!ss.fail(),
              "GF2mat_sparse_alist::read(): Wrong alist data (num_nlist("
              << i << "))");
    it_assert((num_nlist(i) >= 0) && (num_nlist(i) <= M),
              "GF2mat_sparse_alist::read(): Wrong alist data (num_nlist("
              << i << "))");
  }
  ss.seekg(0, std::ios::end);
  ss.clear();

  // parse weight of each row m
  num_mlist.set_size(M);
  num_mlist.clear();
  getline(file, line);
  ss << line;
  for (int i = 0; i < M; i++) {
    ss >> num_mlist(i);
    it_assert(!ss.fail(),
              "GF2mat_sparse_alist::read(): Wrong alist data (num_mlist("
              << i << "))");
    it_assert((num_mlist(i) >= 0) && (num_mlist(i) <= N),
              "GF2mat_sparse_alist::read(): Wrong alist data (num_mlist("
              << i << "))");
  }
  ss.seekg(0, std::ios::end);
  ss.clear();

  // parse coordinates in the n direction with non-zero entries
  nlist.set_size(N, max_num_n);
  nlist.clear();
  for (int i = 0; i < N; i++) {
    getline(file, line);
    ss << line;
    for (int j = 0; j < num_nlist(i); j++) {
      ss >> nlist(i, j);
      it_assert(!ss.fail(),
                "GF2mat_sparse_alist::read(): Wrong alist data (nlist("
                << i << "," << j << "))");
      it_assert((nlist(i, j) > 0) && (nlist(i, j) <= M),
                "GF2mat_sparse_alist::read(): Wrong alist data (nlist("
                << i << "," << j << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
  }

  // parse coordinates in the m direction with non-zero entries
  mlist.set_size(M, max_num_m);
  mlist.clear();
  for (int i = 0; i < M; i++) {
    getline(file, line);
    ss << line;
    for (int j = 0; j < num_mlist(i); j++) {
      ss >> mlist(i, j);
      it_assert(!ss.fail(),
                "GF2mat_sparse_alist::read(): Wrong alist data (mlist("
                << i << "," << j << "))");
      it_assert((mlist(i, j) > 0) && (mlist(i, j) <= N),
                "GF2mat_sparse_alist::read(): Wrong alist data (mlist("
                << i << "," << j << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
  }

  file.close();
  data_ok = true;
}

void GF2mat_sparse_alist::write(const std::string &fname) const
{
  it_assert(data_ok,
            "GF2mat_sparse_alist::write(): alist data not ready for writing");

  std::ofstream file(fname.c_str(), std::ofstream::out);
  it_assert(file.is_open(),
            "GF2mat_sparse_alist::write(): Could not open file \""
            << fname << "\" for writing");

  file << N << " " << M << std::endl;
  file << max_num_n << " " << max_num_m << std::endl;

  for (int i = 0; i < num_nlist.length() - 1; i++)
    file << num_nlist(i) << " ";
  file << num_nlist(num_nlist.length() - 1) << std::endl;

  for (int i = 0; i < num_mlist.length() - 1; i++)
    file << num_mlist(i) << " ";
  file << num_mlist(num_mlist.length() - 1) << std::endl;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < num_nlist(i) - 1; j++)
      file << nlist(i, j) << " ";
    file << nlist(i, num_nlist(i) - 1) << std::endl;
  }

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < num_mlist(i) - 1; j++)
      file << mlist(i, j) << " ";
    file << mlist(i, num_mlist(i) - 1) << std::endl;
  }

  file.close();
}


GF2mat_sparse GF2mat_sparse_alist::to_sparse(bool transpose) const
{
  GF2mat_sparse sbmat(M, N, max_num_m);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < num_mlist(i); j++) {
      sbmat.set_new(i, mlist(i, j) - 1, bin(1));
    }
  }
  sbmat.compact();

  if (transpose) {
    return sbmat.transpose();
  }
  else {
    return sbmat;
  }
}


// ----------------------------------------------------------------------
// WARNING: This method is very slow. Sparse_Mat has to be extended with
// some extra functions to improve the performance of this.
// ----------------------------------------------------------------------
void GF2mat_sparse_alist::from_sparse(const GF2mat_sparse &sbmat,
                                      bool transpose)
{
  if (transpose) {
    from_sparse(sbmat.transpose(), false);
  }
  else {
    // check matrix dimension
    M = sbmat.rows();
    N = sbmat.cols();

    num_mlist.set_size(M);
    num_nlist.set_size(N);

    // fill mlist matrix, num_mlist vector and max_num_m
    mlist.set_size(0, 0);
    int tmp_cols = 0; // initial number of allocated columns
    for (int i = 0; i < M; i++) {
      ivec temp_row(0);
      for (int j = 0; j < N; j++) {
        if (sbmat(i, j) == bin(1)) {
          temp_row = concat(temp_row, j + 1);
        }
      }
      int trs = temp_row.size();
      if (trs > tmp_cols) {
        tmp_cols = trs;
        mlist.set_size(M, tmp_cols, true);
      }
      else if (trs < tmp_cols) {
        temp_row.set_size(tmp_cols, true);
      }
      mlist.set_row(i, temp_row);
      num_mlist(i) = trs;
    }
    max_num_m = max(num_mlist);

    // fill nlist matrix, num_nlist vector and max_num_n
    nlist.set_size(0, 0);
    tmp_cols = 0; // initial number of allocated columns
    for (int j = 0; j < N; j++) {
      ivec temp_row = sbmat.get_col(j).get_nz_indices() + 1;
      int trs = temp_row.size();
      if (trs > tmp_cols) {
        tmp_cols = trs;
        nlist.set_size(N, tmp_cols, true);
      }
      else if (trs < tmp_cols) {
        temp_row.set_size(tmp_cols, true);
      }
      nlist.set_row(j, temp_row);
      num_nlist(j) = trs;
    }
    max_num_n = max(num_nlist);

    data_ok = true;
  }
}


// ----------------------------------------------------------------------
// Implementation of a dense GF2 matrix class
// ----------------------------------------------------------------------

GF2mat::GF2mat(int i, int j): nrows(i), ncols(j),
    nwords((j >> shift_divisor) + 1)
{
  data.set_size(nrows, nwords);
  data.clear();
}

GF2mat::GF2mat(): nrows(1), ncols(1), nwords(1)
{
  data.set_size(nrows, nwords);
  data.clear();
}

GF2mat::GF2mat(const bvec &x, bool is_column)
{
  if (is_column) {     // create column vector
    nrows = length(x);
    ncols = 1;
    nwords = 1;
    data.set_size(nrows, nwords);
    data.clear();
    for (int i = 0; i < nrows; i++) {
      set(i, 0, x(i));
    }
  }
  else {    // create row vector
    nrows = 1;
    ncols = length(x);
    nwords = (ncols >> shift_divisor) + 1;
    data.set_size(nrows, nwords);
    data.clear();
    for (int i = 0; i < ncols; i++) {
      set(0, i, x(i));
    }
  }
}


GF2mat::GF2mat(const bmat &X): nrows(X.rows()), ncols(X.cols())
{
  nwords = (ncols >> shift_divisor) + 1;
  data.set_size(nrows, nwords);
  data.clear();
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      set(i, j, X(i, j));
    }
  }
}


GF2mat::GF2mat(const GF2mat_sparse &X)
{
  nrows = X.rows();
  ncols = X.cols();
  nwords = (ncols >> shift_divisor) + 1;
  data.set_size(nrows, nwords);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      data(i, j) = 0;
    }
  }

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < X.get_col(j).nnz(); i++)   {
      bin b = X.get_col(j).get_nz_data(i);
      set(X.get_col(j).get_nz_index(i), j, b);
    }
  }
}

GF2mat::GF2mat(const GF2mat_sparse &X, int m1, int n1, int m2, int n2)
{
  it_assert(X.rows() > m2, "GF2mat(): indexes out of range");
  it_assert(X.cols() > n2, "GF2mat(): indexes out of range");
  it_assert(m1 >= 0 && n1 >= 0 && m2 >= m1 && n2 >= n1,
            "GF2mat::GF2mat(): indexes out of range");

  nrows = m2 - m1 + 1;
  ncols = n2 - n1 + 1;
  nwords = (ncols >> shift_divisor) + 1;
  data.set_size(nrows, nwords);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      data(i, j) = 0;
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++)   {
      bin b = X(i + m1, j + n1);
      set(i, j, b);
    }
  }
}


GF2mat::GF2mat(const GF2mat_sparse &X, const ivec &columns)
{
  it_assert(X.cols() > max(columns),
            "GF2mat::GF2mat(): index out of range");
  it_assert(min(columns) >= 0,
            "GF2mat::GF2mat(): column index must be positive");

  nrows = X.rows();
  ncols = length(columns);
  nwords = (ncols >> shift_divisor) + 1;
  data.set_size(nrows, nwords);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      data(i, j) = 0;
    }
  }

  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < X.get_col(columns(j)).nnz(); i++)   {
      bin b = X.get_col(columns(j)).get_nz_data(i);
      set(X.get_col(columns(j)).get_nz_index(i), j, b);
    }
  }
}


void GF2mat::set_size(int m, int n, bool copy)
{
  nrows = m;
  ncols = n;
  nwords = (ncols >> shift_divisor) + 1;
  data.set_size(nrows, nwords, copy);
  if (!copy)
    data.clear();
}


GF2mat_sparse GF2mat::sparsify() const
{
  GF2mat_sparse Z(nrows, ncols);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (get(i, j) == 1) {
        Z.set(i, j, 1);
      }
    }
  }

  return Z;
}

bvec GF2mat::bvecify() const
{
  it_assert(nrows == 1 || ncols == 1,
            "GF2mat::bvecify() matrix must be a vector");
  int n = (nrows == 1 ? ncols : nrows);
  bvec result(n);
  if (nrows == 1) {
    for (int i = 0; i < n; i++) {
      result(i) = get(0, i);
    }
  }
  else {
    for (int i = 0; i < n; i++) {
      result(i) = get(i, 0);
    }
  }
  return result;
}


void GF2mat::set_row(int i, bvec x)
{
  it_assert(length(x) == ncols,
            "GF2mat::set_row(): dimension mismatch");
  for (int j = 0; j < ncols; j++) {
    set(i, j, x(j));
  }
}

void GF2mat::set_col(int j, bvec x)
{
  it_assert(length(x) == nrows,
            "GF2mat::set_col(): dimension mismatch");
  for (int i = 0; i < nrows; i++) {
    set(i, j, x(i));
  }
}


GF2mat gf2dense_eye(int m)
{
  GF2mat Z(m, m);
  for (int i = 0; i < m; i++) {
    Z.set(i, i, 1);
  }
  return Z;
}

GF2mat GF2mat::get_submatrix(int m1, int n1, int m2, int n2) const
{
  it_assert(m1 >= 0 && n1 >= 0 && m2 >= m1 && n2 >= n1
            && m2 < nrows && n2 < ncols,
            "GF2mat::get_submatrix() index out of range");
  GF2mat result(m2 - m1 + 1, n2 - n1 + 1);

  for (int i = m1; i <= m2; i++) {
    for (int j = n1; j <= n2; j++) {
      result.set(i - m1, j - n1, get(i, j));
    }
  }

  return result;
}


GF2mat GF2mat::concatenate_vertical(const GF2mat &X) const
{
  it_assert(X.ncols == ncols,
            "GF2mat::concatenate_vertical(): dimension mismatch");
  it_assert(X.nwords == nwords,
            "GF2mat::concatenate_vertical(): dimension mismatch");

  GF2mat result(nrows + X.nrows, ncols);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      result.data(i, j) = data(i, j);
    }
  }

  for (int i = 0; i < X.nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      result.data(i + nrows, j) = X.data(i, j);
    }
  }

  return result;
}

GF2mat GF2mat::concatenate_horizontal(const GF2mat &X) const
{
  it_assert(X.nrows == nrows,
            "GF2mat::concatenate_horizontal(): dimension mismatch");

  GF2mat result(nrows, X.ncols + ncols);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      result.set(i, j, get(i, j));
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < X.ncols; j++) {
      result.set(i, j + ncols, X.get(i, j));
    }
  }

  return result;
}

bvec GF2mat::get_row(int i) const
{
  bvec result(ncols);
  for (int j = 0; j < ncols; j++) {
    result(j) = get(i, j);
  }

  return result;
}

bvec GF2mat::get_col(int j) const
{
  bvec result(nrows);
  for (int i = 0; i < nrows; i++) {
    result(i) = get(i, j);
  }

  return result;
}


int GF2mat::T_fact(GF2mat &T, GF2mat &U, ivec &perm) const
{
  T = gf2dense_eye(nrows);
  U = *this;

  perm = zeros_i(ncols);
  for (int i = 0; i < ncols; i++) {
    perm(i) = i;
  }

  if (nrows > 250) {   // avoid cluttering output ...
    it_info_debug("Performing T-factorization of GF(2) matrix...  rows: "
                  << nrows << " cols: "  << ncols << " .... " << std::endl);
  }
  int pdone = 0;
  for (int j = 0; j < nrows; j++) {
    // Now working on diagonal element j,j
    // First try find a row with a 1 in column i
    int i1, j1;
    for (i1 = j; i1 < nrows; i1++) {
      for (j1 = j; j1 < ncols; j1++) {
        if (U.get(i1, j1) == 1) { goto found; }
      }
    }

    return j;

  found:
    U.swap_rows(i1, j);
    T.swap_rows(i1, j);
    U.swap_cols(j1, j);

    int temp = perm(j);
    perm(j) = perm(j1);
    perm(j1) = temp;

    // now subtract row i from remaining rows
    for (int i1 = j + 1; i1 < nrows; i1++)  {
      if (U.get(i1, j) == 1) {
        int ptemp = floor_i(100.0 * (i1 + j * nrows) / (nrows * nrows));
        if (nrows > 250 && ptemp > pdone + 10) {
          it_info_debug(ptemp << "% done.");
          pdone = ptemp;
        }
        U.add_rows(i1, j);
        T.add_rows(i1, j);
      }
    }
  }
  return nrows;
}


int GF2mat::T_fact_update_bitflip(GF2mat &T, GF2mat &U,
                                  ivec &perm, int rank,
                                  int r, int c) const
{
  // First, update U (before re-triangulization)
  int c0 = c;
  for (c = 0; c < ncols; c++) {
    if (perm(c) == c0) {
      goto foundidx;
    }
  }
  it_error("GF2mat::T_fact_update_bitflip() - internal error");

foundidx:
  for (int i = 0; i < nrows; i++) {
    if (T.get(i, r) == 1) {
      U.addto_element(i, c, 1);
    }
  }

  // first move column c to the end
  bvec lastcol = U.get_col(c);
  int temp_perm = perm(c);
  for (int j = c; j < ncols - 1; j++) {
    U.set_col(j, U.get_col(j + 1));
    perm(j) = perm(j + 1);
  }
  U.set_col(ncols - 1, lastcol);
  perm(ncols - 1) = temp_perm;

  // then, if the matrix is tall, also move row c to the end
  if (nrows >= ncols) {
    bvec lastrow_U = U.get_row(c);
    bvec lastrow_T = T.get_row(c);
    for (int i = c; i < nrows - 1; i++) {
      U.set_row(i, U.get_row(i + 1));
      T.set_row(i, T.get_row(i + 1));
    }
    U.set_row(nrows - 1, lastrow_U);
    T.set_row(nrows - 1, lastrow_T);

    // Do Gaussian elimination on the last row
    for (int j = c; j < ncols; j++)  {
      if (U.get(nrows - 1, j) == 1) {
        U.add_rows(nrows - 1, j);
        T.add_rows(nrows - 1, j);
      }
    }
  }

  // Now, continue T-factorization from the point (rank-1,rank-1)
  for (int j = rank - 1; j < nrows; j++) {
    int i1, j1;
    for (i1 = j; i1 < nrows; i1++) {
      for (j1 = j; j1 < ncols; j1++) {
        if (U.get(i1, j1) == 1) {
          goto found;
        }
      }
    }

    return j;

  found:
    U.swap_rows(i1, j);
    T.swap_rows(i1, j);
    U.swap_cols(j1, j);

    int temp = perm(j);
    perm(j) = perm(j1);
    perm(j1) = temp;

    for (int i1 = j + 1; i1 < nrows; i1++)  {
      if (U.get(i1, j) == 1) {
        U.add_rows(i1, j);
        T.add_rows(i1, j);
      }
    }
  }

  return nrows;
}

bool GF2mat::T_fact_update_addcol(GF2mat &T, GF2mat &U,
                                  ivec &perm, bvec newcol) const
{
  int i0 = T.rows();
  int j0 = U.cols();
  it_assert(j0 > 0, "GF2mat::T_fact_update_addcol(): dimension mismatch");
  it_assert(i0 == T.cols(),
            "GF2mat::T_fact_update_addcol(): dimension mismatch");
  it_assert(i0 == U.rows(),
            "GF2mat::T_fact_update_addcol(): dimension mismatch");
  it_assert(length(perm) == j0,
            "GF2mat::T_fact_update_addcol(): dimension mismatch");
  it_assert(U.get(j0 - 1, j0 - 1) == 1,
            "GF2mat::T_fact_update_addcol(): dimension mismatch");
  // The following test is VERY TIME-CONSUMING
  it_assert_debug(U.row_rank() == j0,
                  "GF2mat::T_fact_update_addcol(): factorization has incorrect rank");

  bvec z = T * newcol;
  GF2mat Utemp = U.concatenate_horizontal(GF2mat(z, 1));

  // start working on position (j0,j0)
  int i;
  for (i = j0; i < i0; i++) {
    if (Utemp.get(i, j0) == 1) {
      goto found;
    }
  }
  return (false); // adding the new column would not improve the rank

found:
  perm.set_length(j0 + 1, true);
  perm(j0) = j0;

  Utemp.swap_rows(i, j0);
  T.swap_rows(i, j0);

  for (int i1 = j0 + 1; i1 < i0; i1++)  {
    if (Utemp.get(i1, j0) == 1) {
      Utemp.add_rows(i1, j0);
      T.add_rows(i1, j0);
    }
  }

  U = Utemp;
  return (true); // the new column was successfully added
}




GF2mat GF2mat::inverse() const
{
  it_assert(nrows == ncols, "GF2mat::inverse(): Matrix must be square");

  // first compute the T-factorization
  GF2mat T, U;
  ivec perm;
  int rank = T_fact(T, U, perm);
  it_assert(rank == ncols, "GF2mat::inverse(): Matrix is not full rank");

  // backward substitution
  for (int i = ncols - 2; i >= 0; i--) {
    for (int j = ncols - 1; j > i; j--) {
      if (U.get(i, j) == 1) {
        U.add_rows(i, j);
        T.add_rows(i, j);
      }
    }
  }
  T.permute_rows(perm, 1);
  return T;
}

int GF2mat::row_rank() const
{
  GF2mat T, U;
  ivec perm;
  int rank = T_fact(T, U, perm);
  return rank;
}

bool GF2mat::is_zero() const
{
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      if (data(i, j) != 0) {
        return false;
      }
    }
  }
  return true;
}

bool GF2mat::operator==(const GF2mat &X) const
{
  if (X.nrows != nrows) { return false; }
  if (X.ncols != ncols) { return false; }
  it_assert(X.nwords == nwords, "GF2mat::operator==() dimension mismatch");

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nwords; j++) {
      //      if (X.get(i,j)!=get(i,j)) {
      if (X.data(i, j) != data(i, j)) {
        return false;
      }
    }
  }
  return true;
}


void GF2mat::add_rows(int i, int j)
{
  it_assert(i >= 0 && i < nrows, "GF2mat::add_rows(): out of range");
  it_assert(j >= 0 && j < nrows, "GF2mat::add_rows(): out of range");
  for (int k = 0; k < nwords; k++)   {
    data(i, k) ^= data(j, k);
  }
}

void GF2mat::swap_rows(int i, int j)
{
  it_assert(i >= 0 && i < nrows, "GF2mat::swap_rows(): index out of range");
  it_assert(j >= 0 && j < nrows, "GF2mat::swap_rows(): index out of range");
  for (int k = 0; k < nwords; k++) {
    unsigned char temp = data(i, k);
    data(i, k) = data(j, k);
    data(j, k) = temp;
  }
}

void GF2mat::swap_cols(int i, int j)
{
  it_assert(i >= 0 && i < ncols, "GF2mat::swap_cols(): index out of range");
  it_assert(j >= 0 && j < ncols, "GF2mat::swap_cols(): index out of range");
  bvec temp = get_col(i);
  set_col(i, get_col(j));
  set_col(j, temp);
}


void GF2mat::operator=(const GF2mat &X)
{
  nrows = X.nrows;
  ncols = X.ncols;
  nwords = X.nwords;
  data = X.data;
}

GF2mat operator*(const GF2mat &X, const GF2mat &Y)
{
  it_assert(X.ncols == Y.nrows, "GF2mat::operator*(): dimension mismatch");
  it_assert(X.nwords > 0, "Gfmat::operator*(): dimension mismatch");
  it_assert(Y.nwords > 0, "Gfmat::operator*(): dimension mismatch");

  /*
  // this can be done more efficiently?
  GF2mat result(X.nrows,Y.ncols);
  for (int i=0; i<X.nrows; i++) {
  for (int j=0; j<Y.ncols; j++) {
  bin b=0;
  for (int k=0; k<X.ncols; k++) {
  bin x = X.get(i,k);
  bin y = Y.get(k,j);
  b ^= (x&y);
  }
  result.set(i,j,b);
  }
  }
  return result; */

  // is this better?
  return mult_trans(X, Y.transpose());
}

bvec operator*(const GF2mat &X, const bvec &y)
{
  it_assert(length(y) == X.ncols, "GF2mat::operator*(): dimension mismatch");
  it_assert(X.nwords > 0, "Gfmat::operator*(): dimension mismatch");

  /*
  // this can be done more efficiently?
  bvec result = zeros_b(X.nrows);
  for (int i=0; i<X.nrows; i++) {
  // do the nwords-1 data columns first
  for (int j=0; j<X.nwords-1; j++) {
  int ind = j<<shift_divisor;
  unsigned char r=X.data(i,j);
  while (r) {
  result(i) ^= (r & y(ind));
  r >>= 1;
  ind++;
  }
  }
  // do the last column separately
  for (int j=(X.nwords-1)<<shift_divisor; j<X.ncols; j++) {
  result(i) ^= (X.get(i,j) & y(j));
  }
  }
  return result; */

  // is this better?
  return (mult_trans(X, GF2mat(y, 0))).bvecify();
}

GF2mat mult_trans(const GF2mat &X, const GF2mat &Y)
{
  it_assert(X.ncols == Y.ncols, "GF2mat::mult_trans(): dimension mismatch");
  it_assert(X.nwords > 0, "GF2mat::mult_trans(): dimension mismatch");
  it_assert(Y.nwords > 0, "GF2mat::mult_trans(): dimension mismatch");
  it_assert(X.nwords == Y.nwords, "GF2mat::mult_trans(): dimension mismatch");

  GF2mat result(X.nrows, Y.nrows);

  for (int i = 0; i < X.nrows; i++) {
    for (int j = 0; j < Y.nrows; j++) {
      bin b = 0;
      int kindx = i;
      int kindy = j;
      for (int k = 0; k < X.nwords; k++) {
        unsigned char r = X.data(kindx) & Y.data(kindy);
        /* The following can be speeded up by using a small lookup
           table for the number of ones and shift only a few times (or
           not at all if table is large) */
        while (r) {
          b ^= r & 1;
          r >>= 1;
        };
        kindx += X.nrows;
        kindy += Y.nrows;
      }
      result.set(i, j, b);
    }
  }
  return result;
}

GF2mat GF2mat::transpose() const
{
  // CAN BE SPEEDED UP
  GF2mat result(ncols, nrows);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      result.set(j, i, get(i, j));
    }
  }
  return result;
}

GF2mat operator+(const GF2mat &X, const GF2mat &Y)
{
  it_assert(X.nrows == Y.nrows, "GF2mat::operator+(): dimension mismatch");
  it_assert(X.ncols == Y.ncols, "GF2mat::operator+(): dimension mismatch");
  it_assert(X.nwords == Y.nwords, "GF2mat::operator+(): dimension mismatch");
  GF2mat result(X.nrows, X.ncols);

  for (int i = 0; i < X.nrows; i++) {
    for (int j = 0; j < X.nwords; j++) {
      result.data(i, j) = X.data(i, j) ^ Y.data(i, j);
    }
  }

  return result;
}

void GF2mat::permute_cols(ivec &perm, bool I)
{
  it_assert(length(perm) == ncols,
            "GF2mat::permute_cols(): dimensions do not match");

  GF2mat temp = (*this);
  for (int j = 0; j < ncols; j++) {
    if (I == 0) {
      set_col(j, temp.get_col(perm(j)));
    }
    else {
      set_col(perm(j), temp.get_col(j));
    }
  }
}

void GF2mat::permute_rows(ivec &perm, bool I)
{
  it_assert(length(perm) == nrows,
            "GF2mat::permute_rows(): dimensions do not match");

  GF2mat temp = (*this);
  for (int i = 0; i < nrows; i++) {
    if (I == 0) {
      for (int j = 0; j < nwords; j++) {
        data(i, j) = temp.data(perm(i), j);
      }
    }
    else {
      for (int j = 0; j < nwords; j++) {
        data(perm(i), j) = temp.data(i, j);
      }
    }
  }
}


std::ostream &operator<<(std::ostream &os, const GF2mat &X)
{
  int i, j;
  os << "---- GF(2) matrix of dimension " << X.nrows << "*" << X.ncols
  << " -- Density: " << X.density() << " ----" << std::endl;

  for (i = 0; i < X.nrows; i++) {
    os << "      ";
    for (j = 0; j < X.ncols; j++) {
      os << X.get(i, j) << " ";
    }
    os << std::endl;
  }

  return os;
}

double GF2mat::density() const
{
  int no_of_ones = 0;

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      no_of_ones += (get(i, j) == 1 ? 1 : 0);
    }
  }

  return ((double) no_of_ones) / (nrows*ncols);
}


it_file &operator<<(it_file &f, const GF2mat &X)
{
  // 3 64-bit unsigned words for: nrows, ncols and nwords + rest for char data
  uint64_t bytecount = 3 * sizeof(uint64_t)
                       + X.nrows * X.nwords * sizeof(char);
  f.write_data_header("GF2mat", bytecount);

  f.low_level_write(static_cast<uint64_t>(X.nrows));
  f.low_level_write(static_cast<uint64_t>(X.ncols));
  f.low_level_write(static_cast<uint64_t>(X.nwords));
  for (int i = 0; i < X.nrows; i++) {
    for (int j = 0; j < X.nwords; j++) {
      f.low_level_write(static_cast<char>(X.data(i, j)));
    }
  }
  return f;
}

it_ifile &operator>>(it_ifile &f, GF2mat &X)
{
  it_file::data_header h;

  f.read_data_header(h);
  if (h.type == "GF2mat") {
    uint64_t tmp;
    f.low_level_read(tmp);
    X.nrows = static_cast<int>(tmp);
    f.low_level_read(tmp);
    X.ncols = static_cast<int>(tmp);
    f.low_level_read(tmp);
    X.nwords = static_cast<int>(tmp);
    X.data.set_size(X.nrows, X.nwords);
    for (int i = 0; i < X.nrows; i++) {
      for (int j = 0; j < X.nwords; j++) {
        char r;
        f.low_level_read(r);
        X.data(i, j) = static_cast<unsigned char>(r);
      }
    }
  }
  else {
    it_error("it_ifile &operator>>() - internal error");
  }

  return f;
}

} // namespace itpp

