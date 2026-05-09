/*!
 * \file
 * \brief diagonal Mixture of Gaussians class - source file
 * \author Conrad Sanderson
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

#include <itpp/base/math/log_exp.h>
#include <itpp/stat/mog_diag.h>
#include <cstdlib>


namespace itpp
{

double MOG_diag::log_lhood_single_gaus_internal(const double * c_x_in, const int k) const
{

  const double * c_mean = c_means[k];
  const double * c_diag_cov_inv_etc = c_diag_covs_inv_etc[k];

  double acc = 0.0;

  for (int d = 0; d < D; d++) {
    double tmp_val = c_x_in[d] - c_mean[d];
    acc += (tmp_val * tmp_val) * c_diag_cov_inv_etc[d];
  }
  return(c_log_det_etc[k] - acc);
}


double MOG_diag::log_lhood_single_gaus_internal(const vec &x_in, const int k) const
{
  return log_lhood_single_gaus_internal(x_in._data(), k);
}


double MOG_diag::log_lhood_single_gaus(const double * c_x_in, const int k) const
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::log_lhood_single_gaus(): model not valid");
    it_assert(((k >= 0) && (k < K)), "MOG::log_lhood_single_gaus(): k specifies a non-existant Gaussian");
  }
  return log_lhood_single_gaus_internal(c_x_in, k);
}


double MOG_diag::log_lhood_single_gaus(const vec &x_in, const int k) const
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::log_lhood_single_gaus(): model not valid");
    it_assert(check_size(x_in), "MOG_diag::log_lhood_single_gaus(): x has wrong dimensionality");
    it_assert(((k >= 0) && (k < K)), "MOG::log_lhood_single_gaus(): k specifies a non-existant Gaussian");
  }
  return log_lhood_single_gaus_internal(x_in._data(), k);
}


double MOG_diag::log_lhood_internal(const double * c_x_in)
{

  bool danger = paranoid;

  for (int k = 0;k < K;k++)  {
    double tmp = c_log_weights[k] + log_lhood_single_gaus_internal(c_x_in, k);
    c_tmpvecK[k] = tmp;

    if (tmp >= log_max_K)  danger = true;
  }


  if (danger) {
    double log_sum = c_tmpvecK[0];
    for (int k = 1; k < K; k++)  log_sum = log_add(log_sum, c_tmpvecK[k]);
    return(log_sum);
  }
  else {
    double sum = 0.0;
    for (int k = 0;k < K;k++) sum += std::exp(c_tmpvecK[k]);
    return(std::log(sum));
  }
}


double MOG_diag::log_lhood_internal(const vec &x_in)
{
  return log_lhood_internal(x_in._data());
}


double MOG_diag::log_lhood(const vec &x_in)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::log_lhood(): model not valid");
    it_assert(check_size(x_in), "MOG_diag::log_lhood(): x has wrong dimensionality");
  }
  return log_lhood_internal(x_in._data());
}


double MOG_diag::log_lhood(const double * c_x_in)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::log_lhood(): model not valid");
    it_assert((c_x_in != 0), "MOG_diag::log_lhood(): c_x_in is a null pointer");
  }

  return log_lhood_internal(c_x_in);
}


double MOG_diag::lhood_internal(const double * c_x_in)
{

  bool danger = paranoid;

  for (int k = 0;k < K;k++)  {
    double tmp = c_log_weights[k] + log_lhood_single_gaus_internal(c_x_in, k);
    c_tmpvecK[k] = tmp;

    if (tmp >= log_max_K)  danger = true;
  }


  if (danger) {
    double log_sum = c_tmpvecK[0];
    for (int k = 1; k < K; k++)  log_sum = log_add(log_sum, c_tmpvecK[k]);
    return(trunc_exp(log_sum));
  }
  else {
    double sum = 0.0;
    for (int k = 0;k < K;k++) sum += std::exp(c_tmpvecK[k]);
    return(sum);
  }
}

double MOG_diag::lhood_internal(const vec &x_in) { return lhood_internal(x_in._data()); }

double MOG_diag::lhood(const vec &x_in)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::lhood(): model not valid");
    it_assert(check_size(x_in), "MOG_diag::lhood(): x has wrong dimensionality");
  }
  return lhood_internal(x_in._data());
}


double MOG_diag::lhood(const double * c_x_in)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::lhood(): model not valid");
    it_assert((c_x_in != 0), "MOG_diag::lhood(): c_x_in is a null pointer");
  }

  return lhood_internal(c_x_in);
}


double MOG_diag::avg_log_lhood(const double ** c_x_in, const int N)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::avg_log_lhood(): model not valid");
    it_assert((c_x_in != 0), "MOG_diag::avg_log_lhood(): c_x_in is a null pointer");
    it_assert((N >= 0), "MOG_diag::avg_log_lhood(): N is zero or negative");
  }

  double acc = 0.0;
  for (int n = 0;n < N;n++) acc += log_lhood_internal(c_x_in[n]);
  return(acc / N);
}


double MOG_diag::avg_log_lhood(const Array<vec> &X_in)
{
  if (do_checks) {
    it_assert(valid, "MOG_diag::avg_log_lhood(): model not valid");
    it_assert(check_size(X_in), "MOG_diag::avg_log_lhood(): X is empty or at least one vector has the wrong dimensionality");
  }
  const int N = X_in.size();
  double acc = 0.0;
  for (int n = 0;n < N;n++)  acc += log_lhood_internal(X_in(n)._data());
  return(acc / N);
}

void MOG_diag::zero_all_ptrs()
{
  c_means             = 0;
  c_diag_covs         = 0;
  c_diag_covs_inv_etc = 0;
  c_weights           = 0;
  c_log_weights       = 0;
  c_log_det_etc       = 0;
  c_tmpvecK           = 0;
}


void MOG_diag::free_all_ptrs()
{
  c_means             = disable_c_access(c_means);
  c_diag_covs         = disable_c_access(c_diag_covs);
  c_diag_covs_inv_etc = disable_c_access(c_diag_covs_inv_etc);
  c_weights           = disable_c_access(c_weights);
  c_log_weights       = disable_c_access(c_log_weights);
  c_log_det_etc       = disable_c_access(c_log_det_etc);
  c_tmpvecK           = disable_c_access(c_tmpvecK);
}


void MOG_diag::setup_means()
{
  MOG_generic::setup_means();
  disable_c_access(c_means);
  c_means = enable_c_access(means);
}


void MOG_diag::setup_covs()
{
  MOG_generic::setup_covs();
  if (full) return;

  disable_c_access(c_diag_covs);
  disable_c_access(c_diag_covs_inv_etc);
  disable_c_access(c_log_det_etc);

  c_diag_covs         = enable_c_access(diag_covs);
  c_diag_covs_inv_etc = enable_c_access(diag_covs_inv_etc);
  c_log_det_etc       = enable_c_access(log_det_etc);
}


void MOG_diag::setup_weights()
{
  MOG_generic::setup_weights();

  disable_c_access(c_weights);
  disable_c_access(c_log_weights);

  c_weights = enable_c_access(weights);
  c_log_weights = enable_c_access(log_weights);
}


void MOG_diag::setup_misc()
{
  disable_c_access(c_tmpvecK);
  tmpvecK.set_size(K);
  c_tmpvecK = enable_c_access(tmpvecK);

  MOG_generic::setup_misc();
  if (full) convert_to_diag_internal();
}


void MOG_diag::load(const std::string &name_in)
{
  MOG_generic::load(name_in);
  if (full) convert_to_diag();
}


double ** MOG_diag::enable_c_access(Array<vec> & A_in)
{
  int rows = A_in.size();
  double ** A = (double **)std::malloc(rows * sizeof(double *));
  if (A)  for (int row = 0;row < rows;row++)  A[row] = A_in(row)._data();
  return(A);
}

int ** MOG_diag::enable_c_access(Array<ivec> & A_in)
{
  int rows = A_in.size();
  int ** A = (int **)std::malloc(rows * sizeof(int *));
  if (A)  for (int row = 0;row < rows;row++)  A[row] = A_in(row)._data();
  return(A);
}

double ** MOG_diag::disable_c_access(double ** A_in) { if (A_in) std::free(A_in); return(0); }
int ** MOG_diag::disable_c_access(int ** A_in) { if (A_in) std::free(A_in); return(0); }

double * MOG_diag::enable_c_access(vec & v_in) { return v_in._data(); }
int * MOG_diag::enable_c_access(ivec & v_in) { return v_in._data(); }

double * MOG_diag::disable_c_access(double *) { return(0); }
int * MOG_diag::disable_c_access(int *) { return(0); }

}
