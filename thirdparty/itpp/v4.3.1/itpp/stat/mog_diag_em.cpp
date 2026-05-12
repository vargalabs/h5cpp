/*!
 * \file
 * \brief Expectation Maximisation based optimisers for Mixture of Gaussians - source file
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


#include <itpp/stat/mog_diag_em.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/timing.h>

#include <iostream>
#include <iomanip>

namespace itpp
{

//! update log versions of parameters and any necessary constants
void inline MOG_diag_EM_sup::update_internals()
{

  double Ddiv2_log_2pi = D / 2.0 * std::log(m_2pi);

  for (int k = 0;k < K;k++)  c_log_weights[k] = std::log(c_weights[k]);

  for (int k = 0;k < K;k++) {
    double acc = 0.0;
    double * c_diag_cov = c_diag_covs[k];
    double * c_diag_cov_inv_etc = c_diag_covs_inv_etc[k];

    for (int d = 0;d < D;d++) {
      double tmp = c_diag_cov[d];
      c_diag_cov_inv_etc[d] = 1.0 / (2.0 * tmp);
      acc += std::log(tmp);
    }

    c_log_det_etc[k] = -Ddiv2_log_2pi - 0.5 * acc;
  }

}


//! for helping to avoid numerical instability
void inline MOG_diag_EM_sup::sanitise_params()
{

  double acc = 0.0;
  for (int k = 0;k < K;k++) {
    if (c_weights[k] < weight_floor)  c_weights[k] = weight_floor;
    if (c_weights[k] > 1.0)  c_weights[k] = 1.0;
    acc += c_weights[k];
  }
  for (int k = 0;k < K;k++)  c_weights[k] /= acc;

  for (int k = 0;k < K;k++)
    for (int d = 0;d < D;d++)
      if (c_diag_covs[k][d] < var_floor)  c_diag_covs[k][d] = var_floor;

}

//! update parameters using the Maximum Likelihood version of the EM algorithm
double MOG_diag_EM_sup::ml_update_params()
{

  double acc_loglhood = 0.0;

  for (int k = 0;k < K;k++)  {
    c_acc_loglhood_K[k] = 0.0;

    double * c_acc_mean = c_acc_means[k];
    double * c_acc_cov  = c_acc_covs[k];

    for (int d = 0;d < D;d++) { c_acc_mean[d] = 0.0; c_acc_cov[d] = 0.0; }
  }

  for (int n = 0;n < N;n++) {
    double * c_x =  c_X[n];

    bool danger = paranoid;
    for (int k = 0;k < K;k++)  {
      double tmp = c_log_weights[k] + MOG_diag::log_lhood_single_gaus_internal(c_x, k);
      c_tmpvecK[k] = tmp;
      if (tmp >= log_max_K)  danger = true;
    }

    if (danger) {

      double log_sum = c_tmpvecK[0];
      for (int k = 1;k < K;k++)  log_sum = log_add(log_sum, c_tmpvecK[k]);
      acc_loglhood += log_sum;

      for (int k = 0;k < K;k++) {

        double * c_acc_mean = c_acc_means[k];
        double * c_acc_cov = c_acc_covs[k];

        double tmp_k = trunc_exp(c_tmpvecK[k] - log_sum);
        acc_loglhood_K[k] += tmp_k;

        for (int d = 0;d < D;d++) {
          double tmp_x = c_x[d];
          c_acc_mean[d] +=  tmp_k * tmp_x;
          c_acc_cov[d] += tmp_k * tmp_x * tmp_x;
        }
      }
    }
    else {

      double sum = 0.0;
      for (int k = 0;k < K;k++) { double tmp = std::exp(c_tmpvecK[k]); c_tmpvecK[k] = tmp; sum += tmp; }
      acc_loglhood += std::log(sum);

      for (int k = 0;k < K;k++) {

        double * c_acc_mean = c_acc_means[k];
        double * c_acc_cov = c_acc_covs[k];

        double tmp_k = c_tmpvecK[k] / sum;
        c_acc_loglhood_K[k] += tmp_k;

        for (int d = 0;d < D;d++) {
          double tmp_x = c_x[d];
          c_acc_mean[d] +=  tmp_k * tmp_x;
          c_acc_cov[d] += tmp_k * tmp_x * tmp_x;
        }
      }
    }
  }

  for (int k = 0;k < K;k++) {

    double * c_mean = c_means[k];
    double * c_diag_cov = c_diag_covs[k];

    double * c_acc_mean = c_acc_means[k];
    double * c_acc_cov = c_acc_covs[k];

    double tmp_k = c_acc_loglhood_K[k];

    c_weights[k] = tmp_k / N;

    for (int d = 0;d < D;d++) {
      double tmp_mean = c_acc_mean[d] / tmp_k;
      c_mean[d] = tmp_mean;
      c_diag_cov[d] = c_acc_cov[d] / tmp_k - tmp_mean * tmp_mean;
    }
  }

  return(acc_loglhood / N);

}


void MOG_diag_EM_sup::ml_iterate()
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::showpos;
  using std::noshowpos;
  using std::scientific;
  using std::fixed;
  using std::flush;
  using std::setprecision;

  double avg_log_lhood_old = -1.0 * std::numeric_limits<double>::max();

  Real_Timer tt;

  if (verbose) {
    cout << "MOG_diag_EM_sup::ml_iterate()" << endl;
    cout << setw(14) << "iteration";
    cout << setw(14) << "avg_loglhood";
    cout << setw(14) << "delta";
    cout << setw(10) << "toc";
    cout << endl;
  }

  for (int i = 0; i < max_iter; i++) {
    sanitise_params();
    update_internals();

    if (verbose) tt.tic();
    double avg_log_lhood_new = ml_update_params();

    if (verbose) {
      double delta = avg_log_lhood_new - avg_log_lhood_old;

      cout << noshowpos << fixed;
      cout << setw(14) << i;
      cout << showpos << scientific << setprecision(3);
      cout << setw(14) << avg_log_lhood_new;
      cout << setw(14) << delta;
      cout << noshowpos << fixed;
      cout << setw(10) << tt.toc();
      cout << endl << flush;
    }

    if (avg_log_lhood_new <= avg_log_lhood_old)  break;

    avg_log_lhood_old = avg_log_lhood_new;
  }
}


void MOG_diag_EM_sup::ml(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in, double var_floor_in, double weight_floor_in, bool verbose_in)
{

  it_assert(model_in.is_valid(), "MOG_diag_EM_sup::ml(): initial model not valid");
  it_assert(check_array_uniformity(X_in), "MOG_diag_EM_sup::ml(): 'X' is empty or contains vectors of varying dimensionality");
  it_assert((max_iter_in > 0), "MOG_diag_EM_sup::ml(): 'max_iter' needs to be greater than zero");

  verbose = verbose_in;

  N = X_in.size();

  Array<vec> means_in = model_in.get_means();
  Array<vec> diag_covs_in = model_in.get_diag_covs();
  vec weights_in = model_in.get_weights();

  init(means_in, diag_covs_in, weights_in);

  means_in.set_size(0);
  diag_covs_in.set_size(0);
  weights_in.set_size(0);

  if (K > N) {
    it_warning("MOG_diag_EM_sup::ml(): WARNING: K > N");
  }
  else {
    if (K > N / 10) {
      it_warning("MOG_diag_EM_sup::ml(): WARNING: K > N/10");
    }
  }

  var_floor = var_floor_in;
  weight_floor = weight_floor_in;

  const double tiny = std::numeric_limits<double>::min();
  if (var_floor < tiny) var_floor = tiny;
  if (weight_floor < tiny) weight_floor = tiny;
  if (weight_floor > 1.0 / K) weight_floor = 1.0 / K;

  max_iter = max_iter_in;

  tmpvecK.set_size(K);
  tmpvecD.set_size(D);
  acc_loglhood_K.set_size(K);

  acc_means.set_size(K);
  for (int k = 0;k < K;k++) acc_means(k).set_size(D);
  acc_covs.set_size(K);
  for (int k = 0;k < K;k++) acc_covs(k).set_size(D);

  c_X = enable_c_access(X_in);
  c_tmpvecK = enable_c_access(tmpvecK);
  c_tmpvecD = enable_c_access(tmpvecD);
  c_acc_loglhood_K = enable_c_access(acc_loglhood_K);
  c_acc_means = enable_c_access(acc_means);
  c_acc_covs = enable_c_access(acc_covs);

  ml_iterate();

  model_in.init(means, diag_covs, weights);

  disable_c_access(c_X);
  disable_c_access(c_tmpvecK);
  disable_c_access(c_tmpvecD);
  disable_c_access(c_acc_loglhood_K);
  disable_c_access(c_acc_means);
  disable_c_access(c_acc_covs);


  tmpvecK.set_size(0);
  tmpvecD.set_size(0);
  acc_loglhood_K.set_size(0);
  acc_means.set_size(0);
  acc_covs.set_size(0);

  cleanup();

}

void MOG_diag_EM_sup::map(MOG_diag &, MOG_diag &, Array<vec> &, int, double,
                          double, double, bool)
{
  it_error("MOG_diag_EM_sup::map(): not implemented yet");
}


//
// convenience functions

void MOG_diag_ML(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in, double var_floor_in, double weight_floor_in, bool verbose_in)
{
  MOG_diag_EM_sup EM;
  EM.ml(model_in, X_in, max_iter_in, var_floor_in, weight_floor_in, verbose_in);
}

void MOG_diag_MAP(MOG_diag &, MOG_diag &, Array<vec> &, int, double, double,
                  double, bool)
{
  it_error("MOG_diag_MAP(): not implemented yet");
}

}

