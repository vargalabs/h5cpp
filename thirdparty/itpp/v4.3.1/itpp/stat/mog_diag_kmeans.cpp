/*!
 * \file
 * \brief kmeans based optimiser for Mixture of Gaussians - source file
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


#include <itpp/stat/mog_diag_kmeans.h>
#include <iostream>

namespace itpp
{

inline double MOG_diag_kmeans_sup::dist(const double * x, const double * y) const
{
  double acc = 0.0;
  for (int d = 0;d < D;d++) { double tmp = x[d] - y[d]; acc += tmp * tmp; }
  return(acc);
}

void MOG_diag_kmeans_sup::assign_to_means()
{

  for (int k = 0;k < K;k++) c_count[k] = 0;

  for (int n = 0;n < N;n++) {

    int k = 0;
    double min_dist = dist(c_means[k], c_X[n]);
    int k_winner = k;

    for (int k = 1;k < K;k++) {
      double tmp_dist = dist(c_means[k], c_X[n]);
      if (tmp_dist < min_dist) { min_dist = tmp_dist; k_winner = k; }
    }

    c_partitions[ k_winner ][ count[k_winner] ] = n;
    c_count[k_winner]++;
  }
}


void MOG_diag_kmeans_sup::recalculate_means()
{

  for (int k = 0;k < K;k++) {

    for (int d = 0;d < D;d++)  c_tmpvec[d] = 0.0;

    int Nk = c_count[k];

    for (int n = 0;n < Nk;n++) {
      double * x = c_X[ c_partitions[k][n] ];
      for (int d = 0;d < D;d++)  c_tmpvec[d] += x[d];
    }

    if (Nk > 0) {
      double * c_mean = c_means[k];
      for (int d = 0;d < D;d++)  c_mean[d] = c_tmpvec[d] / Nk;
    }
  }

}


bool MOG_diag_kmeans_sup::dezombify_means()
{

  static int counter = 0;

  bool zombie_mean = false;

  int k = 0;
  int max_count = count[k];
  int k_hog = k;

  for (int k = 1;k < K;k++)  if (c_count[k] > max_count) { max_count = c_count[k]; k_hog = k; }

  for (int k = 0;k < K;k++) {
    if (c_count[k] == 0) {

      zombie_mean = true;
      if (verbose) {
        it_warning("MOG_diag_kmeans_sup::dezombify_means(): detected zombie mean");
      }
      if (k_hog == k) {
        it_warning("MOG_diag_kmeans_sup::dezombify_means(): weirdness: k_hog == k");
        return(false);
      }

      if (counter >= c_count[k_hog])  counter = 0;

      double * c_mean = c_means[k];
      double * c_x = c_X[ c_partitions[k_hog][counter] ];

      for (int d = 0;d < D;d++)  c_mean[d] = 0.5 * (c_means[k_hog][d] + c_x[d]);
      counter++;
    }

  }

  if (zombie_mean)  assign_to_means();

  return(true);
}


double MOG_diag_kmeans_sup::measure_change() const
{

  double tmp_dist = 0.0;
  for (int k = 0;k < K;k++)   tmp_dist += dist(c_means[k], c_means_old[k]);
  return(tmp_dist);
}


void MOG_diag_kmeans_sup::initial_means()
{

  for (int d = 0;d < D;d++) c_tmpvec[d] = 0.0;

  for (int n = 0;n < N;n++) {
    double * c_x = c_X[n];
    for (int d = 0;d < D;d++)  c_tmpvec[d] += c_x[d];
  }

  for (int d = 0;d < D;d++) c_tmpvec[d] /= N;

  int step = int(floor(double(N) / double(K)));
  for (int k = 0;k < K;k++) {
    double * c_mean = c_means[k];
    double * c_x = c_X[k*step];

    for (int d = 0;d < D;d++)  c_mean[d] = 0.5 * (c_tmpvec[d] + c_x[d]);
  }
}


void MOG_diag_kmeans_sup::iterate()
{

  for (int k = 0;k < K;k++)  for (int d = 0;d < D;d++)  c_means_old[k][d] = c_means[k][d];

  for (int i = 0;i < max_iter;i++) {

    assign_to_means();
    if (!dezombify_means()) return;
    recalculate_means();

    double change = measure_change();

    if (verbose) std::cout << "MOG_diag_kmeans_sup::iterate(): iteration = " << i << "  change = " << change << std::endl;
    if (change == 0) break;

    for (int k = 0;k < K;k++)  for (int d = 0;d < D;d++)   c_means_old[k][d] = c_means[k][d];
  }

}


void MOG_diag_kmeans_sup::calc_means()
{
  initial_means();
  iterate();
}


void MOG_diag_kmeans_sup::calc_covs()
{

  for (int k = 0;k < K;k++) {
    int Nk = c_count[k];

    if (Nk >= 2) {
      double * c_mean = c_means[k];

      for (int d = 0;d < D;d++) c_tmpvec[d] = 0.0;

      for (int n = 0;n < Nk;n++) {
        double * c_x = c_X[ c_partitions[k][n] ];
        for (int d = 0;d < D;d++) { double tmp = c_x[d] - c_mean[d];  c_tmpvec[d] += tmp * tmp; }
      }

      for (int d = 0;d < D;d++)  c_diag_covs[k][d] = trust * (c_tmpvec[d] / (Nk - 1.0)) + (1.0 - trust) * (1.0);
    }
    else {
      for (int d = 0;d < D;d++)  c_diag_covs[k][d] = 1.0;
    }
  }

}


void MOG_diag_kmeans_sup::calc_weights()
{
  for (int k = 0;k < K;k++)  c_weights[k] = trust * (c_count[k] / double(N)) + (1.0 - trust) * (1.0 / K);
}


void MOG_diag_kmeans_sup::normalise_vectors()
{

  for (int d = 0;d < D;d++) {
    double acc = 0.0;
    for (int n = 0;n < N;n++)  acc += c_X[n][d];
    c_norm_mu[d] = acc / N;
  }

  for (int d = 0;d < D;d++) {
    double acc = 0.0;
    for (int n = 0;n < N;n++) { double tmp = c_X[n][d] - c_norm_mu[d];  acc += tmp * tmp; }
    c_norm_sd[d] = std::sqrt(acc / (N - 1));
  }

  for (int n = 0;n < N;n++) for (int d = 0;d < D;d++) {
      c_X[n][d] -= c_norm_mu[d];
      if (c_norm_sd[d] > 0.0)  c_X[n][d] /= c_norm_sd[d];
    }
}


void MOG_diag_kmeans_sup::unnormalise_vectors()
{

  for (int n = 0;n < N;n++)  for (int d = 0;d < D;d++) {
      if (c_norm_sd[d] > 0.0)  c_X[n][d] *= c_norm_sd[d];
      c_X[n][d] += c_norm_mu[d];
    }
}


void MOG_diag_kmeans_sup::unnormalise_means()
{

  for (int k = 0;k < K;k++) for (int d = 0;d < D;d++) {
      if (norm_sd[d] > 0.0) c_means[k][d] *= c_norm_sd[d];
      c_means[k][d] += norm_mu[d];
    }
}


void MOG_diag_kmeans_sup::run(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in, double trust_in, bool normalise_in, bool verbose_in)
{

  it_assert(model_in.is_valid(), "MOG_diag_kmeans_sup::run(): given model is not valid");
  it_assert((max_iter_in > 0), "MOG_diag_kmeans_sup::run(): 'max_iter' needs to be greater than zero");
  it_assert(((trust_in >= 0.0) && (trust_in <= 1.0)), "MOG_diag_kmeans_sup::run(): 'trust' must be between 0 and 1 (inclusive)");

  verbose = verbose_in;

  Array<vec> means_in = model_in.get_means();
  Array<vec> diag_covs_in = model_in.get_diag_covs();
  vec weights_in = model_in.get_weights();

  init(means_in, diag_covs_in, weights_in);

  means_in.set_size(0);
  diag_covs_in.set_size(0);
  weights_in.set_size(0);

  it_assert(check_size(X_in), "MOG_diag_kmeans_sup::run(): 'X' is empty or contains vectors of wrong dimensionality");

  N = X_in.size();

  if (K > N) {
    it_warning("MOG_diag_kmeans_sup::run(): K > N");
  }
  else {
    if (K > N / 10) {
      it_warning("MOG_diag_kmeans_sup::run(): K > N/10");
    }
  }

  max_iter = max_iter_in;
  trust = trust_in;

  means_old.set_size(K);
  for (int k = 0;k < K;k++) means_old(k).set_size(D);
  partitions.set_size(K);
  for (int k = 0;k < K;k++) partitions(k).set_size(N);
  count.set_size(K);
  tmpvec.set_size(D);
  norm_mu.set_size(D);
  norm_sd.set_size(D);

  c_X = enable_c_access(X_in);
  c_means_old = enable_c_access(means_old);
  c_partitions = enable_c_access(partitions);
  c_count = enable_c_access(count);
  c_tmpvec = enable_c_access(tmpvec);
  c_norm_mu = enable_c_access(norm_mu);
  c_norm_sd = enable_c_access(norm_sd);

  if (normalise_in)  normalise_vectors();

  calc_means();
  if (normalise_in) { unnormalise_vectors(); unnormalise_means(); }
  calc_covs();
  calc_weights();

  model_in.init(means, diag_covs, weights);

  disable_c_access(c_X);
  disable_c_access(c_means_old);
  disable_c_access(c_partitions);
  disable_c_access(c_count);
  disable_c_access(c_tmpvec);
  disable_c_access(c_norm_mu);
  disable_c_access(c_norm_sd);

  means_old.set_size(0);
  partitions.set_size(0);
  count.set_size(0);
  tmpvec.set_size(0);
  norm_mu.set_size(0);
  norm_sd.set_size(0);

  cleanup();

}

//
// convenience functions

void MOG_diag_kmeans(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in, double trust_in, bool normalise_in, bool verbose_in)
{
  MOG_diag_kmeans_sup km;
  km.run(model_in, X_in, max_iter_in, trust_in, normalise_in, verbose_in);
}


}

