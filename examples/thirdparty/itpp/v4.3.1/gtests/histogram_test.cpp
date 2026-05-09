/*!
 * \file
 * \brief Histogram class test program
 * \author Andy Panov and Adam Piatyszek
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

#include <itpp/itstat.h>
#include <iomanip>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
void assert_ivec(const ivec &expected, const ivec &actual)
{
  ASSERT_EQ(expected.length(), actual.length());
  for (int n = 0; n < expected.length(); ++n)
  {
    ASSERT_EQ(expected[n], actual[n]);
  }
}
static
void assert_vec(const vec &expected, const vec &actual)
{
  static const double eps = 1e-12;
  ASSERT_EQ(expected.length(), actual.length());
  for (int n = 0; n < expected.length(); ++n)
  {
    ASSERT_NEAR(expected[n], actual[n], eps);
  }
}

TEST (Histogram, All)
{
  // Histogram tests

  // create histogram
  Histogram<double> hist(-3, 3, 21);

  // matrix dimension for statistical test
  int mat_dim = 100;

  // compute histogram for a random matrix
  RNG_reset(0);
  hist.update(randn(mat_dim, mat_dim));

  ivec bins = hist.get_bins();
  vec exp_pdf = hist.get_pdf();
  double pdf_max = max(exp_pdf);
  ASSERT_DOUBLE_EQ(0.1242, pdf_max);
  ivec bins_ref = "29 39 75 118 225 406 622 734 1035 1135 1242 1119 973 781 584 370 247 143 75 21 27";
  assert_ivec(bins_ref, bins);
  vec exp_pdf_ref = "0.0029 0.0039 0.0075 0.0118 0.0225 0.0406 0.0622 0.0734 0.1035 0.1135 0.1242 0.1119 "
  "0.0973 0.0781 0.0584 0.037 0.0247 0.0143 0.0075 0.0021 0.0027";
  assert_vec(exp_pdf_ref, exp_pdf);
  ASSERT_EQ(10000, sum(hist.get_bins()));
  ASSERT_DOUBLE_EQ(1.0, sum(exp_pdf));

  // reset histogram, so we can start next experiment
  hist.reset();

  // compute histogram for a random vector
  int num_stat_trials = 50000;

  // compute histogram for random vector
  hist.update(randn(num_stat_trials));
  bins = hist.get_bins();
  exp_pdf = hist.get_pdf();
  pdf_max = max(exp_pdf);
  ASSERT_DOUBLE_EQ(0.11752, pdf_max);
  bins_ref = "115 169 319 648 1246 2009 2917 3960 5000 5724 5876 5683 4936 4126 2880 1962 1182 679 318 157 94";
  assert_ivec(bins_ref, bins);
  exp_pdf_ref = "0.0023 0.00338 0.00638 0.01296 0.02492 0.04018 0.05834 0.0792 0.1 0.11448 0.11752 0.11366 "
  "0.09872 0.08252 0.0576 0.03924 0.02364 0.01358 0.00636 0.00314 0.00188";
  assert_vec(exp_pdf_ref, exp_pdf);
  ASSERT_EQ(50000, sum(hist.get_bins()));
  ASSERT_DOUBLE_EQ(1.0, sum(exp_pdf));

  // compute CDF. CDF is computed vs. right bin boundaries
  vec exp_cdf = hist.get_cdf();
  vec bin_right = hist.get_bin_rights();
  vec exp_cdf_ref = "0.0023 0.00568 0.01206 0.02502 0.04994 0.09012 0.14846 0.22766 0.32766 0.44214 0.55966 "
  "0.67332 0.77204 0.85456 0.91216 0.9514 0.97504 0.98862 0.99498 0.99812 1";
  vec bin_right_ref = "-2.85 -2.55 -2.25 -1.95 -1.65 -1.35 -1.05 -0.75 -0.45 -0.15 0.15 0.45 0.75 1.05 1.35 "
  "1.65 1.95 2.25 2.55 2.85 3.15";
  assert_vec(exp_cdf_ref, exp_cdf);
  assert_vec(bin_right_ref, bin_right);
}
