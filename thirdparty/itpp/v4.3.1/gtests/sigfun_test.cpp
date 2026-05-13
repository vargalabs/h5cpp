/*!
 * \file
 * \brief Signal processing routines test program
 * \author Tony Ottosson, Thomas Eriksson, Pal Frenger, Tobias Ringstrom
 *         and Adam Piatyszek
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

#include <itpp/itsignal.h>
#include "gtest/gtest.h"

using namespace itpp;

static
void assert_vec_p(const vec &ref, const vec &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol);
  }
}
#define assert_vec(ref, act) assert_vec_p(ref, act, __LINE__)

TEST(SigFun, All)
{
  RNG_reset(12345);

  vec x = randn(10);
  vec y = randn(10);
  int max_lag;

  // Testing cross correlation
  vec ref = "0.343673 2.34251 -0.276029 -1.69426 0.575064 0.236971 -1.39007 -2.51791 3.36735 3.37605 -0.551727 "
  "-1.91208 1.41728 1.65607 -3.27706 -0.802629 1.02004 0.841397 0.149574";
  assert_vec(ref, xcorr(x, y));
  max_lag = 0;
  ref = "3.37605";
  assert_vec(ref, xcorr(x, y, max_lag));
  max_lag = 5;
  ref = "0.575064 0.236971 -1.39007 -2.51791 3.36735 3.37605 -0.551727 -1.91208 1.41728 1.65607 -3.27706";
  assert_vec(ref, xcorr(x, y, max_lag));
  max_lag = -1;
  ref = "0.343673 2.34251 -0.276029 -1.69426 0.575064 0.236971 -1.39007 -2.51791 3.36735 3.37605 -0.551727 "
  "-1.91208 1.41728 1.65607 -3.27706 -0.802629 1.02004 0.841397 0.149574";
  assert_vec(ref, xcorr(x, y, max_lag));

  max_lag = 0;
  ref = "0.337605";
  assert_vec(ref, xcorr(x, y, max_lag, "biased"));
  max_lag = 5;
  ref = "0.0575064 0.0236971 -0.139007 -0.251791 0.336735 0.337605 -0.0551727 -0.191208 0.141728 0.165607 -0.327706";
  assert_vec(ref, xcorr(x, y, max_lag, "biased"));
  max_lag = -1;
  ref = "0.0343673 0.234251 -0.0276029 -0.169426 0.0575064 0.0236971 -0.139007 -0.251791 0.336735 0.337605 -0.0551727 "
  "-0.191208 0.141728 0.165607 -0.327706 -0.0802629 0.102004 0.0841397 0.0149574";
  assert_vec(ref, xcorr(x, y, max_lag, "biased"));

  max_lag = 0;
  ref = "0.337605";
  assert_vec(ref, xcorr(x, y, max_lag, "unbiased"));
  max_lag = 5;
  ref = "0.115013 0.0394951 -0.198582 -0.314739 0.37415 0.337605 -0.061303 -0.23901 0.202469 0.276012 -0.655413";
  assert_vec(ref, xcorr(x, y, max_lag, "unbiased"));
  max_lag = -1;
  ref = "0.343673 1.17126 -0.0920096 -0.423564 0.115013 0.0394951 -0.198582 -0.314739 0.37415 0.337605 -0.061303 "
  "-0.23901 0.202469 0.276012 -0.655413 -0.200657 0.340014 0.420698 0.149574";
  assert_vec(ref, xcorr(x, y, max_lag, "unbiased"));

  max_lag = 0;
  ref = "0.453329";
  assert_vec(ref, xcorr(x, y, max_lag, "coeff"));
  max_lag = 5;
  ref = "0.0772185 0.03182 -0.186656 -0.3381 0.452161 0.453329 -0.0740849 -0.25675 0.19031 0.222374 -0.440038";
  assert_vec(ref, xcorr(x, y, max_lag, "coeff"));
  max_lag = -1;
  ref = "0.0461477 0.314548 -0.0370646 -0.227502 0.0772185 0.03182 -0.186656 -0.3381 0.452161 0.453329 -0.0740849 "
  "-0.25675 0.19031 0.222374 -0.440038 -0.107776 0.136969 0.112981 0.0200845";
  assert_vec(ref, xcorr(x, y, max_lag, "coeff"));

  // Testing auto correlation
  ref = "0.170114 1.41078 1.78954 -0.36962 -5.83 3.57721 -2.23494 -1.06105 -3.00721 14.5 -3.00721 -1.06105 -2.23494 "
  "3.57721 -5.83 -0.36962 1.78954 1.41078 0.170114";
  assert_vec(ref, xcorr(x));
  max_lag = 0;
  ref = "14.5";
  assert_vec(ref, xcorr(x, max_lag));
  max_lag = 5;
  ref = "-5.83 3.57721 -2.23494 -1.06105 -3.00721 14.5 -3.00721 -1.06105 -2.23494 3.57721 -5.83";
  assert_vec(ref, xcorr(x, max_lag));
  max_lag = -1;
  ref = "0.170114 1.41078 1.78954 -0.36962 -5.83 3.57721 -2.23494 -1.06105 -3.00721 14.5 -3.00721 -1.06105 -2.23494 "
  "3.57721 -5.83 -0.36962 1.78954 1.41078 0.170114";
  assert_vec(ref, xcorr(x, max_lag));

  max_lag = 0;
  ref = "1.45";
  assert_vec(ref, xcorr(x, max_lag, "biased"));
  max_lag = 5;
  ref = "-0.583 0.357721 -0.223494 -0.106105 -0.300721 1.45 -0.300721 -0.106105 -0.223494 0.357721 -0.583";
  assert_vec(ref, xcorr(x, max_lag, "biased"));
  max_lag = -1;
  ref = "0.0170114 0.141078 0.178954 -0.036962 -0.583 0.357721 -0.223494 -0.106105 -0.300721 1.45 -0.300721 -0.106105 "
  "-0.223494 0.357721 -0.583 -0.036962 0.178954 0.141078 0.0170114";
  assert_vec(ref, xcorr(x, max_lag, "biased"));

  max_lag = 0;
  ref = "1.45";
  assert_vec(ref, xcorr(x, max_lag, "unbiased"));
  max_lag = 5;
  ref = "-1.166 0.596201 -0.319278 -0.132632 -0.334135 1.45 -0.334135 -0.132632 -0.319278 0.596201 -1.166";
  assert_vec(ref, xcorr(x, max_lag, "unbiased"));
  max_lag = -1;
  ref = "0.170114 0.70539 0.596515 -0.092405 -1.166 0.596201 -0.319278 -0.132632 -0.334135 1.45 -0.334135 -0.132632 "
  "-0.319278 0.596201 -1.166 -0.092405 0.596515 0.70539 0.170114";
  assert_vec(ref, xcorr(x, max_lag, "unbiased"));

  max_lag = 0;
  ref = "1";
  assert_vec(ref, xcorr(x, max_lag, "coeff"));
  max_lag = 5;
  ref = "-0.402069 0.246704 -0.154134 -0.0731761 -0.207394 1 -0.207394 -0.0731761 -0.154134 0.246704 -0.402069";
  assert_vec(ref, xcorr(x, max_lag, "coeff"));
  max_lag = -1;
  ref = "0.011732 0.0972952 0.123417 -0.025491 -0.402069 0.246704 -0.154134 -0.0731761 -0.207394 1 -0.207394 "
  "-0.0731761 -0.154134 0.246704 -0.402069 -0.025491 0.123417 0.0972952 0.011732";
  assert_vec(ref, xcorr(x, max_lag, "coeff"));
}
