/*!
 * \file
 * \brief Matrix inversion routines test program
 * \author Tony Ottosson and Adam Piatyszek
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
using namespace std;

static
void assert_mat(const mat &exp, const mat &act)
{
  static const double eps = 1e-5;
  ASSERT_EQ(exp.rows(), act.rows());
  ASSERT_EQ(exp.cols(), act.cols());
  for (int n = 0; n < exp.rows(); ++n)
  {
    for (int k = 0; k < exp.cols(); ++k)
    {
      ASSERT_NEAR(exp(n,k), act(n,k), eps);
    }
  }
}
static
void assert_cmat(const cmat &exp, const cmat &act)
{
  static const double eps = 1e-5;
  ASSERT_EQ(exp.rows(), act.rows());
  ASSERT_EQ(exp.cols(), act.cols());
  for (int n = 0; n < exp.rows(); ++n)
  {
    for (int k = 0; k < exp.cols(); ++k)
    {
      ASSERT_NEAR(exp(n,k).real(), act(n,k).real(), eps);
      ASSERT_NEAR(exp(n,k).imag(), act(n,k).imag(), eps);
    }
  }
}

TEST (Inv, All)
{
  // Test of Matrix inversion routines
  RNG_reset(0);
  {
    // Real matrix

    mat X = randn(5, 5), Y;
    Y = inv(X);
    mat Y_ref = "1.02914 0.086169 2.46109 -0.54101 2.55964;"
                "6.03751 0.193059 7.44752 -1.88403 12.0398;"
                "0.636496 -0.231074 0.457221 -0.233538 1.23596;"
                "-4.4763 -0.0890091 -4.95943 1.47648 -7.56113;"
                "3.87823 -0.0536348 4.65722 -0.486929 7.74427";
    assert_mat(Y_ref, Y);

    X = randn(5, 5);
    Y = inv(X);
    Y_ref = "3.18059 -3.97463 2.92744 0.223873 2.60524;"
            "2.31265 -2.53672 1.42136 0.543372 1.29475;"
            "-1.10003 0.582699 -0.380368 -0.218256 -0.819337;"
            "-1.18654 0.94806 -0.743228 0.142985 -0.664563;"
            "-2.37153 3.2677 -1.79564 -0.228335 -1.83916";
    assert_mat(Y_ref, Y);
  }
  {
    // Complex matrix

    cmat X = randn_c(5, 5), Y;
    Y = inv(X);
    cmat Y_ref = "-0.233219-0.564957i 0.289735+0.149868i 0.470503+0.172316i 0.461592+0.104654i -0.341534-0.150273i;"
                 "-0.458316-0.459834i -0.0783873+0.212039i 0.266618-0.511581i 0.0741277-0.555698i -0.325349+0.120556i;"
                 "-0.199788-0.443033i 0.168852+0.46391i 0.881415-0.198862i 0.891107-0.532156i -1.02218+0.518232i;"
                 "-0.513695-0.0187012i 0.357327-0.185935i 0.0794018-0.0419093i -0.240815-0.251631i -0.301137+0.182979i;"
                 "-0.823673+0.556726i 0.116893+0.0113572i -1.17761-0.221457i -1.08623-0.382027i 1.03425+0.763366i";
    assert_cmat(Y_ref, Y);

    X = randn_c(5, 5);
    Y = inv(X);
    Y_ref = "-0.612452-1.18464i -0.491647-1.60476i 0.242101+0.273005i 0.354203+0.915628i -0.775818+0.913778i;"
            "-0.0387323-0.968675i -0.0942719-0.428252i -0.0993426+0.234267i -0.153272+0.210774i -0.504668-0.046527i;"
            "-0.78825+0.122621i -0.145274+0.16973i 0.264336+0.14154i 0.142947-0.44944i -0.297705+0.396744i;"
            "-0.998674+0.323841i -1.24939+1.84618i -0.364234+0.311316i 0.681211-1.20705i 0.512001+0.157226i;"
            "-0.729248+1.85844i -0.952761+2.86892i 0.410097+0.0181813i 0.0702624-1.47188i 1.2297+0.186874i";
    assert_cmat(Y_ref, Y);
  }
}
