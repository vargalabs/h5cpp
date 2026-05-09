/*!
 * \file
 * \brief Linear systems of equations solving test program
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

static const double tol = 1e-4;

static
void assert_vec(const vec &ref, const vec &act)
{
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref[n], act[n], tol);
  }
}

static
void assert_cvec(const cvec &ref, const cvec &act)
{
  ASSERT_EQ(ref.length(), act.length());
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref[n].real(), act[n].real(), tol);
    ASSERT_NEAR(ref[n].imag(), act[n].imag(), tol);
  }
}

static
void assert_mat(const mat &ref, const mat &act)
{
  ASSERT_EQ(ref.rows(), act.rows());
  ASSERT_EQ(ref.cols(), act.cols());
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol);
    }
  }
}

static
void assert_cmat(const cmat &ref, const cmat &act)
{
  ASSERT_EQ(ref.rows(), act.rows());
  ASSERT_EQ(ref.cols(), act.cols());
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k).real(), act(n,k).real(), tol);
      ASSERT_NEAR(ref(n,k).imag(), act(n,k).imag(), tol);
    }
  }
}

TEST(LS, All)
{
  RNG_reset(0);

  // Solving linear systems of equations

  {
    // Real systems
    mat A, B, X;
    vec b, x;

    A = randn(4, 4);
    b = randn(4);
    x = ls_solve(A, b);
    vec x_ref = "-0.312176 0.00252645 -0.393292 0.493941";
    assert_vec(x_ref, x);

    A = randn(4, 4);
    B = randn(4, 2);
    X = ls_solve(A, B);
    mat X_ref = "-1.55189 -1.50623;"
            "0.185151 2.98571;"
            "-0.051986 -0.872471;"
            "0.0339548 0.426711";
    assert_mat(X_ref, X);

    A = randn(4, 4);
    A = A.transpose() * A;
    b = randn(4);
    x = ls_solve(A, b);
    x_ref = "6.96571 3.4157 -15.9777 1.67972";
    assert_vec(x_ref, x);

    A = randn(4, 4);
    A = A.transpose() * A;
    B = randn(4, 2);
    X = ls_solve(A, B);
    X_ref = "0.150362 0.405445;"
            "-3.54916 -11.7528;"
            "-0.245019 -1.77316;"
            "5.35996 14.6563";
    assert_mat(X_ref, X);

    // Overdetermined system
    A = randn(4, 2);
    b = randn(4);
    x = ls_solve_od(A, b);
    x_ref = "-0.652584 0.39863";
    assert_vec(x_ref, x);

    // Overdetermined system
    A = randn(4, 2);
    B = randn(4, 3);
    X = ls_solve_od(A, B);
    X_ref = "-1.05062 -0.44978 -0.756015;"
            "0.937086 0.249577 -0.225774";
    assert_mat(X_ref, X);

    // Underdetermined system
    A = randn(2, 4);
    b = randn(2);
    x = ls_solve_ud(A, b);
    x_ref = "1.84586 -0.741848 -1.51662 1.72224";
    assert_vec(x_ref, x);

    // Underdetermined system
    A = randn(2, 4);
    B = randn(2, 3);
    X = ls_solve_ud(A, B);
    X_ref = "-1.33908 -0.190497 -0.653759;"
            "1.84025 0.443054 0.70934;"
            "0.396185 -0.331404 0.597969;"
            "-0.442725 0.0149668 -0.297468";
    assert_mat(X_ref, X);
  }


  {
    // Complex systems
    cmat A, B, X;
    cvec b, x;

    // Square system
    A = randn_c(4, 4);
    b = randn_c(4);
    x = ls_solve(A, b);
    cvec x_ref = "-0.665037-0.657617i -0.258391+0.204449i -1.31542-0.160235i -0.490128+1.08786i";
    assert_cvec(x_ref, x);

    // Square system
    A = randn_c(4, 4);
    B = randn_c(4, 2);
    X = ls_solve(A, B);
    cmat X_ref = "1.06773+0.57688i 1.72069-0.274636i;"
                 "0.642428-0.433694i 1.77196+0.66473i;"
                 "0.465216-2.37791i 0.302972-2.9492i;"
                 "0.188992-0.408485i -0.849173-0.00855581i";
    assert_cmat(X_ref, X);

    // Square system (chol)
    A = randn_c(4, 4);
    A = A.transpose() * A;
    b = randn_c(4);
    x = ls_solve(A, b);
    x_ref = "0.656667-0.696584i 0.448907+0.212332i -0.73524-1.06406i 0.991158-0.396795i";
    assert_cvec(x_ref, x);

    // Square system (Chol)
    A = randn_c(4, 4);
    A = A.transpose() * A;
    B = randn_c(4, 2);
    X = ls_solve(A, B);
    X_ref = "0.102859-2.48984i 5.31307-0.364868i;"
            "-0.460786+1.98873i -5.60124+0.500538i;"
            "-1.26409+0.955385i -2.9779-2.32294i;"
            "-1.13999-0.738824i 1.12245-4.405i";
    assert_cmat(X_ref, X);

    // Overdetermined system
    A = randn_c(4, 2);
    b = randn_c(4);
    x = ls_solve_od(A, b);
    x_ref = "0.0737946+0.819157i 0.0650055+0.0621781i";
    assert_cvec(x_ref, x);

    // Overdetermined system
    A = randn_c(4, 2);
    B = randn_c(4, 3);
    X = ls_solve_od(A, B);
    X_ref = "-0.796269+0.064283i -0.57216-0.000698553i 0.365305-0.067377i;"
            "0.845856-0.416823i -0.264445-0.73061i -0.853229+0.231852i";
    assert_cmat(X_ref, X);

    // Underdetermined system
    A = randn_c(2, 4);
    b = randn_c(2);
    x = ls_solve_ud(A, b);
    x_ref = "-0.0430224+0.0471777i -0.0200806-0.0467243i 0.234193-0.2424i 0.421536+0.260783i";
    assert_cvec(x_ref, x);

    // Underdetermined system
    A = randn_c(2, 4);
    B = randn_c(2, 3);
    X = ls_solve_ud(A, B);
    X_ref = "-0.314656-0.0734181i 0.110243-0.333847i 0.084258+0.0813639i;"
            "-0.234372-0.0964481i -0.0627429-0.292354i 0.571474+0.0219301i;"
            "0.0852421-0.00435951i 0.0685549+0.00636267i -0.129461+0.265599i;"
            "-0.23233-0.188923i -0.0575991+0.0723313i 0.472732-0.915352i";
    assert_cmat(X_ref, X);
  }
}
