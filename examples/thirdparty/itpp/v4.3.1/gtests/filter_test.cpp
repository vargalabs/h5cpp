/*!
 * \file
 * \brief Filter classes test program
 * \author Hakan Eriksson, Thomas Eriksson, Tony Ottosson and Adam Piatyszek
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
#include <itpp/itsignal.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST (FilterTest, All)
{
  static const double eps = 1e-4;
  int i;

  RNG_reset(0);

  // Test signals
  vec x = randn(20), x2 = randn(20);
  cvec cx = randn_c(20), cx2 = randn_c(20);

  // Filter coefficients
  vec b(10);
  b.ones();
  b(0) += 0.1;
  cvec cb(10);
  cb.ones();
  cb(0) += 0.1;

  cvec ca(2);
  ca(0) = 1.0;
  ca(1) = -1.0;
  vec a(2);
  a(0) = 1.0;
  a(1) = -1.0;

  vec y, y2;
  vec s1, s2;
  vec y_ref, y2_ref;
  vec s1_ref, s2_ref;

  cvec cy, cy2;
  cvec cs1, cs2;
  cvec cy_ref, cy2_ref;
  cvec cs1_ref, cs2_ref;

  // MA Filter

  MA_Filter<double, double, double> H(b);
  MA_Filter<complex<double>, double, complex<double> > CH(b);
  MA_Filter<complex<double>, complex<double>, complex<double> > C(cb);

  y = H(x);
  y_ref = "-0.592547 -1.94108 -0.682537 -0.925824 -1.31505 -2.24478 -0.737816 -1.60975 -2.31362 -1.35149 -1.47907 "
  "-4.25692 -5.03592 -5.93111 -5.19747 -5.74787 -6.48426 -6.27389 -5.5178 -5.3668";
  s1 = H.get_state();
  s1_ref = "-0.531968 0.875791 0.0146015 -0.387522 0.387353 -1.27628 0.249247 -0.940147 -0.112755 -3.7327";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }
  for (i = 0; i < s1.length(); ++i) {
    ASSERT_NEAR(s1_ref[i], s1[i], eps);
  }
  y2 = H(x2);
  y2_ref = "-4.50163 -1.47862 -0.920061 1.41155 0.817598 1.60014 1.96435 3.88003 4.36433 3.70982 4.04877 4.65587 "
  "4.90407 4.96362 3.4177 3.86524 2.53372 2.13073 0.89627 -1.23503";
  s2 = H.get_state();
  s2_ref = "0.679012 -1.75261 -0.509272 0.891211 -0.641614 -0.150657 -1.46646 1.28292 0.54896 0.0587334";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }
  H.set_state(s1);
  y2 = H(x2);
  y2_ref = "-4.50163 -1.47862 -0.920061 1.41155 0.817598 1.60014 1.96435 3.88003 4.36433 3.70982 4.04877 4.65587 "
  "4.90407 4.96362 3.4177 3.86524 2.53372 2.13073 0.89627 -1.23503";
  s2 = H.get_state();
  s2_ref = "0.679012 -1.75261 -0.509272 0.891211 -0.641614 -0.150657 -1.46646 1.28292 0.54896 0.0587334";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }

  cy = CH(cx);
  cy_ref = "-1.11701-0.115333i -0.884597+1.64759i -0.00639818+1.80179i 0.914903+3.31264i 1.71424+2.51272i "
  "0.676698+3.40542i 1.56528+3.18382i 1.84761+3.91019i 2.40857+1.97671i 3.26292+1.80654i 3.50956+1.99636i "
  "2.23868-0.661574i 3.27804-0.355509i 0.733909-1.68614i 0.491808-0.360612i 1.35977-1.79598i 0.390552-2.44626i "
  "0.579894-3.85374i 0.0423402-1.57992i -0.921652-1.93196i";
  cs1 = CH.get_state();
  cs1_ref = "-0.623863+0.0491432i -0.121881-0.559813i 0.0421585+0.453723i 0.446095-0.75948i -0.21549-0.76228i "
  "0.034294-0.555875i 0.393913+0.668983i -1.34092+0.103413i 1.58013+0.44975i -1.10391-0.963544i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }
  for (i = 0; i < cs1.length(); ++i) {
    ASSERT_NEAR(cs1_ref[i].real(), cs1[i].real(), eps);
    ASSERT_NEAR(cs1_ref[i].imag(), cs1[i].imag(), eps);
  }
  cy2 = CH(cx2);
  cs2 = CH.get_state();
  cy2_ref = "-0.102301-1.48432i 1.13288+0.526827i -0.277328-1.46198i -0.44344-0.488853i -1.85361-1.42378i "
  "-2.43148-0.667385i -2.82985+0.0317305i -3.51517-0.266065i -3.17257-0.201244i -2.63849-0.197425i "
  "-1.57769-0.050627i -1.24298-1.703i -1.07173-1.34898i 1.16861-2.30151i 2.36235-4.17338i 2.2602-2.5629i "
  "2.23618-2.54142i 3.479-0.439856i 1.15184-0.162033i 1.71127+1.05949i";
  cs2_ref = "1.1526+0.454975i 0.717425+0.745976i -1.74527+0.696227i 0.824019+1.04348i -0.633288+0.112685i "
  "-0.609585+1.44666i 0.209481-1.86573i 0.837287-0.168026i 0.355532-0.91947i 0.531331-0.56189i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }

  CH.set_state(cs1);
  cy2 = CH(cx2);
  cs2 = CH.get_state();
  cy2_ref = "-0.102301-1.48432i 1.13288+0.526827i -0.277328-1.46198i -0.44344-0.488853i -1.85361-1.42378i "
  "-2.43148-0.667385i -2.82985+0.0317305i -3.51517-0.266065i -3.17257-0.201244i -2.63849-0.197425i "
  "-1.57769-0.050627i -1.24298-1.703i -1.07173-1.34898i 1.16861-2.30151i 2.36235-4.17338i 2.2602-2.5629i "
  "2.23618-2.54142i 3.479-0.439856i 1.15184-0.162033i 1.71127+1.05949i";
  cs2_ref = "1.1526+0.454975i 0.717425+0.745976i -1.74527+0.696227i 0.824019+1.04348i -0.633288+0.112685i "
  "-0.609585+1.44666i 0.209481-1.86573i 0.837287-0.168026i 0.355532-0.91947i 0.531331-0.56189i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }

  cy = C(cx);
  cy_ref = "-1.11701-0.115333i -0.884597+1.64759i -0.00639818+1.80179i 0.914903+3.31264i 1.71424+2.51272i "
  "0.676698+3.40542i 1.56528+3.18382i 1.84761+3.91019i 2.40857+1.97671i 3.26292+1.80654i 3.50956+1.99636i "
  "2.23868-0.661574i 3.27804-0.355509i 0.733909-1.68614i 0.491808-0.360612i 1.35977-1.79598i 0.390552-2.44626i "
  "0.579894-3.85374i 0.0423402-1.57992i -0.921652-1.93196i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  y =  filter(b, 1, x);
  y_ref = "-0.592547 -1.94108 -0.682537 -0.925824 -1.31505 -2.24478 -0.737816 -1.60975 -2.31362 -1.35149 "
  "-1.47907 -4.25692 -5.03592 -5.93111 -5.19747 -5.74787 -6.48426 -6.27389 -5.5178 -5.3668";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }

  cy =  filter(b, 1, cx);
  cy_ref = "-1.11701-0.115333i -0.884597+1.64759i -0.00639818+1.80179i 0.914903+3.31264i 1.71424+2.51272i "
  "0.676698+3.40542i 1.56528+3.18382i 1.84761+3.91019i 2.40857+1.97671i 3.26292+1.80654i 3.50956+1.99636i "
  "2.23868-0.661574i 3.27804-0.355509i 0.733909-1.68614i 0.491808-0.360612i 1.35977-1.79598i 0.390552-2.44626i "
  "0.579894-3.85374i 0.0423402-1.57992i -0.921652-1.93196i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  cy =  filter(cb, 1, cx);
  cy_ref = "-1.11701-0.115333i -0.884597+1.64759i -0.00639818+1.80179i 0.914903+3.31264i 1.71424+2.51272i "
  "0.676698+3.40542i 1.56528+3.18382i 1.84761+3.91019i 2.40857+1.97671i 3.26292+1.80654i 3.50956+1.99636i "
  "2.23868-0.661574i 3.27804-0.355509i 0.733909-1.68614i 0.491808-0.360612i 1.35977-1.79598i 0.390552-2.44626i "
  "0.579894-3.85374i 0.0423402-1.57992i -0.921652-1.93196i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  // AR Filter
  AR_Filter<double, double, double> HAR(a);
  AR_Filter<complex<double>, double, complex<double> > CHAR(a);
  AR_Filter<complex<double>, complex<double>, complex<double> > CAR(ca);

  y = HAR(x);
  y_ref = "-0.538679 -1.81359 -0.78536 -0.913054 -1.2785 -2.15694 -0.866827 -1.54221 -2.24349 "
  "-1.43258 -1.96455 -5.69725 -5.81 -6.75015 -6.5009 -7.77718 -7.38983 -7.77735 -7.76275 -6.88696";
  s1 = HAR.get_state();
  s1_ref = "-6.88696";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }
  for (i = 0; i < s1.length(); ++i) {
    ASSERT_NEAR(s1_ref[i], s1[i], eps);
  }
  y2 = HAR(x2);
  y2_ref = "-6.50443 -7.11482 -6.76504 -5.46828 -5.66376 -6.13037 -5.48956 -4.04207 -3.45694 "
  "-3.20257 -2.52356 -2.46483 -1.91587 -0.632949 -2.09941 -2.25006 -2.89168 -2.00047 -2.50974 -4.26235";
  s2 = HAR.get_state();
  s2_ref = "-4.26235";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }
  HAR.set_state(s1);
  y2 = HAR(x2);
  y2_ref = "-6.50443 -7.11482 -6.76504 -5.46828 -5.66376 -6.13037 -5.48956 -4.04207 -3.45694 -3.20257 "
  "-2.52356 -2.46483 -1.91587 -0.632949 -2.09941 -2.25006 -2.89168 -2.00047 -2.50974 -4.26235";
  s2 = HAR.get_state();
  s2_ref = "-4.26235";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }

  cy = CHAR(cx);
  cs1 = CHAR.get_state();
  cy_ref = "-1.01547-0.104848i -0.896494+1.48828i -0.087316+1.77329i 0.823792+3.1727i 1.63329+2.57271i "
  "0.763661+3.32972i 1.49241+3.19708i 1.81532+3.84536i 2.35464+2.14659i 3.18034+1.83746i 2.55648+1.8866i "
  "1.45258+0.923057i 3.03271+1.37281i 1.69179+1.47622i 2.08571+2.1452i 2.12+1.58933i 1.90451+0.827048i "
  "2.3506+0.0675678i 2.39276+0.521291i 2.27088-0.0385215i";
  cs1_ref = "2.27088-0.0385215i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }
  for (i = 0; i < cs1.length(); ++i) {
    ASSERT_NEAR(cs1_ref[i].real(), cs1[i].real(), eps);
    ASSERT_NEAR(cs1_ref[i].imag(), cs1[i].imag(), eps);
  }
  cy2 = CHAR(cx2);
  cs2 = CHAR.get_state();
  cy2_ref = "2.43752+0.36221i 2.57201+1.351i 2.73871+0.041755i 1.38384+0.901402i 0.336797+0.737782i "
  "-0.252542+0.9052i -0.864179+0.862998i -1.13726-0.101997i -0.8123+0.28168i -0.408032-0.18889i "
  "0.744568+0.266086i 1.2759-0.295804i 1.63143-1.21527i 2.46872-1.3833i 2.6782-3.24903i 2.06861-1.80237i "
  "1.43533-1.68969i 2.25935-0.646202i 0.514072+0.0500244i 1.2315+0.796001i";
  cs2_ref = "1.2315+0.796001i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }
  CHAR.set_state(cs1);
  cy2 = CHAR(cx2);
  cs2 = CHAR.get_state();
  cy2_ref = "2.43752+0.36221i 2.57201+1.351i 2.73871+0.041755i 1.38384+0.901402i 0.336797+0.737782i "
  "-0.252542+0.9052i -0.864179+0.862998i -1.13726-0.101997i -0.8123+0.28168i -0.408032-0.18889i "
  "0.744568+0.266086i 1.2759-0.295804i 1.63143-1.21527i 2.46872-1.3833i 2.6782-3.24903i 2.06861-1.80237i "
  "1.43533-1.68969i 2.25935-0.646202i 0.514072+0.0500244i 1.2315+0.796001i";
  cs2_ref = "1.2315+0.796001i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }

  cy = CAR(cx);
  cy_ref = "-1.01547-0.104848i -0.896494+1.48828i -0.087316+1.77329i 0.823792+3.1727i "
  "1.63329+2.57271i 0.763661+3.32972i 1.49241+3.19708i 1.81532+3.84536i 2.35464+2.14659i "
  "3.18034+1.83746i 2.55648+1.8866i 1.45258+0.923057i 3.03271+1.37281i 1.69179+1.47622i "
  "2.08571+2.1452i 2.12+1.58933i 1.90451+0.827048i 2.3506+0.0675678i 2.39276+0.521291i 2.27088-0.0385215i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  y =  filter(1, a, x);
  y_ref = "-0.538679 -1.81359 -0.78536 -0.913054 -1.2785 -2.15694 -0.866827 -1.54221 -2.24349 -1.43258 "
  "-1.96455 -5.69725 -5.81 -6.75015 -6.5009 -7.77718 -7.38983 -7.77735 -7.76275 -6.88696";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }

  cy =  filter(1, a, cx);
  cy_ref = "-1.01547-0.104848i -0.896494+1.48828i -0.087316+1.77329i 0.823792+3.1727i 1.63329+2.57271i "
  "0.763661+3.32972i 1.49241+3.19708i 1.81532+3.84536i 2.35464+2.14659i 3.18034+1.83746i 2.55648+1.8866i "
  "1.45258+0.923057i 3.03271+1.37281i 1.69179+1.47622i 2.08571+2.1452i 2.12+1.58933i 1.90451+0.827048i "
  "2.3506+0.0675678i 2.39276+0.521291i 2.27088-0.0385215i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  cy =  filter(1, ca, cx);
  cy_ref = "-1.01547-0.104848i -0.896494+1.48828i -0.087316+1.77329i 0.823792+3.1727i 1.63329+2.57271i "
  "0.763661+3.32972i 1.49241+3.19708i 1.81532+3.84536i 2.35464+2.14659i 3.18034+1.83746i 2.55648+1.8866i "
  "1.45258+0.923057i 3.03271+1.37281i 1.69179+1.47622i 2.08571+2.1452i 2.12+1.58933i 1.90451+0.827048i "
  "2.3506+0.0675678i 2.39276+0.521291i 2.27088-0.0385215i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  // ARMA Filter

  ARMA_Filter<double, double, double> HARMA(b, a);
  ARMA_Filter<complex<double>, double, complex<double> > CHARMA(b, a);
  ARMA_Filter<complex<double>, complex<double>, complex<double> > CARMA(cb, ca);

  y = HARMA(x);
  s1 = HARMA.get_state();
  y_ref = "-0.592547 -2.53363 -3.21617 -4.14199 -5.45704 -7.70182 -8.43964 -10.0494 -12.363 "
  "-13.7145 -15.1936 -19.4505 -24.4864 -30.4175 -35.615 -41.3629 -47.8471 -54.121 -59.6388 -65.0056";
  s1_ref = "-6.88696 -7.76275 -7.77735 -7.38983 -7.77718 -6.5009 -6.75015 -5.81 -5.69725";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }
  for (i = 0; i < s1.length(); ++i) {
    ASSERT_NEAR(s1_ref[i], s1[i], eps);
  }
  y2 = HARMA(x2);
  s2 = HARMA.get_state();
  y2_ref = "-69.5072 -70.9858 -71.9059 -70.4944 -69.6768 -68.0766 -66.1123 -62.2322 -57.8679 -54.1581 "
  "-50.1093 -45.4535 -40.5494 -35.5858 -32.1681 -28.3028 -25.7691 -23.6384 -22.7421 -23.9771";
  s2_ref = "-4.26235 -2.50974 -2.00047 -2.89168 -2.25006 -2.09941 -0.632949 -1.91587 -2.46483";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }
  HARMA.set_state(s1);
  y2 = HARMA(x2);
  s2 = HARMA.get_state();
  y2_ref = "-69.5072 -70.9858 -71.9059 -70.4944 -69.6768 -68.0766 -66.1123 -62.2322 -57.8679 -54.1581 "
  "-50.1093 -45.4535 -40.5494 -35.5858 -32.1681 -28.3028 -25.7691 -23.6384 -22.7421 -23.9771";
  s2_ref = "-4.26235 -2.50974 -2.00047 -2.89168 -2.25006 -2.09941 -0.632949 -1.91587 -2.46483";
  for (i = 0; i < y2.length(); ++i) {
    ASSERT_NEAR(y2_ref[i], y2[i], eps);
  }
  for (i = 0; i < s2.length(); ++i) {
    ASSERT_NEAR(s2_ref[i], s2[i], eps);
  }

  cy = CHARMA(cx);
  cs1 = CHARMA.get_state();
  cy_ref = "-1.11701-0.115333i -2.00161+1.53226i -2.00801+3.33405i -1.09311+6.64669i 0.621131+9.1594i "
  "1.29783+12.5648i 2.86311+15.7486i 4.71072+19.6588i 7.11929+21.6355i 10.3822+23.4421i 13.8918+25.4384i "
  "16.1304+24.7769i 19.4085+24.4214i 20.1424+22.7352i 20.6342+22.3746i 21.994+20.5786i 22.3845+18.1324i "
  "22.9644+14.2786i 23.0068+12.6987i 22.0851+10.7667i";
  cs1_ref = "2.27088-0.0385215i 2.39276+0.521291i 2.3506+0.0675678i 1.90451+0.827048i 2.12+1.58933i "
  "2.08571+2.1452i 1.69179+1.47622i 3.03271+1.37281i 1.45258+0.923057i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }
  for (i = 0; i < cs1.length(); ++i) {
    ASSERT_NEAR(cs1_ref[i].real(), cs1[i].real(), eps);
    ASSERT_NEAR(cs1_ref[i].imag(), cs1[i].imag(), eps);
  }

  cy2 = CHARMA(cx2);
  cs2 = CHARMA.get_state();
  cy2_ref = "21.9828+9.28243i 23.1157+9.80926i 22.8384+8.34728i 22.3949+7.85843i 20.5413+6.43465i "
  "18.1098+5.76726i 15.28+5.79899i 11.7648+5.53293i 8.59225+5.33168i 5.95376+5.13426i 4.37607+5.08363i "
  "3.1331+3.38063i 2.06137+2.03165i 3.22998-0.269851i 5.59233-4.44324i 7.85252-7.00614i 10.0887-9.54756i "
  "13.5677-9.98741i 14.7195-10.1494i 16.4308-9.08996i";
  cs2_ref = "1.2315+0.796001i 0.514072+0.0500244i 2.25935-0.646202i 1.43533-1.68969i 2.06861-1.80237i "
  "2.6782-3.24903i 2.46872-1.3833i 1.63143-1.21527i 1.2759-0.295804i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }

  CHARMA.set_state(cs1);
  cy2 = CHARMA(cx2);
  cs2 = CHARMA.get_state();
  cy2_ref = "21.9828+9.28243i 23.1157+9.80926i 22.8384+8.34728i 22.3949+7.85843i 20.5413+6.43465i "
  "18.1098+5.76726i 15.28+5.79899i 11.7648+5.53293i 8.59225+5.33168i 5.95376+5.13426i 4.37607+5.08363i "
  "3.1331+3.38063i 2.06137+2.03165i 3.22998-0.269851i 5.59233-4.44324i 7.85252-7.00614i 10.0887-9.54756i "
  "13.5677-9.98741i 14.7195-10.1494i 16.4308-9.08996i";
  cs2_ref = "1.2315+0.796001i 0.514072+0.0500244i 2.25935-0.646202i 1.43533-1.68969i 2.06861-1.80237i "
  "2.6782-3.24903i 2.46872-1.3833i 1.63143-1.21527i 1.2759-0.295804i";
  for (i = 0; i < cy2.length(); ++i) {
    ASSERT_NEAR(cy2_ref[i].real(), cy2[i].real(), eps);
    ASSERT_NEAR(cy2_ref[i].imag(), cy2[i].imag(), eps);
  }
  for (i = 0; i < cs2.length(); ++i) {
    ASSERT_NEAR(cs2_ref[i].real(), cs2[i].real(), eps);
    ASSERT_NEAR(cs2_ref[i].imag(), cs2[i].imag(), eps);
  }

  cy = CARMA(cx);
  cy_ref = "-1.11701-0.115333i -2.00161+1.53226i -2.00801+3.33405i -1.09311+6.64669i 0.621131+9.1594i "
  "1.29783+12.5648i 2.86311+15.7486i 4.71072+19.6588i 7.11929+21.6355i 10.3822+23.4421i 13.8918+25.4384i "
  "16.1304+24.7769i 19.4085+24.4214i 20.1424+22.7352i 20.6342+22.3746i 21.994+20.5786i 22.3845+18.1324i "
  "22.9644+14.2786i 23.0068+12.6987i 22.0851+10.7667i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  y =  filter(b, a, x);
  y_ref = "-0.592547 -2.53363 -3.21617 -4.14199 -5.45704 -7.70182 -8.43964 -10.0494 -12.363 -13.7145 "
  "-15.1936 -19.4505 -24.4864 -30.4175 -35.615 -41.3629 -47.8471 -54.121 -59.6388 -65.0056";
  for (i = 0; i < y.length(); ++i) {
    ASSERT_NEAR(y_ref[i], y[i], eps);
  }

  cy =  filter(b, a, cx);
  cy_ref = "-1.11701-0.115333i -2.00161+1.53226i -2.00801+3.33405i -1.09311+6.64669i 0.621131+9.1594i "
  "1.29783+12.5648i 2.86311+15.7486i 4.71072+19.6588i 7.11929+21.6355i 10.3822+23.4421i 13.8918+25.4384i "
  "16.1304+24.7769i 19.4085+24.4214i 20.1424+22.7352i 20.6342+22.3746i 21.994+20.5786i 22.3845+18.1324i "
  "22.9644+14.2786i 23.0068+12.6987i 22.0851+10.7667i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }

  cy =  filter(cb, ca, cx);
  cy_ref = "-1.11701-0.115333i -2.00161+1.53226i -2.00801+3.33405i -1.09311+6.64669i 0.621131+9.1594i "
  "1.29783+12.5648i 2.86311+15.7486i 4.71072+19.6588i 7.11929+21.6355i 10.3822+23.4421i 13.8918+25.4384i "
  "16.1304+24.7769i 19.4085+24.4214i 20.1424+22.7352i 20.6342+22.3746i 21.994+20.5786i 22.3845+18.1324i "
  "22.9644+14.2786i 23.0068+12.6987i 22.0851+10.7667i";
  for (i = 0; i < cy.length(); ++i) {
    ASSERT_NEAR(cy_ref[i].real(), cy[i].real(), eps);
    ASSERT_NEAR(cy_ref[i].imag(), cy[i].imag(), eps);
  }
}
