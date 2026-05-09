/*!
 * \file
 * \brief FastICA test program
 * \author Francois Cayre, Teddy Furon and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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

/*----------------------------------------------------------------------------------------*
 * FastICA for IT++                                                                       *
 *----------------------------------------------------------------------------------------*
 * This code is Copyright (C) 2004 by Francois CAYRE and Teddy FURON                      *
 *                                    TEMICS Project                                      *
 *                                    INRIA/Rennes (IRISA)                                *
 *                                    Campus Universitaire de Beaulieu                    *
 *                                    35042 RENNES cedex FRANCE                           *
 *                                                                                        *
 * Email : firstname.lastname@irisa.fr                                                    *
 *                                                                                        *
 * This is the IT++ implementation of the original Matlab package FastICA.                *
 *                                                                                        *
 * Matlab package is Copyright (C) 1998 by Jarmo HURRI, Hugo GAVERT, Jaakko SARELA and    *
 *                                         Aapo HYVARINEN                                 *
 *                                         Laboratory of Information and Computer Science *
 *                                         Helsinki University of Technology              *
 *                                                                                        *
 * URL : http://www.cis.hut.fi/projects/ica/fastica/about.shtml                           *
 *                                                                                        *
 * If you use results given by this FastICA software in an article for a scientific       *
 * journal, conference proceedings or similar, please include the following original      *
 * reference in the bibliography :                                                        *
 *                                                                                        *
 *     A. Hyvarinen. Fast and Robust Fixed-Point Algorithms for Independent Component     *
 *     Analysis. IEEE Transactions on Neural Networks 10(3):626-634, 1999.                *
 *----------------------------------------------------------------------------------------*
 *                                           DISCLAIMER                                   *
 *                                                                                        *
 * This software package is free software ; you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by the Free Software    *
 * Foundation ; either version 2 of the License, or any later version.                    *
 *                                                                                        *
 * The software package is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR     *
 * A PARTICULAR PURPOSE.                                                                  *
 * See the GNU General Public License for more details.                                   *
 *                                                                                        *
 *----------------------------------------------------------------------------------------*
 * Differences with the original Matlab implementation :                                  *
 * - no GUI                                                                               *
 * - return something even in the case of a convergence problem                           *
 * - optimization of SVD decomposition (performed 2 times in Matlab, only 1 time in IT++) *
 * - default approach is SYMM wit non-linearity POW3                                      *
 *----------------------------------------------------------------------------------------*/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <itpp/itsignal.h>
#include <cstdio>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST(Fastica, All)
{
#if defined(FASTICA_TEST_FILE)

  const double eps = 1e-12;
  const mat mixing_matrix_expect = "0.20144255910862657 0.37985941512204879; -0.38969201072547255 0.31918961738834795";
  const mat separating_matrix_expect = "1.503295812136356 -1.789033968666887; 1.835342805132183 0.94873936681213422";
  const mat independent_components_expect = "0.081523541458408566 0.74717136112026661 1.2323127422446056 1.4178953014486217 "
                                            "1.2582094210524861 0.79207900492589467 0.13335632609814074 -0.55695312005441711 "
                                            "-1.1101100814898155 -1.3909552319850542;"
                                            "-1.7138482673446558 -1.5321013719589991 -1.352666385717076 -1.1770680527274013 "
                                            "-1.0058918059024045 -0.83864043309954472 -0.67385567733571405 -0.50947550907352301 "
                                            "-0.34333866810304298 -0.17371418886128326";//sampled independent comp matrix

  FILE * fpin = NULL;
  float tmp = 0.0;

  // Separate nrIC independent components in nrSamples samples
  int nrSamples = 0, nrIC = 0;

  fpin = fopen(FASTICA_TEST_FILE, "r");
  ASSERT_TRUE(NULL != fpin) << "Error: Could not open FASTICA_TEST_FILE for reading";

  //Test program for FastICA / IT++

  int ret = fscanf(fpin, "%d", &nrSamples);
  ret = fscanf(fpin, "%d", &nrIC);

  mat X = zeros(nrIC, nrSamples);

  for (int i = 0; i < nrSamples; i++)
    for (int j = 0; j < nrIC; j++) {
      ret = fscanf(fpin , "%f", &tmp);
      X(j, i) = tmp;
    }

  fclose(fpin);

  RNG_reset(0);//separate uses randu

  // Instantiate an ICA object with default parameters : SYMM approach and
  // POW3 non-linearity
  // Be sure that :
  // - nrSamples = number of samples = nb of columns of the input matrix
  // - nrIC = number of sensors = nb of rows of the input matrix
  Fast_ICA my_fastica(X);

  // Set number of independent components to separate :
  // By default, this value is taken from the dimension of
  // the input data. This line is for illustration purposes.
  // May help in some cases.
  my_fastica.set_nrof_independent_components(nrIC);

  // Perform ICA
  bool result = my_fastica.separate();
  ASSERT_TRUE(result) << "Algorithm failed";

  // Get results
  mat mixing_matrix = my_fastica.get_mixing_matrix();
  mat separating_matrix = my_fastica.get_separating_matrix();
  mat independent_components = my_fastica.get_independent_components();

  for (int i = 0; i < mixing_matrix.rows(); ++i) {
    for (int j = 0; j < mixing_matrix.cols(); ++j) {
      ASSERT_NEAR(mixing_matrix_expect(i,j), mixing_matrix(i,j), eps);
      ASSERT_NEAR(separating_matrix_expect(i,j), separating_matrix(i,j), eps);
    }
  }
  ASSERT_EQ(nrIC, independent_components.rows());
  ASSERT_EQ(nrSamples, independent_components.cols());
  for (int i = 0; i < nrIC; ++i) {
    for (int j = 0; j < nrSamples/100; ++j) {
      ASSERT_NEAR(independent_components_expect(i,j), independent_components(i,j*nrSamples/100), eps);
    }
  }

  // Another test with other parameters
  Fast_ICA my_fastica2(X);

  // Set GAUSS non-linearity
  my_fastica2.set_non_linearity(FICA_NONLIN_GAUSS);

  // Use deflation approach : IC are computed one by one
  my_fastica2.set_approach(FICA_APPROACH_DEFL);

  // Perform ICA
  result = my_fastica2.separate();
  ASSERT_TRUE(result) << "Algorithm failed";

  // Get results
  mixing_matrix = my_fastica.get_mixing_matrix();
  separating_matrix = my_fastica.get_separating_matrix();
  independent_components = my_fastica.get_independent_components();

  for (int i = 0; i < mixing_matrix.rows(); ++i) {
    for (int j = 0; j < mixing_matrix.cols(); ++j) {
      ASSERT_NEAR(mixing_matrix_expect(i,j), mixing_matrix(i,j), eps);
      ASSERT_NEAR(separating_matrix_expect(i,j), separating_matrix(i,j), eps);
    }
  }
  ASSERT_EQ(nrIC, independent_components.rows());
  ASSERT_EQ(nrSamples, independent_components.cols());
  for (int i = 0; i < nrIC; ++i) {
    for (int j = 0; j < nrSamples/100; ++j) {
      ASSERT_NEAR(independent_components_expect(i,j), independent_components(i,j*nrSamples/100), eps);
    }
  }

  // Another test which should fail
  const int rows = 10;
  const int comp = 3;
  RNG_reset(1);
  mat signal = randu(rows, 100);
  mat guess = zeros(rows, comp);

  Fast_ICA my_fastica3(signal);

  // Use deflation approach : IC are computed one by one
  my_fastica3.set_approach(FICA_APPROACH_DEFL);
  my_fastica3.set_nrof_independent_components(comp);
  my_fastica3.set_init_guess(guess);
  my_fastica3.set_max_num_iterations(100);

  // Perform ICA
  result = my_fastica3.separate();
  ASSERT_FALSE(result) << "Algorithm should fail";

#else

  FAIL() << "FASTICA_TEST_FILE not defined. Test skipped.";

#endif // defined(FASTICA_TEST_FILE)
}

