/*!
 * \file
 * \brief FastICA test program
 * \author Francois Cayre, Teddy Furon and Adam Piatyszek
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

#include <itpp/itsignal.h>
#include <cstdio>

using namespace itpp;
using namespace std;


#if defined(FASTICA_TEST_FILE)

int main()
{
  FILE * fpin = NULL;
  float tmp = 0.0;

  // Separate nrIC independent components in nrSamples samples
  int nrSamples = 0, nrIC = 0;

  fpin = fopen(FASTICA_TEST_FILE, "r");
  if (fpin == 0) {
    cerr << "Error: Could not open FASTICA_TEST_FILE for reading" << endl;
    return 1;
  }

  cout << "=====================================" << endl;
  cout << "   Test program for FastICA / IT++   " << endl;
  cout << "=====================================" << endl;

  int ret = fscanf(fpin, "%d", &nrSamples);
  ret = fscanf(fpin, "%d", &nrIC);

  mat X = zeros(nrIC, nrSamples);

  for (int i = 0; i < nrSamples; i++)
    for (int j = 0; j < nrIC; j++) {
      ret = fscanf(fpin , "%f", &tmp);
      X(j, i) = tmp;
    }

  fclose(fpin);

  // Instantiate an ICA object with default parameters : SYMM approach and
  // POW3 non-linearity
  // Be sure that :
  // - nrSamples = number of samples = nb of columns of the input matrix
  // - nrIC = number of sensors = nb of rows of the input matrix
  cout << "\n==========================================================" << endl;
  cout << "Use SYMM approach and POW3 non-linearity :" << endl;
  Fast_ICA my_fastica(X);

  // Set number of independent components to separate :
  // By default, this value is taken from the dimension of
  // the input data. This line is for illustration purposes.
  // May help in some cases.
  my_fastica.set_nrof_independent_components(nrIC);

  // Perform ICA
  bool result = my_fastica.separate();

  if (result)
  {
    // Get results
    cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
    cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
    cout << "Separated independent components = "
         << my_fastica.get_independent_components() << endl;
  } else
  {
	cout << "Algorithm failed" << endl;
  }

  // Another test with other parameters
  cout << "\n==========================================================" << endl;
  cout << "Use Gaussian non-linearity and deflation approach :" << endl;

  Fast_ICA my_fastica2(X);

  // Set GAUSS non-linearity
  my_fastica2.set_non_linearity(FICA_NONLIN_GAUSS);

  // Use deflation approach : IC are computed one by one
  my_fastica2.set_approach(FICA_APPROACH_DEFL);

  // Perform ICA
  result = my_fastica2.separate();

  if (result)
  {
    // Get results
    cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
    cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
    cout << "Separated independent components = "
         << my_fastica.get_independent_components() << endl;
  } else
  {
	cout << "Algorithm failed" << endl;
  }

  // Another test which should fail
  cout << "\n==========================================================" << endl;
  cout << "Use Gaussian non-linearity and deflation approach :" << endl;

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

  if (result)
  {
    // Get results
    cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
    cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
    cout << "Separated independent components = "
         << my_fastica.get_independent_components() << endl;
  } else
  {
	cout << "Algorithm failed" << endl;
  }

  cout << "\nEnd of Fast_ICA test execution. " << endl;

  return 0;
}

#else

int main()
{
  cerr << "FASTICA_TEST_FILE not defined. Test skipped." << endl;
  return 1;
}

#endif // defined(FASTICA_TEST_FILE)
