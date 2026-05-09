/*!
 * \file
 * \brief Definition of FastICA (Independent Component Analysis) for IT++
 * \author Francois Cayre and Teddy Furon
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
 *
 * This is IT++ implementation of the original Matlab package FastICA.
 *
 * This code is Copyright (C) 2004 by:
 *   Francois CAYRE and Teddy FURON
 *   TEMICS Project
 *   INRIA/Rennes (IRISA)
 *   Campus Universitaire de Beaulieu
 *   35042 RENNES cedex FRANCE
 *
 * Email: firstname.lastname@irisa.fr
 *
 * Matlab package is Copyright (C) 1998 by:
 *   Jarmo HURRI, Hugo GAVERT, Jaakko SARELA and Aapo HYVARINEN
 *   Laboratory of Information and Computer Science
 *   Helsinki University of Technology
 *
 * URL: http://www.cis.hut.fi/projects/ica/fastica/about.shtml
 *
 * If you use results given by this FastICA software in an article for
 * a scientific journal, conference proceedings or similar, please
 * include the following original reference in the bibliography:
 *
 *   A. Hyvarinen, Fast and Robust Fixed-Point Algorithms for
 *   Independent Component Analysis, IEEE Transactions on Neural
 *   Networks 10(3):626-634, 1999
 *
 * Differences with the original Matlab implementation:
 * - no GUI
 * - return something even in the case of a convergence problem
 * - optimization of SVD decomposition (performed 2 times in Matlab,
 *   only 1 time in IT++)
 * - default approach is SYMM with non-linearity POW3
 */

#ifndef FASTICA_H
#define FASTICA_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

//! Use deflation approach : compute IC one-by-one in a Gram-Schmidt-like fashion
#define FICA_APPROACH_DEFL 2
//! Use symmetric approach : compute all ICs at a time
#define FICA_APPROACH_SYMM 1

//! Use x^3 non-linearity
#define FICA_NONLIN_POW3 10
//! Use tanh(x) non-linearity
#define FICA_NONLIN_TANH 20
//! Use Gaussian non-linearity
#define FICA_NONLIN_GAUSS 30
//! Use skew non-linearity
#define FICA_NONLIN_SKEW 40

//! Set random start for Fast_ICA
#define FICA_INIT_RAND  0
//! Set predefined start for Fast_ICA
#define FICA_INIT_GUESS 1

//! Eigenvalues of the covariance matrix lower than FICA_TOL are discarded for analysis
#define FICA_TOL 1e-9

namespace itpp
{

/*!
  \addtogroup fastica
*/

//---------------------- FastICA --------------------------------------

/*!
\brief Fast_ICA Fast Independent Component Analysis (Fast ICA)
\ingroup fastica

The software is based upon original FastICA for Matlab from
A. Hyvarinen. Fast and Robust Fixed-Point Algorithms for
Independent Component Analysis.  IEEE Transactions on Neural
Networks, 10(3), pp. 626-634, 1999.

Example:
\code
FastICA fastica(sources);
fastica.set_nrof_independent_components(sources.rows());
fastica.set_non_linearity(  FICA_NONLIN_TANH );
fastica.set_approach( FICA_APPROACH_DEFL );
fastica.separate();
mat ICs = fastica.get_independent_components();
\endcode
*/
class ITPP_EXPORT Fast_ICA
{

public:

  /*!
    \brief Constructor

    Construct a Fast_ICA object with mixed signals to separate.

    \param ma_mixed_sig (Input) Mixed signals to separate
  */
  Fast_ICA(mat ma_mixed_sig);

  /*!
    \brief Explicit launch of main FastICA function

    Explicit launch of the Fast_ICA algorithm.
	\returns true if algorithm converged and false otherwise
  */
  bool separate(void);

  /*!
    \brief Set approach : FICA_APPROACH_DEFL or FICA_APPROACH_SYMM (default)

    Set approach to use : FICA_APPROACH_SYMM (symmetric) or FICA_APPROACH_DEFL (deflation). The symmetric approach computes all ICs at a time, whereas the deflation approach computes them one by one.

    \param in_approach (Input) Type of approach to use
  */
  void set_approach(int in_approach);

  /*!
    \brief Set number of independent components to separate

    Set the number of ICs to compute.

    \param in_nrIC (Input) Number of ICs to compute
  */
  void set_nrof_independent_components(int in_nrIC);

  /*!
    \brief Set non-linearity

    Set non-linearity to use : FICA_NONLIN_POW3 (default),  FICA_NONLIN_TANH, FICA_NONLIN_GAUSS, FICA_NONLIN_SKEW

    \param in_g (Input) Non-linearity. Can be selected from FICA_NONLIN_POW3, FICA_NONLIN_TANH, FICA_NONLIN_GAUSS or FICA_NONLIN_SKEW
  */
  void set_non_linearity(int in_g);

  /*!
    \brief Set fine tuning

    Set fine tuning true or false.

    \param in_finetune (Input) Boolean (true or false)
  */
  void set_fine_tune(bool in_finetune);

  /*!
    \brief Set \f$a_1\f$ parameter

    Set internal parameter \f$a_1\f$ of Fast_ICA (See reference paper).

    \param fl_a1 (Input) Parameter \f$a_1\f$ from reference paper
  */
  void set_a1(double fl_a1);

  /*!
    \brief Set \f$a_2\f$ parameter

    Set internal parameter \f$a_2\f$ of Fast_ICA (See reference paper).

    \param fl_a2 (Input) Parameter \f$a_2\f$ from reference paper
  */
  void set_a2(double fl_a2);

  /*!
    \brief Set \f$\mu\f$ parameter

    Set internal parameter \f$\mu\f$ of Fast_ICA (See reference paper).

    \param fl_mu (Input) Parameter \f$\mu\f$ from reference paper
  */
  void set_mu(double fl_mu);

  /*!
    \brief Set convergence parameter \f$\epsilon\f$

    Set \f$\epsilon\f$ parameter for convergence precision.

    \param fl_epsilon (Input) \f$\epsilon\f$ is convergence precision
  */
  void set_epsilon(double fl_epsilon);

  /*!
    \brief Set sample size

    Set the percentage of samples to take into account at every iteration.

    \param fl_sampleSize (Input) Percentage of data to take into account at every iteration
  */
  void set_sample_size(double fl_sampleSize);

  /*!
    \brief Set stabilization mode true or off

    Set stabilization mode.

    \param in_stabilization (Input) Set stabilization true or false
  */
  void set_stabilization(bool in_stabilization);

  /*!
    \brief Set maximum number of iterations

    Set maximum number of iterations for Fast_ICA.

    \param in_maxNumIterations (Input) Maximum number of iterations to go through
  */
  void set_max_num_iterations(int in_maxNumIterations);

  /*!
    \brief Set maximum number of iterations for fine tuning

    Set maximum numberr of iterations for fine tuning.

    \param in_maxFineTune (Input) Maximum number of iterations for fine tuning stage
  */
  void set_max_fine_tune(int in_maxFineTune);

  /*!
    \brief Set first eigenvalue index to take into account

    Set first eigenvalue index to take into account.

    \param in_firstEig (Input) First eigenvalue index to take into account
  */
  void set_first_eig(int in_firstEig);

  /*!
    \brief Set last eigenvalue index to take into account

    Set last eigenvalue index to take into account.

    \param in_lastEig (Input) Last eigenvalue index to take into account
  */
  void set_last_eig(int in_lastEig);

  /*!
    \brief If true, only perform Principal Component Analysis (default = false)

    Wether to perform PCA only or PCA+ICA.

    \param in_PCAonly (Input) True = PCA only, false = PCA+ICA (default)
  */
  void set_pca_only(bool in_PCAonly);

  /*!
    \brief Set initial guess matrix instead of random (default)

    Set initial matrix instead of random matrix.

    \param ma_initGuess (Input) Initial guess matrix
  */
  void set_init_guess(mat ma_initGuess);


  /*!
    \brief Get mixing matrix

    Return mixing matrix.

    \return Mixing matrix
  */
  mat get_mixing_matrix();

  /*!
    \brief Get separating matrix

    Return separating matrix.

    \return Separating matrix
  */
  mat get_separating_matrix();

  /*!
    \brief Get separated signals

    Return separated signals (Independent Components).

    \return ICs
  */
  mat get_independent_components();

  /*!
    \brief Get number of independent components

    Return number of ICs.

    \return Number of ICs
  */
  int get_nrof_independent_components();

  /*!
    \brief Get nrIC first columns of the de-whitening matrix

    Return principal eigenvectors.

    \return Principal eigenvectors
  */
  mat get_principal_eigenvectors();

  /*!
    \brief Get the whitening matrix

    Return whitening matrix.

    \return Whitening matrix
  */
  mat get_whitening_matrix();

  /*!
    \brief Get the de-whitening matrix

    Return dewhitening matrix.

    \return Dewhitening matrix
  */
  mat get_dewhitening_matrix();

  /*!
    \brief Get whitened signals

    Return whitened signals.

    \return Whitened signals
  */
  mat get_white_sig();

private:

  int approach, numOfIC, g, initState;
  bool finetune, stabilization, PCAonly;
  double a1, a2, mu, epsilon, sampleSize;
  int maxNumIterations, maxFineTune;

  int firstEig, lastEig;

  mat initGuess;

  mat mixedSig, A, W, icasig;

  mat whiteningMatrix;
  mat dewhiteningMatrix;
  mat whitesig;

  mat E, VecPr;
  vec D;

}; // class Fast_ICA

} // namespace itpp


#endif // #ifndef FASTICA_H
