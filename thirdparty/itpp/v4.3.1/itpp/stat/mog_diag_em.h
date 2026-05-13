/*!
 * \file
 * \brief Expectation Maximisation (EM) based optimisers for MOG - header file
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

#ifndef MOG_DIAG_EM_H
#define MOG_DIAG_EM_H

#include <itpp/stat/mog_diag.h>
#include <itpp/itexports.h>
#include <itpp/base/base_exports.h>

namespace itpp
{

/*!
  \brief support class for MOG_diag_ML() and MOG_diag_MAP()
  \author Conrad Sanderson
*/
class ITPP_EXPORT MOG_diag_EM_sup : public MOG_diag
{

public:

  //! Default constructor
  MOG_diag_EM_sup() { verbose = false; }

  //! Default destructor
  ~MOG_diag_EM_sup() { }

  //! ADD DOCUMENTATION HERE
  void ml(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in = 10, double var_floor_in = 0.0, double weight_floor_in = 0.0, bool verbose_in = false);
  //! ADD DOCUMENTATION HERE
  void map(MOG_diag &model_in, MOG_diag &prior_model, Array<vec> &X_in, int max_iter_in = 10, double alpha_in = 0.5, double var_floor_in = 0.0, double weight_floor_in = 0.0, bool verbose_in = false);

protected:

  //! Whether we print the progress
  bool verbose;

  //! number of training vectors
  int N;

  //! Maximum number of iterations
  int max_iter;

  //! 'C' pointers to training vectors
  double ** c_X;

  //! ADD DOCUMENTATION HERE
  double var_floor;
  //! ADD DOCUMENTATION HERE
  double weight_floor;

  //! ADD DOCUMENTATION HERE
  void inline update_internals();
  //! ADD DOCUMENTATION HERE
  void inline sanitise_params();
  //! ADD DOCUMENTATION HERE
  double ml_update_params();
  //! ADD DOCUMENTATION HERE
  void ml_iterate();

private:

  vec tmpvecK;
  vec tmpvecD;
  vec acc_loglhood_K;

  Array<vec> acc_means;
  Array<vec> acc_covs;

  double * c_tmpvecK;
  double * c_tmpvecD;
  double * c_acc_loglhood_K;

  double ** c_acc_means;
  double ** c_acc_covs;


};

//
// convenience functions

/*!
  \ingroup MOG
  \author Conrad Sanderson

  Maximum Likelihood Expectation Maximisation based optimisation of the
  parameters of an instance of the MOG_diag class. The seed values
  (starting points) are typically first obtained via MOG_diag_kmeans().
  See [CSB06] and the references therein for detailed mathematical descriptions.

  - [CSB06]
    <a href="http://ieeexplore.ieee.org/xpl/abs_free.jsp?arNumber=1561601">
    F. Cardinaux, C. Sanderson and S. Bengio,
    "User authentication via adapted statistical models of face images",
    IEEE Transactions on Signal Processing, Vol 54, No. 1, 2006, pp. 361-373.
    </a>

  \param model_in The model to optimise (MOG_diag)
  \param X_in The training data (array of vectors)
  \param max_iter_in Maximum number of iterations. Default is 10.
  \param var_floor_in Variance floor (lowest allowable variance).  Default is 0.0 (but see the note below)
  \param weight_floor_in  Weight floor (lowest allowable weight).  Default is 0.0 (but see the note below)
  \param verbose_in  Whether progress in printed.  Default is false.

  \note The variance and weight floors are set to std::numeric_limits<double>::min()
        if they are below that value. As such, they are machine dependant.
        The largest allowable weight floor is 1/K, where K is the number of Gaussians.
*/
void MOG_diag_ML(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in = 10, double var_floor_in = 0.0, double weight_floor_in = 0.0, bool verbose_in = false);

/*!
  NOT YET IMPLEMENTED.
  Maximum a Posteriori (MAP) Expectation Maximisation optimiser for Mixtures of Gaussians.

  \param model_in The model to optimise (MOG_diag)
  \param prior_model_in The model representing the prior
  \param X_in The training data (array of vectors)
  \param max_iter_in Maximum number of iterations
  \param alpha_in Coefficient for combining the parameters with the prior.  0 <= _alpha <= 1.
  \param var_floor_in Variance floor (lowest allowable variance).  Set to 0.0 to use the default.
  \param weight_floor_in  Weight floor (lowest allowable weight).  Set to 0.0 to use the default.
  \param verbose_in ADD DOCUMENTATION HERE

  \note NOT YET IMPLEMENTED.
  \note The variance and weight floors are set to std::numeric_limits<double>::min()
        if they are below that value.
        The largest allowable weight floor is 1/K, where K is the number of Gaussians.
*/
void MOG_diag_MAP(MOG_diag &model_in, MOG_diag &prior_model_in, Array<vec> &X_in, int max_iter_in = 10, double alpha_in = 0.5, double var_floor_in = 0.0, double weight_floor_in = 0.0, bool verbose_in = false);

}

#endif // #ifndef MOG_DIAG_EM_H

