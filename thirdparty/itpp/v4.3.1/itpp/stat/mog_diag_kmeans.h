/*!
 * \file
 * \brief K-means based optimiser for Mixture of Gaussians - header file
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


#ifndef MOG_DIAG_KMEANS_H
#define MOG_DIAG_KMEANS_H

#include <itpp/stat/mog_diag.h>
#include <itpp/itexports.h>
#include <itpp/base/base_exports.h>

namespace itpp
{

/*!
  \brief support class for MOG_diag_kmeans()
  \author Conrad Sanderson
*/
class ITPP_EXPORT MOG_diag_kmeans_sup : public MOG_diag
{

public:
  //! Default constructor
  MOG_diag_kmeans_sup() { verbose = false; }

  //! Default destructor
  ~MOG_diag_kmeans_sup() { }

  //! ADD DOCUMENTATION HERE
  void run(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in = 10, double trust_in = 0.5, bool normalise_in = true, bool verbose_in = false);

protected:

  //! ADD DOCUMENTATION HERE
  inline double dist(const double * x, const double * y) const;
  //! ADD DOCUMENTATION HERE
  void assign_to_means();
  //! ADD DOCUMENTATION HERE
  void recalculate_means();
  //! ADD DOCUMENTATION HERE
  bool dezombify_means();
  //! ADD DOCUMENTATION HERE
  double measure_change() const;
  //! ADD DOCUMENTATION HERE
  void initial_means();
  //! ADD DOCUMENTATION HERE
  void iterate();
  //! ADD DOCUMENTATION HERE
  void calc_means();
  //! ADD DOCUMENTATION HERE
  void calc_covs();
  //! ADD DOCUMENTATION HERE
  void calc_weights();
  //! ADD DOCUMENTATION HERE
  void normalise_vectors();
  //! ADD DOCUMENTATION HERE
  void unnormalise_vectors();
  //! ADD DOCUMENTATION HERE
  void unnormalise_means();

  //! Maximum number of iterations
  int max_iter;

  /*! \brief trust factor, where 0 <= trust <= 1.
      \note The higher the trust factor, the more we trust the estimates of covariance matrices and weights.
   */
  double trust;

  //! Whether we print the progress
  bool verbose;

  //! number of training vectors
  int N;

  //! 'C' pointers to training vectors
  double ** c_X;

  //! means from the previous iteration, used to measure progress
  Array<vec> means_old;

  //! 'C' pointers to old means
  double ** c_means_old;

  //! contains indices of vectors assigned to each mean
  Array<ivec> partitions;

  //! 'C' pointers to partition vectors
  int ** c_partitions;

  //! keeps a count of the number of vectors assigned to each mean
  ivec count;

  //! 'C' pointer to the count vector
  int * c_count;

private:

  vec norm_mu;
  double * c_norm_mu;

  vec norm_sd;
  double * c_norm_sd;

  vec tmpvec;
  double * c_tmpvec;


};

//
// convenience functions

/*!
  \ingroup MOG
  \author Conrad Sanderson

  K-means based optimisation (training) of the parameters of an instance of the MOG_diag class.
  The obtained parameters are typically used as a seed by MOG_diag_ML().

  \param model_in The model to optimise
  \param X_in The training data
  \param max_iter_in Maximum number of iterations. Default is 10.
  \param trust_in The trust factor, where 0 <= \c trust_in <= 1.  Default is 0.5.
  \param normalise_in Use normalised distance measure (in effect). Default is true.
  \param verbose_in Whether to print progress. Default is false.

  \note The higher the trust factor, the more we trust
  the estimates of covariance matrices and weights.
  Set this to 1.0 only if you have plenty of training data.
  One rule of thumb is to have 10*D vectors per Gaussian,
  where D is the dimensionality of the vectors.
  For smaller amounts of data, a lower trust factor
  will help (but not completely avoid) the EM algorithm
  ( used in MOG_diag_ML() ) from getting stuck in a local minimum.

  \note Setting \c normalise_in to true causes the the training
  data to be normalised to zero mean and unit variance prior
  to running the k-means algorithm.  The data is unnormalised
  before returning.  The normalisation helps clustering when
  the range of values varies greatly between dimensions.
  e.g. dimension 1 may have values in the [-1,+1] interval,
  while dimension 2 may have values in the [-100,+100] interval.
  Without normalisation, the distance between vectors is
  dominated by dimension 2.
*/
void MOG_diag_kmeans(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in = 10, double trust_in = 0.5, bool normalise_in = true, bool verbose_in = false);

}

#endif  // #ifndef MOG_DIAG_KMEANS_H

