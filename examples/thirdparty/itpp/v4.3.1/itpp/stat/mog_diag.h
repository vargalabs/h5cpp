/*!
 * \file
 * \brief Diagonal Mixture of Gaussians class - header file
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

#ifndef MOG_DIAG_H
#define MOG_DIAG_H

#include <itpp/stat/mog_generic.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \ingroup MOG
  \brief Diagonal Mixture of Gaussians (MOG) class
  \author Conrad Sanderson

  Used for representing a statistical distribution as a
  convex combination of multi-variate Gaussian functions.
  Also known as a Gaussian Mixture Model.
  This class allows loading and saving of the MOG's parameters,
  as well as calculation of likelihoods. The parameters
  are set by the user or an optimisation algorithm
  (for example, see the MOG_diag_EM class).

  \note This class is optimised for diagonal covariance matrices.
        For speed reasons it uses C style arrays for direct access to memory.
*/
class ITPP_EXPORT MOG_diag : public MOG_generic
{

public:

  /*! \brief Default constructor
      \note An empty model is created.
            The likelihood functions are not useable
            until the model's parameters are set
  */
  MOG_diag() { zero_all_ptrs(); init(); }

  /*! \brief Construct the MOG_diag object by loading the parameters from a model file
      \param name The model's filename
  */
  MOG_diag(const std::string &name) { zero_all_ptrs(); load(name); }

  /*! \brief construct a default model (all Gaussians have zero mean and unit variance for all dimensions)
      \param K_in Number of Gaussians
      \param D_in Dimensionality
      \param full_in Ignored.  Present for compatability with the MOG_generic class
  */
  MOG_diag(const int &K_in, const int &D_in, bool full_in = false) { zero_all_ptrs(); init(K_in, D_in, full_in); }

  /*! \brief Construct a model using user supplied mean vectors
      \param means_in Array of mean vectors
      \note  The number of mean vectors specifies the number of Gaussians.
      The covariance matrices are in effect set equal to the identity matrix.
      The weights for all Gaussians are the same, equal to 1/K, where K is the number of Gaussians
  */
  MOG_diag(Array<vec> &means_in, bool) { zero_all_ptrs(); init(means_in, false);  }

  /*! \brief Construct a model using user supplied parameters (diagonal covariance version)
      \param means_in Array of mean vectors
      \param diag_covs_in Array of vectors representing diagonal covariances
      \param weights_in vector of weights
      \note  The number of mean vectors, covariance vectors and weights must be the same
  */
  MOG_diag(Array<vec> &means_in, Array<vec> &diag_covs_in, vec &weights_in) { zero_all_ptrs(); init(means_in, diag_covs_in, weights_in); }

  /*! \brief Construct a model using user supplied parameters (full covariance version)
      \param means_in Array of mean vectors
      \param full_covs_in Array of full covariance matrices
      \param weights_in vector of weights
      \note  The full covariance matrices are converted to be diagonal.
             The number of mean vectors, covariance matrices and weights must be the same.
  */
  MOG_diag(Array<vec> &means_in, Array<mat> &full_covs_in, vec &weights_in) { zero_all_ptrs(); init(means_in, full_covs_in, weights_in); convert_to_diag(); }

  //! Default destructor
  ~MOG_diag() { cleanup(); }

  /*! \brief Release memory used by the model. The model will be empty.
      \note The likelihood functions are not useable
            until the model's parameters are re-initialised
  */
  void cleanup() { free_all_ptrs(); MOG_generic::cleanup(); }

  /*! \brief Initialise the model by loading the parameters from a model file.
      \param name_in The model's filename
      \note If the model file contains a full covariance matrix model,
            the covariance matrices will be converted to be diagonal after loading.
  */
  void load(const std::string &name_in);

  //! Do nothing.  Present for compatability with the MOG_generic class.
  void convert_to_full() {};

  //! calculate the log likelihood of C vector \c c_x_in using only Gaussian \c k
  double log_lhood_single_gaus(const double * c_x_in, const int k) const;

  //! calculate the log likelihood of IT++ vector \c x_in using only Gaussian \c k
  double log_lhood_single_gaus(const vec &x_in, const int k) const;

  //! calculate the log likelihood of C vector \c c_x_in
  double log_lhood(const double * c_x_in);

  //! calculate the log likelihood of IT++ vector \c x_in
  double log_lhood(const vec &x_in);

  //! calculate the likelihood of C vector \c c_x_in
  double lhood(const double * c_x_in);

  //! calculate the likelihood of IT++ vector \c x_in
  double lhood(const vec &x_in);

  //! calculate the average log likelihood of an array of C vectors ( \c c_x_in )
  double avg_log_lhood(const double ** c_x_in, int N);

  //! calculate the average log likelihood of an array of IT++ vectors ( \c X_in )
  double avg_log_lhood(const Array<vec> & X_in);

protected:

  void setup_means();
  void setup_covs();
  void setup_weights();
  void setup_misc();

  //! ADD DOCUMENTATION HERE
  double log_lhood_single_gaus_internal(const double * c_x_in, const int k) const;
  //! ADD DOCUMENTATION HERE
  double log_lhood_single_gaus_internal(const vec &x_in, const int k) const;
  //! ADD DOCUMENTATION HERE
  double log_lhood_internal(const double * c_x_in);
  //! ADD DOCUMENTATION HERE
  double log_lhood_internal(const vec &x_in);
  //! ADD DOCUMENTATION HERE
  double lhood_internal(const double * c_x_in);
  //! ADD DOCUMENTATION HERE
  double lhood_internal(const vec &x_in);

  //! Enable C style access to an Array of vectors (vec)
  double ** enable_c_access(Array<vec> & A_in);

  //! Enable C style access to an Array of vectors (ivec)
  int ** enable_c_access(Array<ivec> & A_in);

  //! Enable C style access to a vector (vec)
  double * enable_c_access(vec & v_in);

  //! Enable C style access to a vector (ivec)
  int * enable_c_access(ivec & v_in);

  //! Disable C style access to an Array of vectors (vec)
  double ** disable_c_access(double ** A_in);

  //! Disable C style access to an Array of vectors (ivec)
  int ** disable_c_access(int ** A_in);

  //! Disable C style access to a vector (vec)
  double * disable_c_access(double * v_in);

  //! Disable C style access to a vector (ivec)
  int * disable_c_access(int * v_in);

  //! ADD DOCUMENTATION HERE
  void zero_all_ptrs();
  //! ADD DOCUMENTATION HERE
  void free_all_ptrs();

  //! pointers to the mean vectors
  double ** c_means;

  //! pointers to the covariance vectors
  double ** c_diag_covs;

  //! pointers to the inverted covariance vectors
  double ** c_diag_covs_inv_etc;

  //! pointer to the weight vector
  double * c_weights;

  //! pointer to the log version of the weight vector
  double * c_log_weights;

  //! pointer to the log_det_etc vector
  double * c_log_det_etc;

private:

  vec tmpvecK;
  double * c_tmpvecK;

};

}

#endif // #ifndef MOG_DIAG_H

