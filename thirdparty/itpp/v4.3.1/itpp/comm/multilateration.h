/*!
 * \file
 * \brief Definition of multilateration class for indoor localization
 * \author Bogdan Cristea
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

#ifndef MULTILATERATION_H
#define MULTILATERATION_H

#include <itpp/itbase.h>
#include <itpp/itexports.h>

namespace itpp
{

class Algorithm;
struct Point;

/*!
  \ingroup misccommfunc
  \brief %Multilateration class for 3D indoor localization

  Implements geometry-based methods for indoor localization:
  - spherical multilateration (which uses Time Of Arrival (TOA) ranging techniques),
  - hyperbolic multilateration (using Time Difference Of Arrival (TDOA)) and
  - hybrid multilateration (both TOA and TDOA are used)

  In addition, it allows to compute the theoretical performance of the algorithm based on Cramer Rao Lower Bound (CRLB).

  Geometry-based methods for indoor localization use several Base Stations (BSs), whose position is known, in order
  to compute the position of the Mobile Station (MS). By computing the distance to each BS (TOA ranging) or the difference
  between the distances to two BSs (TDOA) a system of non-linear equations is obtained that allows to compute the MS
  position. At least 4 measures (TOA or TDOA) are needed in order to obtain a determinate equation system.
  The algorithm implemented in this class can handle any number of measures (at least four) by using an asymptotic Maximum
  Likelihood (ML) estimator [1]. The input of the algorithm is represented by a method vector, specifying the type of each
  ranging measure (0 for TOA and 1 for TDOA), a matrix with BSs positions and a vector (for spherical and hybrid) or a matrix
  (for hyperbolic multilateration) with the ranging measures. The output is a vector of length 3 with the position of the MS
  in 3D cartezian coordinates.

  Note that for hybrid multilateration the method vector should have at least a one and a zero, for spherical multilateration the
  method vector is all zeros, while for hyperbolic multilateration is all ones.

  The CRLB is computed as the Euclidean distance between the estimated position of the MS and the true MS position. The noise
  variance is needed as input together with the true MS position. It is assumed that the noise affecting the measures has
  the same variance for all measures.

  Usage example:
  \code
  Multilateration multi;
  bvec method(4);
  method.zeros();//spherical multilateration
  mat bs_pos = randn(3, 4);//four BSs
  multi.setup(method, bs_pos);
  vec measures(4);
  //measures are generated following TOA ranging method (see unit tests file for an example)
  vec ms_pos;//algorithm output
  bool rc = multi.get_pos(ms_pos, measures);//algorithm has failed if the output is false
  \endcode

  Reference:
  [1] Urruela, A. and Riba, J. - Novel closed-form ML position estimator for hyperbolic location, ICASSP'04
*/
class ITPP_EXPORT Multilateration
{
public:
  //! %Multilateration types as detected from user input (method binary vector)
  enum Type {MULTI_FAILURE = -1, //!< the algorithm has failed
             MULTI_SPHERICAL, //!< spherical multilateration
             MULTI_HYPERBOLIC, //!< hyperbolic multilateration
             MULTI_HYBRID //!< hybrid multilateration
            };
  //! %Multilateration class default constructor
  Multilateration() :
    algo_(NULL), nb_fails_part(0), nb_fails_pos(0), type_(MULTI_FAILURE), method_(itpp::bvec()), bs_pos_(NULL), nb_bs_(0)
    {}
  //! %Multilateration class constructor
  /*! The BS positions are specified as a matrix, each BS position can be specified on either rows or columns.
   * The method vector specify the measure type: 0 for TOA or 1 for TDOA. For example, a vector with all zeros represents
   * spherical multilateration, while a vector with all ones represents hyperbolic multilateration.
   *
   * For spherical multilateration the number of BSs and the method length must be equal.
   * For hybrid and hyperbolic multilateration the number of BSs must be the method length plus one.
   */
  Multilateration(const itpp::bvec &method, //!< multilateration method
                  const itpp::mat &bs_pos //!< base station positions in 3D cartezian coordinates
                 ) :
    algo_(NULL), nb_fails_part(0), nb_fails_pos(0), type_(MULTI_FAILURE), method_(itpp::bvec()), bs_pos_(NULL), nb_bs_(0) {
    setup(method, bs_pos);
  }
  //! %Multilateration destructor
  virtual ~Multilateration();
  //! Setup function for specifying the multilateration method and the base station positions
  /*! The BS positions are specified as a matrix, each BS position can be specified on either rows or columns.
   * The method vector specify the measure type: O for TOA or 1 for TDOA. A vector with all zeros represents
   * spherical multilateration, while a vector with all ones represents hyperbolic multilateration.
   *
   * For spherical multilateration the number of BSs and the method lenght must be equal and it must also equal
   * the length of the measures vector.
   *
   * For hybrid and hyperbolic multilateration the number of BSs must be the method length plus one.
   */
  void setup(const itpp::bvec &method, //!< multilateration method
             const itpp::mat &bs_pos //!< base station positions
            ) {
    if((false == set_bs_pos(bs_pos)) || (false == set_method(method))) {
      it_error("cannot init multilateration");
    }
  }
  //! Computes the mobile station position for spherical and hybrid multilateration
  /*! For spherical multilateration the vector of measures should be generated as follows:
   * \f[ measures(i) = dist(bs\_pos(i), ms\_pos) \f]
   * where \f$ dist() \f$ is the Euclidean distance between two points in 3D cartezian coordinates.
   *
   * For hybrid multilateration the vector of measures is generated as:
   * - if \f$ 1 == method(i) \f$ (TDOA ranging measure)
   *   \f[ measures(i) = dist(bs\_pos(i+1), ms\_pos)-dist(bs\_pos(0), ms\_pos) \f]
   * - if \f$ 0 == method(i) \f$ (TOA ranging measure)
   *   \f[ measures(i) = dist(bs\_pos(i), ms\_pos) \f]
   */
  bool get_pos(itpp::vec &ms_pos, //!< output with mobile station position in 3D cartezian coordinates
               const itpp::vec &measures //!< vector with ranging measures
              ) {
    return get_pos(ms_pos, measures._data());
  }
  //! Computes the mobile station position for hyperbolic multilateration
  /*! The matrix of measures is computed as follows:
  * \f[ measures(i,j) = dist(bs\_pos(i), ms\_pos)-dist(bs\_pos(j), ms\_pos) \f]
  * where \f$ dist() \f$ is the Euclidean distance between two points in 3D cartezian coordinates.
  */
  bool get_pos(itpp::vec &ms_pos, //!< output with mobile station position in 3D cartezian coordiates
               const itpp::mat &measures //!< matrix with ranging measures
              ) {
    return get_pos(ms_pos, measures._data());
  }
  //! Gets the number of failures of the partitioning algorithm used internally by the ML-estimator
  unsigned int get_nb_fails_part() const {
    return nb_fails_part;
  }
  //! Gets the number of failures of the positioning algorithm used internally by the ML-estimator
  unsigned int get_nb_fails_pos() const {
    return nb_fails_pos;
  }
  //! Resets the error counters (number of failures for the partitioning and positioning algorithms)
  void reset_err_counters() {
    nb_fails_part = 0;
    nb_fails_pos = 0;
  }
  //! Gets the type of the multilateration method currently used by the ML-estimator
  Type get_type() const {
    return type_;
  }
  //! Computes the Cramer Rao lower bound for the ML-estimator assuming the same noise variance for all measures
  double get_crlb(const vec &ms_pos, //!< true mobile station position
                  double sigma2 //!< noise variance affecting the measures
                 );
private:
  //! Computes MS position using as input an arrays of measures
  bool get_pos(itpp::vec &ms_pos, const double *measures);
  //! Sets multilateation method vector
  bool set_method(const itpp::bvec &method);
  //! Sets BS positions (must be called before set_method)
  bool set_bs_pos(const itpp::mat &bs_pos);
  //! Converts a hybrid multilateration into an equivalent multilateration
  bool hybrid2spherical(Point *bs_pos, double *meas);
  bool partition(unsigned int **subsets_idx, unsigned int *subsets_nb, const Point *bs_pos, unsigned int nb_bs, unsigned int subset_len);
  //! Computes MS position using the ML-estimator
  bool get_ml_pos(Point *ms_pos, const Point *bs_pos, unsigned int nb_bs, const unsigned int *subsets_idx, unsigned int subsets_nb, unsigned int subset_len);
  //! Gets a subset of BSs based on an input subset of indices in the initial array of BSs
  bool get_bs_pos_subset(Point *bs_pos_subset, const Point *bs_pos, unsigned int nb_bs, const unsigned int *subset_idx, unsigned int subset_len);
  //! Computes the product \f$A^T*d*A\f$
  bool prod(double *out, const double *AT, const unsigned int *d, unsigned int cols, unsigned int rows);
  Algorithm *algo_;
  unsigned int nb_fails_part;
  unsigned int nb_fails_pos;
  Type type_;
  itpp::bvec method_;
  Point *bs_pos_;
  unsigned int nb_bs_;
};

}

#endif
