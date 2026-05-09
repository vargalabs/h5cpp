/*!
 * \file
 * \brief Miscellaneous statistics functions and classes - header file
 * \author Tony Ottosson, Johan Bergman and Adam Piatyszek
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

#ifndef MISC_STAT_H
#define MISC_STAT_H

#include <itpp/base/math/min_max.h>
#include <itpp/base/mat.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/matfunc.h>
#include <itpp/itexports.h>


namespace itpp
{

//! \addtogroup statistics
//!@{

/*!
  \brief A class for sampling a signal and calculating statistics
*/
class ITPP_EXPORT Stat
{
public:
  //! Default constructor
  Stat() {clear();}
  //! Destructor
  virtual ~Stat() {}

  //! Clear statistics
  virtual void clear() {
    _n_overflows = 0;
    _n_samples = 0;
    _n_zeros = 0;
    _max = 0.0;
    _min = 0.0;
    _sqr_sum = 0.0;
    _sum = 0.0;
  }

  //! Register a sample and flag for overflow
  virtual void sample(const double s, const bool overflow = false) {
    _n_samples++;
    _sum += s;
    _sqr_sum += s * s;
    if (s < _min) _min = s;
    if (s > _max) _max = s;
    if (overflow) _n_overflows++;
    if (s == 0.0) _n_zeros++;
  }

  //! Number of reported overflows
  int n_overflows() const {return _n_overflows;}
  //! Number of samples
  int n_samples() const {return _n_samples;}
  //! Number of zero samples
  int n_zeros() const {return _n_zeros;}
  //! Average over all samples
  double avg() const {return _sum / _n_samples;}
  //! Maximum sample
  double max() const {return _max;}
  //! Minimum sample
  double min() const {return _min;}
  //! Standard deviation of all samples
  double sigma() const {
    double sigma2 = _sqr_sum / _n_samples - avg() * avg();
    return std::sqrt(sigma2 < 0 ? 0 : sigma2);
  }
  //! Squared sum of all samples
  double sqr_sum() const {return _sqr_sum;}
  //! Sum of all samples
  double sum() const {return _sum;}
  //! Histogram over all samples (not implemented yet)
  vec histogram() const {return vec(0);}

protected:
  //! Number of reported overflows
  int _n_overflows;
  //! Number of samples
  int _n_samples;
  //! Number of zero samples
  int _n_zeros;
  //! Maximum sample
  double _max;
  //! Minimum sample
  double _min;
  //! Squared sum of all samples
  double _sqr_sum;
  //! Sum of all samples
  double _sum;
};


//! The mean value
ITPP_EXPORT double mean(const vec &v);
//! The mean value
ITPP_EXPORT std::complex<double> mean(const cvec &v);
//! The mean value
ITPP_EXPORT double mean(const svec &v);
//! The mean value
ITPP_EXPORT double mean(const ivec &v);
//! The mean value
ITPP_EXPORT double mean(const mat &m);
//! The mean value
ITPP_EXPORT std::complex<double> mean(const cmat &m);
//! The mean value
ITPP_EXPORT double mean(const smat &m);
//! The mean value
ITPP_EXPORT double mean(const imat &m);

//! The geometric mean of a vector
template<class T>
double geometric_mean(const Vec<T> &v)
{
  return std::exp(std::log(static_cast<double>(prod(v))) / v.length());
}

//! The geometric mean of a matrix
template<class T>
double geometric_mean(const Mat<T> &m)
{
  return std::exp(std::log(static_cast<double>(prod(prod(m))))
                  / (m.rows() * m.cols()));
}

//! The median
template<class T>
double median(const Vec<T> &v)
{
  Vec<T> invect(v);
  sort(invect);
  return (double)(invect[(invect.length()-1)/2] + invect[invect.length()/2]) / 2.0;
}

//! Calculate the 2-norm: norm(v)=sqrt(sum(abs(v).^2))
ITPP_EXPORT double norm(const cvec &v);

//! Calculate the 2-norm: norm(v)=sqrt(sum(abs(v).^2))
template<class T>
double norm(const Vec<T> &v)
{
  double E = 0.0;
  for (int i = 0; i < v.size(); i++)
    E += sqr(static_cast<double>(v[i]));

  return std::sqrt(E);
}

//! Calculate the p-norm: norm(v,p)=sum(abs(v).^2)^(1/p)
ITPP_EXPORT double norm(const cvec &v, int p);

//! Calculate the p-norm: norm(v,p)=sum(abs(v).^2)^(1/p)
template<class T>
double norm(const Vec<T> &v, int p)
{
  double E = 0.0;
  for (int i = 0; i < v.size(); i++)
    E += std::pow(fabs(static_cast<double>(v[i])), static_cast<double>(p));

  return std::pow(E, 1.0 / p);
}

//! Calculate the Frobenius norm for s = "fro" (equal to 2-norm)
ITPP_EXPORT double norm(const cvec &v, const std::string &s);

//! Calculate the Frobenius norm for s = "fro" (equal to 2-norm)
template<class T>
double norm(const Vec<T> &v, const std::string &s)
{
  it_assert(s == "fro", "norm(): Unrecognised norm");

  double E = 0.0;
  for (int i = 0; i < v.size(); i++)
    E += sqr(static_cast<double>(v[i]));

  return std::sqrt(E);
}

/*!
 * Calculate the p-norm of a real matrix
 *
 * p = 1: max(svd(m))
 * p = 2: max(sum(abs(X)))
 *
 * Default if no p is given is the 2-norm
 */
ITPP_EXPORT double norm(const mat &m, int p = 2);

/*!
 * Calculate the p-norm of a complex matrix
 *
 * p = 1: max(svd(m))
 * p = 2: max(sum(abs(X)))
 *
 * Default if no p is given is the 2-norm
 */
ITPP_EXPORT double norm(const cmat &m, int p = 2);

//! Calculate the Frobenius norm of a matrix for s = "fro"
ITPP_EXPORT double norm(const mat &m, const std::string &s);

//! Calculate the Frobenius norm of a matrix for s = "fro"
ITPP_EXPORT double norm(const cmat &m, const std::string &s);


//! The variance of the elements in the vector. Normalized with N-1 to be unbiased.
ITPP_EXPORT double variance(const cvec &v);

//! The variance of the elements in the vector. Normalized with N-1 to be unbiased.
template<class T>
double variance(const Vec<T> &v)
{
  int len = v.size();
  const T *p = v._data();
  double sum = 0.0, sq_sum = 0.0;

  for (int i = 0; i < len; i++, p++) {
    sum += *p;
    sq_sum += *p * *p;
  }

  return (double)(sq_sum - sum*sum / len) / (len - 1);
}

//! Calculate the energy: squared 2-norm. energy(v)=sum(abs(v).^2)
template<class T>
double energy(const Vec<T> &v)
{
  return sqr(norm(v));
}


//! Return true if the input value \c x is within the tolerance \c tol of the reference value \c xref
inline bool within_tolerance(double x, double xref, double tol = 1e-14)
{
  return (fabs(x -xref) <= tol) ? true : false;
}

//! Return true if the input value \c x is within the tolerance \c tol of the reference value \c xref
inline bool within_tolerance(std::complex<double> x, std::complex<double> xref, double tol = 1e-14)
{
  return (abs(x -xref) <= tol) ? true : false;
}

//! Return true if the input vector \c x is elementwise within the tolerance \c tol of the reference vector \c xref
inline bool within_tolerance(const vec &x, const vec &xref, double tol = 1e-14)
{
  return (max(abs(x -xref)) <= tol) ? true : false;
}

//! Return true if the input vector \c x is elementwise within the tolerance \c tol of the reference vector \c xref
inline bool within_tolerance(const cvec &x, const cvec &xref, double tol = 1e-14)
{
  return (max(abs(x -xref)) <= tol) ? true : false;
}

//! Return true if the input matrix \c X is elementwise within the tolerance \c tol of the reference matrix \c Xref
inline bool within_tolerance(const mat &X, const mat &Xref, double tol = 1e-14)
{
  return (max(max(abs(X -Xref))) <= tol) ? true : false;
}

//! Return true if the input matrix \c X is elementwise within the tolerance \c tol of the reference matrix \c Xref
inline bool within_tolerance(const cmat &X, const cmat &Xref, double tol = 1e-14)
{
  return (max(max(abs(X -Xref))) <= tol) ? true : false;
}

/*!
  \brief Calculate the central moment of vector x

  The \f$r\f$th sample central moment of the samples in the vector
  \f$ \mathbf{x} \f$ is defined as

  \f[
  m_r = \mathrm{E}[x-\mu]^r = \frac{1}{n} \sum_{i=0}^{n-1} (x_i - \mu)^r
  \f]
  where \f$\mu\f$ is the sample mean.
*/
ITPP_EXPORT double moment(const vec &x, const int r);

/*!
  \brief Calculate the skewness excess of the input vector x

  The skewness is a measure of the degree of asymmetry of distribution. Negative
  skewness means that the distribution is spread more to the left of the mean than to
  the right, and vice versa if the skewness is positive.

  The skewness of the samples in the vector \f$ \mathbf{x} \f$ is
  \f[
  \gamma_1 = \frac{\mathrm{E}[x-\mu]^3}{\sigma^3}
  \f]
  where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.

  The skewness is estimated as
  \f[
  \gamma_1 = \frac{k_3}{{k_2}^{3/2}}
  \f]
  where
  \f[
  k_2 = \frac{n}{n-1} m_2
  \f]
  and
  \f[
  k_3 = \frac{n^2}{(n-1)(n-2)} m_3
  \f]
  Here \f$m_2\f$ is the sample variance and \f$m_3\f$ is the 3rd sample
  central moment.
*/
ITPP_EXPORT double skewness(const vec &x);


/*!
  \brief Calculate the kurtosis excess of the input vector x

  The kurtosis excess is a measure of peakedness of a distribution.
The kurtosis excess is defined as
  \f[
  \gamma_2 = \frac{\mathrm{E}[x-\mu]^4}{\sigma^4} - 3
  \f]
  where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.

  The kurtosis excess is estimated as
  \f[
  \gamma_2 = \frac{k_4}{{k_2}^2}
  \f]
  where
  \f[
  k_2 = \frac{n}{n-1} m_2
  \f]
  and
  \f[
  k_4 = \frac{n^2 [(n+1)m_4 - 3(n-1){m_2}^2]}{(n-1)(n-2)(n-3)}
  \f]
  Here \f$m_2\f$ is the sample variance and \f$m_4\f$ is the 4th sample
  central moment.
*/
ITPP_EXPORT double kurtosisexcess(const vec &x);

/*!
\brief Calculate the kurtosis of the input vector x

The kurtosis is a measure of peakedness of a distribution. The kurtosis
is defined as
\f[
\gamma_2 = \frac{\mathrm{E}[x-\mu]^4}{\sigma^4}
\f]
where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.
For a Gaussian variable, the kurtusis is 3.

See also the definition of kurtosisexcess.
*/
inline double kurtosis(const vec &x) {return kurtosisexcess(x) + 3;}

//!@}

} // namespace itpp

#endif // #ifndef MISC_STAT_H
