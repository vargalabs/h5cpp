/*!
 * \file
 * \brief Miscellaneous statistics functions and classes - source file
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

#include <itpp/base/algebra/svd.h>
#include <itpp/stat/misc_stat.h>


namespace itpp
{

double mean(const vec &v)
{
  return sum(v) / v.length();
}

std::complex<double> mean(const cvec &v)
{
  return sum(v) / double(v.size());
}

double mean(const svec &v)
{
  return (double)sum(v) / v.length();
}

double mean(const ivec &v)
{
  return (double)sum(v) / v.length();
}

double mean(const mat &m)
{
  return sum(sum(m)) / (m.rows()*m.cols());
}

std::complex<double> mean(const cmat &m)
{
  return sum(sum(m)) / static_cast<std::complex<double> >(m.rows()*m.cols());
}

double mean(const smat &m)
{
  return static_cast<double>(sum(sum(m))) / (m.rows()*m.cols());
}

double mean(const imat &m)
{
  return static_cast<double>(sum(sum(m))) / (m.rows()*m.cols());
}


double norm(const cvec &v)
{
  double E = 0.0;
  for (int i = 0; i < v.length(); i++)
    E += std::norm(v[i]);

  return std::sqrt(E);
}

double norm(const cvec &v, int p)
{
  double E = 0.0;
  for (int i = 0; i < v.size(); i++)
    E += std::pow(std::norm(v[i]), p / 2.0); // Yes, 2.0 is correct!

  return std::pow(E, 1.0 / p);
}

double norm(const cvec &v, const std::string &)
{
  return norm(v, 2);
}

/*
 * Calculate the p-norm of a real matrix
 * p = 1: max(svd(m))
 * p = 2: max(sum(abs(X)))
 */
double norm(const mat &m, int p)
{
  it_assert((p == 1) || (p == 2),
            "norm(): Can only calculate a matrix norm of order 1 or 2");

  if (p == 1)
    return max(sum(abs(m)));
  else
    return max(svd(m));
}

/*
 * Calculate the p-norm of a complex matrix
 * p = 1: max(svd(m))
 * p = 2: max(sum(abs(X)))
 */
double norm(const cmat &m, int p)
{
  it_assert((p == 1) || (p == 2),
            "norm(): Can only calculate a matrix norm of order 1 or 2");

  if (p == 1)
    return max(sum(abs(m)));
  else
    return max(svd(m));
}

// Calculate the Frobenius norm of matrix m for s = "fro"
double norm(const mat &m, const std::string &s)
{
  it_assert(s == "fro", "norm(): Unrecognised norm");
  double E = 0.0;
  for (int r = 0; r < m.rows(); ++r) {
    for (int c = 0; c < m.cols(); ++c) {
      E += m(r, c) * m(r, c);
    }
  }
  return std::sqrt(E);
}

// Calculate the Frobenius norm of matrix m for s = "fro"
double norm(const cmat &m, const std::string &s)
{
  it_assert(s == "fro", "norm(): Unrecognised norm");
  double E = 0.0;
  for (int r = 0; r < m.rows(); ++r) {
    for (int c = 0; c < m.cols(); ++c) {
      E += std::norm(m(r, c));
    }
  }
  return std::sqrt(E);
}


double variance(const cvec &v)
{
  int len = v.size();
  double sq_sum = 0.0;
  std::complex<double> sum = 0.0;
  const std::complex<double> *p = v._data();

  for (int i = 0; i < len; i++, p++) {
    sum += *p;
    sq_sum += std::norm(*p);
  }

  return (double)(sq_sum - std::norm(sum) / len) / (len - 1);
}

double moment(const vec &x, const int r)
{
  double m = mean(x), mr = 0;
  int n = x.size();
  double temp;

  switch (r) {
  case 1:
    for (int j = 0; j < n; j++)
      mr += (x(j) - m);
    break;
  case 2:
    for (int j = 0; j < n; j++)
      mr += (x(j) - m) * (x(j) - m);
    break;
  case 3:
    for (int j = 0; j < n; j++)
      mr += (x(j) - m) * (x(j) - m) * (x(j) - m);
    break;
  case 4:
    for (int j = 0; j < n; j++) {
      temp = (x(j) - m) * (x(j) - m);
      temp *= temp;
      mr += temp;
    }
    break;
  default:
    for (int j = 0; j < n; j++)
      mr += std::pow(x(j) - m, double(r));
    break;
  }

  return mr / n;
}


double skewness(const vec &x)
{
  int n = x.size();

  double k2 = variance(x) * n / (n - 1); // 2nd k-statistic
  double k3 = moment(x, 3) * n * n / (n - 1) / (n - 2); //3rd k-statistic

  return k3 / std::pow(k2, 3.0 / 2.0);
}

double kurtosisexcess(const vec &x)
{
  int n = x.size();
  double m2 = variance(x);
  double m4 = moment(x, 4);

  double k2 = m2 * n / (n - 1); // 2nd k-statistic
  double k4 = (m4 * (n + 1) - 3 * (n - 1) * m2 * m2) * n * n / (n - 1) / (n - 2) / (n - 3); //4th k-statistic

  return k4 / (k2*k2);
}

} // namespace itpp
