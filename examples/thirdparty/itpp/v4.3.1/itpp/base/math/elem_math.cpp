/*!
 * \file
 * \brief Elementary mathematical functions - source file
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/base/math/elem_math.h>
#include <itpp/base/itcompat.h>


namespace itpp
{

vec sqr(const cvec &data)
{
  vec temp(data.length());
  for (int i = 0; i < data.length(); i++)
    temp(i) = sqr(data(i));
  return temp;
}

mat sqr(const cmat &data)
{
  mat temp(data.rows(), data.cols());
  for (int i = 0; i < temp.rows(); i++) {
    for (int j = 0; j < temp.cols(); j++) {
      temp(i, j) = sqr(data(i, j));
    }
  }
  return temp;
}

vec abs(const cvec &data)
{
  vec temp(data.length());

  for (int i = 0;i < data.length();i++)
    temp[i] = std::abs(data[i]);

  return temp;
}

mat abs(const cmat &data)
{
  mat temp(data.rows(), data.cols());

  for (int i = 0;i < temp.rows();i++) {
    for (int j = 0;j < temp.cols();j++) {
      temp(i, j) = std::abs(data(i, j));
    }
  }

  return temp;
}

// Deprecated gamma function. Will be changed to tgamma().
double gamma(double x) { return tgamma(x); }
vec gamma(const vec &x) { return apply_function<double>(tgamma, x); }
mat gamma(const mat &x) { return apply_function<double>(tgamma, x); }

// Calculates factorial coefficient for index <= 170.
double fact(int index)
{
  it_error_if(index > 170, "fact(int index): Function overflows if index > 170.");
  it_error_if(index < 0, "fact(int index): index must be non-negative integer");
  double prod = 1;
  for (int i = 1; i <= index; i++)
    prod *= static_cast<double>(i);
  return prod;
}

// Calculates binomial coefficient "n over k".
double binom(int n, int k)
{
  it_assert(k <= n, "binom(n, k): k can not be larger than n");
  it_assert((n >= 0) && (k >= 0), "binom(n, k): n and k must be non-negative integers");
  k = ((n - k) < k) ? n - k : k;

  double out = 1.0;
  for (int i = 1; i <= k; ++i) {
    out *= (i + n - k);
    out /= i;
  }
  return out;
}

// Calculates binomial coefficient "n over k".
int binom_i(int n, int k)
{
  it_assert(k <= n, "binom_i(n, k): k can not be larger than n");
  it_assert((n >= 0) && (k >= 0), "binom_i(n, k): n and k must be non-negative integers");
  k = ((n - k) < k) ? n - k : k;

  int out = 1;
  for (int i = 1; i <= k; ++i) {
    out *= (i + n - k);
    out /= i;
  }
  return out;
}

// Calculates the base 10-logarithm of the binomial coefficient "n over k".
double log_binom(int n, int k)
{
  it_assert(k <= n, "log_binom(n, k): k can not be larger than n");
  it_assert((n >= 0) && (k >= 0), "log_binom(n, k): n and k must be non-negative integers");
  k = ((n - k) < k) ? n - k : k;

  double out = 0.0;
  for (int i = 1; i <= k; i++)
    out += log10(static_cast<double>(i + n - k))
           - log10(static_cast<double>(i));

  return out;
}

// Calculates the greatest common divisor
int gcd(int a, int b)
{
  it_assert((a >= 0) && (b >= 0), "gcd(a, b): a and b must be non-negative integers");
  int v, u, t, q;

  u = a;
  v = b;
  while (v > 0) {
    q = u / v;
    t = u - v * q;
    u = v;
    v = t;
  }
  return u;
}


vec real(const cvec &data)
{
  vec temp(data.length());

  for (int i = 0;i < data.length();i++)
    temp[i] = data[i].real();

  return temp;
}

mat real(const cmat &data)
{
  mat temp(data.rows(), data.cols());

  for (int i = 0;i < temp.rows();i++) {
    for (int j = 0;j < temp.cols();j++) {
      temp(i, j) = data(i, j).real();
    }
  }

  return temp;
}

vec imag(const cvec &data)
{
  vec temp(data.length());

  for (int i = 0;i < data.length();i++)
    temp[i] = data[i].imag();
  return temp;
}

mat imag(const cmat &data)
{
  mat temp(data.rows(), data.cols());

  for (int i = 0;i < temp.rows();i++) {
    for (int j = 0;j < temp.cols();j++) {
      temp(i, j) = data(i, j).imag();
    }
  }

  return temp;
}

vec arg(const cvec &data)
{
  vec temp(data.length());

  for (int i = 0;i < data.length();i++)
    temp[i] = std::arg(data[i]);

  return temp;
}

mat arg(const cmat &data)
{
  mat temp(data.rows(), data.cols());

  for (int i = 0;i < temp.rows();i++) {
    for (int j = 0;j < temp.cols();j++) {
      temp(i, j) = std::arg(data(i, j));
    }
  }

  return temp;
}

#ifdef _MSC_VER
cvec conj(const cvec &x)
{
  cvec temp(x.size());

  for (int i = 0; i < x.size(); i++) {
    temp(i) = std::conj(x(i));
  }

  return temp;
}

cmat conj(const cmat &x)
{
  cmat temp(x.rows(), x.cols());

  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      temp(i, j) = std::conj(x(i, j));
    }
  }

  return temp;
}
#endif

} // namespace itpp
