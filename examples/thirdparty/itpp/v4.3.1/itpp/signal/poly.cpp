/*!
 * \file
 * \brief Polynomial functions
 * \author Tony Ottosson, Kumar Appaiah and Adam Piatyszek
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

#include <itpp/base/itcompat.h>
#include <itpp/signal/poly.h>
#include <itpp/base/converters.h>
#include <itpp/base/algebra/eigen.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>


namespace itpp
{

void poly(const vec &r, vec &p)
{
  int n = r.size();

  p.set_size(n + 1, false);
  p.zeros();
  p(0) = 1.0;

  for (int i = 0; i < n; i++)
    p.set_subvector(1, p(1, i + 1) - r(i)*p(0, i));
}

void poly(const cvec &r, cvec &p)
{
  int n = r.size();

  p.set_size(n + 1, false);
  p.zeros();
  p(0) = 1.0;

  for (int i = 0; i < n; i++)
    p.set_subvector(1, p(1, i + 1) - r(i)*p(0, i));
}



void roots(const vec &p, cvec &r)
{
  int n = p.size(), m, l;
  ivec f = find(p != 0.0);
  m = f.size();
  vec v = p;
  mat A;

  if (m > 0 && n > 1) {
    v = v(f(0), f(m - 1));
    l = v.size();

    if (l > 1) {

      A = diag(ones(l - 2), -1);
      A.set_row(0, -v(1, l - 1) / v(0));
      r = eig(A);
      cvec d;
      cmat V;
      eig(A, d , V);

      if (f(m - 1) < n)
        r = concat(r, zeros_c(n - f(m - 1) - 1));
    }
    else {
      r.set_size(n - f(m - 1) - 1, false);
      r.zeros();
    }
  }
  else
    r.set_size(0, false);
}

void roots(const cvec &p, cvec &r)
{
  int n = p.size(), m, l;
  ivec f;

  // find all non-zero elements
  for (int i = 0; i < n; i++)
    if (p(i) != 0.0)
      f = concat(f, i);


  m = f.size();
  cvec v = p;
  cmat A;

  if (m > 0 && n > 1) {
    v = v(f(0), f(m - 1));
    l = v.size();

    if (l > 1) {
      A = diag(ones_c(l - 2), -1);
      A.set_row(0, -v(1, l - 1) / v(0));
      r = eig(A);
      if (f(m - 1) < n)
        r = concat(r, zeros_c(n - f(m - 1) - 1));
    }
    else {
      r.set_size(n - f(m - 1) - 1, false);
      r.zeros();
    }
  }
  else
    r.set_size(0, false);
}


vec polyval(const vec &p, const vec &x)
{
  it_error_if(p.size() == 0, "polyval: size of polynomial is zero");
  it_error_if(x.size() == 0, "polyval: size of input value vector is zero");

  vec out(x.size());

  out = p(0);

  for (int i = 1; i < p.size(); i++)
    out = p(i) + elem_mult(x, out);

  return out;
}

cvec polyval(const vec &p, const cvec &x)
{
  it_error_if(p.size() == 0, "polyval: size of polynomial is zero");
  it_error_if(x.size() == 0, "polyval: size of input value vector is zero");

  cvec out(x.size());

  out = p(0);

  for (int i = 1; i < p.size(); i++)
    out = std::complex<double>(p(i)) + elem_mult(x, out);

  return out;
}

cvec polyval(const cvec &p, const vec &x)
{
  it_error_if(p.size() == 0, "polyval: size of polynomial is zero");
  it_error_if(x.size() == 0, "polyval: size of input value vector is zero");

  cvec out(x.size());

  out = p(0);

  for (int i = 1; i < p.size(); i++)
    out = std::complex<double>(p(i)) + elem_mult(to_cvec(x), out);

  return out;
}

cvec polyval(const cvec &p, const cvec &x)
{
  it_error_if(p.size() == 0, "polyval: size of polynomial is zero");
  it_error_if(x.size() == 0, "polyval: size of input value vector is zero");

  cvec out(x.size());

  out = p(0);

  for (int i = 1; i < p.size(); i++)
    out = p(i) + elem_mult(x, out);

  return out;
}

double cheb(int n, double x)
{
  it_assert((n >= 0), "cheb(): need a non-negative order n!");

  if (x < 1.0 && x > -1.0) {
    return std::cos(n * std::acos(x));
  }
  else if (x <= -1) {
    return (is_even(n) ? std::cosh(n * ::acosh(-x))
            : -std::cosh(n * ::acosh(-x)));
  }
  return std::cosh(n * ::acosh(x));
}

vec cheb(int n, const vec &x)
{
  it_assert_debug(x.size() > 0, "cheb(): empty vector");

  vec out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out(i) = cheb(n, x(i));
  }
  return out;
}

mat cheb(int n, const mat &x)
{
  it_assert_debug((x.rows() > 0) && (x.cols() > 0), "cheb(): empty matrix");

  mat out(x.rows(), x.cols());
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) {
      out(i, j) = cheb(n, x(i, j));
    }
  }
  return out;
}

} // namespace itpp
