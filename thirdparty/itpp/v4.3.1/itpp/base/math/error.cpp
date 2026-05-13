/*!
 * \file
 * \brief Error functions - source file
 * \author Tony Ottosson, Pal Frenger and Adam Piatyszek
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

#include <itpp/base/math/error.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/itcompat.h>


namespace itpp
{

/*!
 * Abramowitz and Stegun: Eq. (7.1.14) gives this continued fraction
 * for erfc(z)
 *
 * erfc(z) = sqrt(pi).exp(-z^2).  1   1/2   1   3/2   2   5/2
 *                               ---  ---  ---  ---  ---  --- ...
 *                               z +  z +  z +  z +  z +  z +
 *
 * This is evaluated using Lentz's method, as described in the
 * narative of Numerical Recipes in C.
 *
 * The continued fraction is true providing real(z) > 0. In practice
 * we like real(z) to be significantly greater than 0, say greater
 * than 0.5.
 */
std::complex<double> cerfc_continued_fraction(const std::complex<double>& z)
{
  const double tiny = std::numeric_limits<double>::min();

  // first calculate z+ 1/2   1
  //                    ---  --- ...
  //                    z +  z +
  std::complex<double> f(z);
  std::complex<double> C(f);
  std::complex<double> D(0.0);
  std::complex<double> delta;
  double a;

  a = 0.0;
  do {
    a += 0.5;
    D = z + a * D;
    C = z + a / C;
    if ((D.real() == 0.0) && (D.imag() == 0.0))
      D = tiny;
    D = 1.0 / D;
    delta = C * D;
    f = f * delta;
  }
  while (abs(1.0 - delta) > eps);

  // Do the first term of the continued fraction
  f = 1.0 / f;

  // and do the final scaling
  f = f * exp(-z * z) / std::sqrt(pi);

  return f;
}

//! Complementary function to \c cerfc_continued_fraction
std::complex<double> cerf_continued_fraction(const std::complex<double>& z)
{
  if (z.real() > 0)
    return 1.0 - cerfc_continued_fraction(z);
  else
    return -1.0 + cerfc_continued_fraction(-z);
}

/*!
 * Abramawitz and Stegun: Eq. (7.1.5) gives a series for erf(z) good
 * for all z, but converges faster for smallish abs(z), say abs(z) < 2.
 */
std::complex<double> cerf_series(const std::complex<double>& z)
{
  const double tiny = std::numeric_limits<double>::min();
  std::complex<double> sum(0.0);
  std::complex<double> term(z);
  std::complex<double> z2(z*z);

  for (int n = 0; (n < 3) || (abs(term) > abs(sum) * tiny); n++) {
    sum += term / static_cast<double>(2 * n + 1);
    term *= -z2 / static_cast<double>(n + 1);
  }

  return sum * 2.0 / std::sqrt(pi);
}

/*!
 * Numerical Recipes quotes a formula due to Rybicki for evaluating
 * Dawson's Integral:
 *
 * exp(-x^2) integral exp(t^2).dt = 1/sqrt(pi) lim  sum  exp(-(z-n.h)^2) / n
 *            0 to x                           h->0 n odd
 *
 * This can be adapted to erf(z).
 */
std::complex<double> cerf_rybicki(const std::complex<double>& z)
{
  double h = 0.2; // numerical experiment suggests this is small enough

  // choose an even n0, and then shift z->z-n0.h and n->n-h.
  // n0 is chosen so that real((z-n0.h)^2) is as small as possible.
  int n0 = 2 * static_cast<int>(z.imag() / (2 * h) + 0.5);

  std::complex<double> z0(0.0, n0 * h);
  std::complex<double> zp(z - z0);
  std::complex<double> sum(0.0, 0.0);

  // limits of sum chosen so that the end sums of the sum are
  // fairly small. In this case exp(-(35.h)^2)=5e-22
  for (int np = -35; np <= 35; np += 2) {
    std::complex<double> t(zp.real(), zp.imag() - np * h);
    std::complex<double> b(exp(t * t) / static_cast<double>(np + n0));
    sum += b;
  }

  sum *= 2.0 * exp(-z * z) / pi;

  return std::complex<double>(-sum.imag(), sum.real());
}

/*
 * This function calculates a well known error function erf(z) for
 * complex z. Three methods are implemented. Which one is used
 * depends on z.
 */
std::complex<double> erf(const std::complex<double>& z)
{
  // Use the method appropriate to size of z -
  // there probably ought to be an extra option for NaN z, or infinite z
  if (abs(z) < 2.0)
    return cerf_series(z);
  else {
    if (std::abs(z.real()) < 0.5)
      return cerf_rybicki(z);
    else
      return cerf_continued_fraction(z);
  }
}


double erfinv(double P)
{
  double Y, A, B, X, Z, W, WI, SN, SD, F, Z2, SIGMA;
  double A1 = -.5751703, A2 = -1.896513, A3 = -.5496261E-1;
  double B0 = -.1137730, B1 = -3.293474, B2 = -2.374996, B3 = -1.187515;
  double C0 = -.1146666, C1 = -.1314774, C2 = -.2368201, C3 = .5073975e-1;
  double D0 = -44.27977, D1 = 21.98546, D2 = -7.586103;
  double E0 = -.5668422E-1, E1 = .3937021, E2 = -.3166501, E3 = .6208963E-1;
  double F0 = -6.266786, F1 = 4.666263, F2 = -2.962883;
  double G0 = .1851159E-3, G1 = -.2028152E-2, G2 = -.1498384, G3 = .1078639E-1;
  double H0 = .9952975E-1, H1 = .5211733, H2 = -.6888301E-1;
  // double RINFM=1.7014E+38;

  X = P;
  SIGMA = sign(X);
  it_error_if(X < -1 || X > 1, "erfinv : argument out of bounds");
  Z = fabs(X);
  if (Z > .85) {
    A = 1 - Z;
    B = Z;
    W = std::sqrt(-log(A + A * B));
    if (W >= 2.5) {
      if (W >= 4.) {
        WI = 1. / W;
        SN = ((G3 * WI + G2) * WI + G1) * WI;
        SD = ((WI + H2) * WI + H1) * WI + H0;
        F = W + W * (G0 + SN / SD);
      }
      else {
        SN = ((E3 * W + E2) * W + E1) * W;
        SD = ((W + F2) * W + F1) * W + F0;
        F = W + W * (E0 + SN / SD);
      }
    }
    else {
      SN = ((C3 * W + C2) * W + C1) * W;
      SD = ((W + D2) * W + D1) * W + D0;
      F = W + W * (C0 + SN / SD);
    }
  }
  else {
    Z2 = Z * Z;
    F = Z + Z * (B0 + A1 * Z2 / (B1 + Z2 + A2 / (B2 + Z2 + A3 / (B3 + Z2))));
  }
  Y = SIGMA * F;
  return Y;
}

double Qfunc(double x)
{
  return (0.5 * ::erfc(x / 1.41421356237310));
}


// Error function
vec erf(const vec &x) { return apply_function<double>(::erf, x); }
mat erf(const mat &x) { return apply_function<double>(::erf, x); }
cvec erf(const cvec &x)
{
  return apply_function<std::complex<double> >(erf, x);
}
cmat erf(const cmat &x)
{
  return apply_function<std::complex<double> >(erf, x);
}

// Inverse of error function
vec erfinv(const vec &x) { return apply_function<double>(erfinv, x); }
mat erfinv(const mat &x) { return apply_function<double>(erfinv, x); }

// Complementary error function
vec erfc(const vec &x) { return apply_function<double>(::erfc, x); }
mat erfc(const mat &x) { return apply_function<double>(::erfc, x); }

// Q-function
vec Qfunc(const vec &x) { return apply_function<double>(Qfunc, x); }
mat Qfunc(const mat &x) { return apply_function<double>(Qfunc, x); }

} // namespace itpp
