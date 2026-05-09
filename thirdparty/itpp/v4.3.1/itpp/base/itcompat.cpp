/*!
 * \file
 * \brief IT++ compatibility types and functions
 * \author Adam Piatyszek
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
#include <itpp/base/itassert.h>
#include <limits>

//! \cond

#ifndef HAVE_TGAMMA
// "True" gamma function
double tgamma(double x)
{
  double s = (2.50662827510730 + 190.9551718930764 / (x + 1)
              - 216.8366818437280 / (x + 2) + 60.19441764023333
              / (x + 3) - 3.08751323928546 / (x + 4) + 0.00302963870525
              / (x + 5) - 0.00001352385959072596 / (x + 6)) / x;
  if (s < 0)
    return -exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(-s));
  else
    return exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(s));
}
#endif

#if !defined(HAVE_LGAMMA) || (HAVE_DECL_SIGNGAM != 1)
// The sign of the Gamma function is returned in the external integer
// signgam declared in <math.h>. It is 1 when the Gamma function is positive
// or zero, -1 when it is negative. However, MinGW definition of lgamma()
// function does not use the global signgam variable.
int signgam;
// Logarithm of an absolute value of gamma function
double lgamma(double x)
{
  double gam = tgamma(x);
  signgam = (gam < 0) ? -1 : 1;
  return log(fabs(gam));
}
#endif

#ifndef HAVE_CBRT
// Cubic root
double cbrt(double x) { return std::pow(x, 1.0 / 3.0); }
#endif


#ifndef HAVE_EXPM1
#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212146 // ln(2)
#endif
// Implementation taken from GSL (http://www.gnu.org/software/gsl/)
double expm1(double x)
{
  /* FIXME: this should be improved */
  if (std::fabs(x) < M_LN2) {
    /* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */
    double i = 1.0;
    double sum = x;
    double term = x / 1.0;
    do {
      i++;
      term *= x / i;
      sum += term;
    }
    while (std::fabs(term) > (std::fabs(sum)
                              * std::numeric_limits<double>::epsilon()));
    return sum;
  }
  else {
    return std::exp(x) - 1.0;
  }
}
#endif // HAVE_EXPM1


#ifndef HAVE_ERFC
// Complementary error function
double erfc(double Y)
{
  int  ISW, I;
  double P[4], Q[3], P1[6], Q1[5], P2[4], Q2[3];
  double XMIN, XLARGE, SQRPI;
  double X, RES, XSQ, XNUM, XDEN, XI, XBIG, ERFCret;
  P[1] = 0.3166529;
  P[2] = 1.722276;
  P[3] = 21.38533;
  Q[1] = 7.843746;
  Q[2] = 18.95226;
  P1[1] = 0.5631696;
  P1[2] = 3.031799;
  P1[3] = 6.865018;
  P1[4] = 7.373888;
  P1[5] = 4.318779e-5;
  Q1[1] = 5.354217;
  Q1[2] = 12.79553;
  Q1[3] = 15.18491;
  Q1[4] = 7.373961;
  P2[1] = 5.168823e-2;
  P2[2] = 0.1960690;
  P2[3] = 4.257996e-2;
  Q2[1] = 0.9214524;
  Q2[2] = 0.1509421;
  XMIN = 1.0E-5;
  XLARGE = 4.1875E0;
  XBIG = 9.0;
  SQRPI = 0.5641896;
  X = Y;
  ISW = 1;
  if (X < 0) {
    ISW = -1;
    X = -X;
  }
  if (X < 0.477) {
    if (X >= XMIN) {
      XSQ = X * X;
      XNUM = (P[1] * XSQ + P[2]) * XSQ + P[3];
      XDEN = (XSQ + Q[1]) * XSQ + Q[2];
      RES = X * XNUM / XDEN;
    }
    else RES = X * P[3] / Q[2];
    if (ISW == -1) RES = -RES;
    RES = 1.0 - RES;
    goto slut;
  }
  if (X > 4.0) {
    if (ISW > 0) goto ulf;
    if (X < XLARGE) goto eva;
    RES = 2.0;
    goto slut;
  }
  XSQ = X * X;
  XNUM = P1[5] * X + P1[1];
  XDEN = X + Q1[1];
  for (I = 2;I <= 4;I++) {
    XNUM = XNUM * X + P1[I];
    XDEN = XDEN * X + Q1[I];
  }
  RES = XNUM / XDEN;
  goto elin;
ulf:
  if (X > XBIG) goto fred;
eva:
  XSQ = X * X;
  XI = 1.0 / XSQ;
  XNUM = (P2[1] * XI + P2[2]) * XI + P2[3];
  XDEN = XI + Q2[1] * XI + Q2[2];
  RES = (SQRPI + XI * XNUM / XDEN) / X;
elin:
  RES = RES * exp(-XSQ);
  if (ISW == -1) RES = 2.0 - RES;
  goto slut;
fred:
  RES = 0.0;
slut:
  ERFCret = RES;
  return ERFCret;
}
#endif


#ifndef HAVE_ASINH
// Arcus sinhyp
double asinh(double x)
{
  return ((x >= 0) ? std::log(x + std::sqrt(x * x + 1))
          : -std::log(-x + std::sqrt(x * x + 1)));
}
#endif

#ifndef HAVE_ACOSH
// Arcus coshyp
double acosh(double x)
{
  it_error_if(x < 1, "acosh(): Argument out of range");
  return std::log(x + std::sqrt(x * x - 1.0));
}
#endif

#ifndef HAVE_ATANH
// Arcus tanhyp
double atanh(double x)
{
  it_error_if(std::fabs(x) >= 1, "atanh(): Argument out of range");
  return 0.5 * std::log((x + 1) / (x - 1));
}
#endif


#ifndef HAVE_RINT
double rint(double x)
{
  // zero or NaN case
  if ((x == 0.0) || (x != x))
    return x;

  // negative case
  bool neg = false;
  if (x < 0.0) {
    x = -x;
    neg = true;
  }

  double y = std::floor(x + 0.5);
  int i = static_cast<int>(y);
  if ((y - x >= 0.5) && (i & 1))
    --y;

  return neg ? -y : y;
}
#endif // HAVE_RINT

//! \endcond
