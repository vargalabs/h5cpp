/*!
 * \file
 * \brief Implementation of window functions
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger, Adam Piatyszek
 *         and Kumar Appaiah
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
#include <itpp/signal/window.h>
#include <itpp/signal/poly.h>
#include <itpp/base/specmat.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/trig_hyp.h>
#include <itpp/signal/transforms.h>
#include <itpp/base/operators.h>

namespace itpp
{


vec hamming(int n)
{
  vec t(n);

  if (n == 1)
    t(0) = 0.08;
  else
    for (int i = 0;i < n;i++)
      t[i] = (0.54 - 0.46 * std::cos(2.0 * pi * i / (n - 1)));

  return t;
}

vec hanning(int n)
{
  vec t(n);

  for (int i = 0;i < n;i++)
    t(i) = 0.5 * (1.0 - std::cos(2.0 * pi * (i + 1) / (n + 1)));

  return t;
}

// matlab version
vec hann(int n)
{
  vec t(n);

  for (int i = 0;i < n;i++)
    t(i) = 0.5 * (1.0 - std::cos(2.0 * pi * i / (n - 1)));

  return t;
}

vec blackman(int n)
{
  vec t(n);

  for (int i = 0;i < n;i++)
    t(i) = 0.42 - 0.5 * std::cos(2.0 * pi * i / (n - 1)) + 0.08 * std::cos(4.0 * pi * i / (n - 1));

  return t;
}

vec triang(int n)
{
  vec t(n);

  if (n % 2) { // Odd
    for (int i = 0; i < n / 2; i++)
      t(i) = t(n - i - 1) = 2.0 * (i + 1) / (n + 1);
    t(n / 2) = 1.0;
  }
  else
    for (int i = 0; i < n / 2; i++)
      t(i) = t(n - i - 1) = (2.0 * i + 1) / n;

  return t;
}

vec sqrt_win(int n)
{
  vec t(n);

  if (n % 2) { // Odd
    for (int i = 0; i < n / 2; i++)
      t(i) = t(n - i - 1) = std::sqrt(2.0 * (i + 1) / (n + 1));
    t(n / 2) = 1.0;
  }
  else
    for (int i = 0; i < n / 2; i++)
      t(i) = t(n - i - 1) = std::sqrt((2.0 * i + 1) / n);

  return t;
}

vec chebwin(int n, double at)
{
  it_assert((n > 0), "chebwin(): need a positive order n!");

  if (n == 1) {
    return vec("1");
  }

  at = at < 0 ? -at : at;
  // compute the parameter beta
  double beta = std::cosh(::acosh(pow10(at / 20)) / (n - 1));
  vec k = (pi / n) * linspace(0, n - 1, n);
  vec cos_k = cos(k);
  // find the window's DFT coefficients
  vec p = cheb(n - 1, beta * cos_k);

  vec w(n); // the window vector
  // Appropriate IDFT and filling up depending on even/odd n
  if (is_even(n)) {
    w = ifft_real(to_cvec(elem_mult(p, cos_k), elem_mult(p, -sin(k))));
    int half_length = n / 2 + 1;
    w = w.left(half_length) / w(1);
    w = concat(reverse(w), w.right(n - half_length));
  }
  else {
    w = ifft_real(to_cvec(p));
    int half_length = (n + 1) / 2;
    w = w.left(half_length) / w(0);
    w = concat(reverse(w), w.right(n - half_length));
  }
  return w;
}


} // namespace itpp
