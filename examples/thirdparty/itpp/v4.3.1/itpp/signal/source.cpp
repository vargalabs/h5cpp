/*!
 * \file
 * \brief Deterministic sources - source file
 * \author Tobias Ringstrom and Tony Ottosson
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

#include <itpp/signal/source.h>


namespace itpp
{

///////////////////////////////////////////////
// Sine_Source
///////////////////////////////////////////////

Sine_Source::Sine_Source(double freq, double mean, double ampl, double inphase)
{
  A = ampl;
  m = mean;
  theta = inphase;
  dtheta = 2.0 * pi * freq;
}

double Sine_Source::sample()
{
  double samp = m + A * sin(theta);

  theta += dtheta;
  if (theta >= 2.0 * pi)
    theta -= 2.0 * pi;

  return samp;
}

vec Sine_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Sine_Source::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Square_Source
///////////////////////////////////////////////

Square_Source::Square_Source(double freq, double mean, double ampl, double inphase)
{
  A = ampl;
  m = mean;
  theta = inphase / (2.0 * pi);
  dtheta = freq;
}

double Square_Source::sample()
{
  double samp = theta < 0.5 ? 1.0 : -1.0;

  theta += dtheta;
  if (theta >= 1.0)
    theta -= 1.0;

  return samp;
}

vec Square_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Square_Source::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Triangle_Source
///////////////////////////////////////////////

Triangle_Source::Triangle_Source(double freq, double mean, double ampl, double inphase)
{
  A = ampl;
  m = mean;
  theta = inphase / (2.0 * pi);
  dtheta = freq;
}

double Triangle_Source::sample()
{
  double samp = m + 4.0 * A * (theta < 0.25 ? theta : 0.5 - theta);

  theta += dtheta;
  if (theta >= 0.75)
    theta -= 1.0;

  return samp;
}

vec Triangle_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Triangle_Source::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Sawtooth_Source
///////////////////////////////////////////////

Sawtooth_Source::Sawtooth_Source(double freq, double mean, double ampl, double inphase)
{
  A = ampl;
  m = mean;
  theta = inphase / (2.0 * pi);
  dtheta = freq;
}

double Sawtooth_Source::sample()
{
  double samp = 2.0 * A * theta;

  theta += dtheta;
  if (theta >= 0.5)
    theta -= 1.0;

  return samp;
}

vec Sawtooth_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Sawtooth_Source::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Impulse_Source
///////////////////////////////////////////////

Impulse_Source::Impulse_Source(double freq, double ampl, double inphase)
{
  A = ampl;
  pos = inphase / (2.0 * pi);
  dtheta = freq;
}

double Impulse_Source::sample()
{
  double samp;

  if (pos >= 1.0) {
    samp = A;
    pos -= 1.0;
  }
  else {
    samp = 0.0;
    pos += dtheta;
  }

  return samp;
}

vec Impulse_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Impulse_Source::operator()(int h, int w)
{
  mat m(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      m(i, j) = sample();

  return m;
}

///////////////////////////////////////////////
// Pattern_Source
///////////////////////////////////////////////

Pattern_Source::Pattern_Source(const vec &pattern, int start_pos)
{
  pat = pattern;
  pos = start_pos;

  // Calculate the mean and variance.  Note that the variance shall
  // be normalied by N and not N-1 in this case
  mean = var = 0.0;
  for (int i = pat.size() - 1; i >= 0; i--) {
    mean += pat(i);
    var += pat(i) * pat(i);
  }
  mean /= pat.size();
  var /= pat.size();
  var -= mean * mean;
}

double Pattern_Source::sample()
{
  double samp = pat(pos);

  if (pos >= pat.size() - 1)
    pos = 0;
  else
    pos++;

  return samp;
}

vec Pattern_Source::operator()(int n)
{
  vec v(n);

  for (int i=0; i < n; i++)
    v(i) = sample();

  return v;
}

mat Pattern_Source::operator()(int h, int w)
{
  mat m(h, w);
  int i, j;

  for (i = 0; i < h; i++)
    for (j = 0; j < w; j++)
      m(i, j) = sample();

  return m;
}

} // namespace itpp
