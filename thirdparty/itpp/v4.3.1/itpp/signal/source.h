/*!
 * \file
 * \brief Deterministic sources - header file
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

#ifndef SOURCE_H
#define SOURCE_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>

namespace itpp
{

//! \addtogroup detsource

/*!
  \brief Sine-wave source
  \ingroup detsource
*/
class ITPP_EXPORT Sine_Source
{
public:
  //! Constructor. Set frequency, mean, amplitude, and start phase
  Sine_Source(double freq, double mean = 0.0, double ampl = 1.0, double inphase = 0.0);
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  double m, A, theta, dtheta;
};

/*!
  \brief Square-wave source
  \ingroup detsource
*/
class ITPP_EXPORT Square_Source
{
public:
  //! Constructor. Set frequency, mean, amplitude, and start phase
  Square_Source(double freq, double mean = 0.0, double ampl = 1.0, double inphase = 0.0);
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  double m, A, theta, dtheta;
};

/*!
  \brief Triangle-wave source
  \ingroup detsource
*/
class ITPP_EXPORT Triangle_Source
{
public:
  //! Constructor. Set frequency, mean, amplitude and start phase
  Triangle_Source(double freq, double mean = 0.0, double ampl = 1.0, double inphase = 0.0);
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  double m, A, theta, dtheta;
};

/*!
  \brief Sawtooth-wave source
  \ingroup detsource
*/
class ITPP_EXPORT Sawtooth_Source
{
public:
  //! Constructor. Set frequency, mean, amplitude, and start phase
  Sawtooth_Source(double freq, double mean = 0.0, double ampl = 1.0, double inphase = 0.0);
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  double m, A, theta, dtheta;
};

/*!
  \brief Impulse source
  \ingroup detsource
*/
class ITPP_EXPORT Impulse_Source
{
public:
  //! Constructor. Set frequency, amplitude and start phase
  Impulse_Source(double freq, double ampl = 1.0, double inphase = 0.0);
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  double A, pos, dtheta;
};

/*!
  \brief Pattern source
  \ingroup detsource
*/
class ITPP_EXPORT Pattern_Source
{
public:
  //! Constructor. Set pattern and start position
  Pattern_Source(const vec &pattern, int start_pos = 0);
  //! Destructor
  virtual ~Pattern_Source() { }
  //! Return a single sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
protected:
private:
  double sample();
  int pos;
  vec pat;
  double mean, var;
};

} //namespace itpp

#endif // #ifndef SOURCE_H
