/*!
 * \file
 * \brief Pulse shaping classes - header file
 * \author Tony Ottosson, Hakan Eriksson and Adam Piatyszek
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

#ifndef PULSE_SHAPE_H
#define PULSE_SHAPE_H

#include <itpp/base/vec.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/math/trig_hyp.h>
#include <itpp/signal/filter.h>
#include <itpp/signal/resampling.h>


namespace itpp
{

/*!
  \brief General FIR Pulse Shape

  Upsamples and shapes symbols according to a given FIR filter.
  Observe that since the shaping is done with a FIR filter, the
  first samples in the output are zero or small before the memory of
  the filter is filled.

  The class is templated as follows:
  <ul>
  <li> \c T1 is the type of the input samples</li>
  <li> \c T2 is the type of the filter coefficients</li>
  <li> \c T3 is the type of the output samples</li>
  </ul>

  An example of usage is:

  \code
  #include "itpp/itcomm.h"

  vec filter_response;
  filter_response ="0.7 0.3 0.6";
  Pulse_Shape<double,double,double> shaper(filter_response, 4);
  BPSK bpsk;
  vec symbols, samples;

  symbols = bpsk.modulate_bits(randb(20));
  samples = shaper.shape_symbols(symbols);
  \endcode
*/
template<class T1, class T2, class T3>
class Pulse_Shape
{
public:
  //! Constructor
  Pulse_Shape();
  //! Constructor
  Pulse_Shape(const Vec<T2> &impulse_response, int upsampling_factor);
  //! Destructor
  virtual ~Pulse_Shape() {}
  /*!
    \brief Set the general impulse response of the FIR filter

    Observe that the pulse shape must have a duration of an integer
    number of symbols. Thus the length of the impulse response-1
    modulo over sampling is an integer.
  */
  void set_pulse_shape(const Vec<T2> &impulse_response, int upsampling_factor);
  //! Get the pulse shape
  Vec<T2> get_pulse_shape(void) const;
  //! Get the over sampling factor
  int get_upsampling_factor() const;
  //! Get the length of the pulse in number of symbols
  int get_pulse_length() const;
  //! Get the length of the internal FIR filter
  int get_filter_length() const;

  //! Shape the input symbols performing upsampling
  void shape_symbols(const Vec<T1> &input, Vec<T3> &output);
  //! Shape the input symbols performing upsampling
  Vec<T3> shape_symbols(const Vec<T1> &input);

  //! Shape the input samples already upsampled
  void shape_samples(const Vec<T1> &input, Vec<T3> &output);
  //! Shape the input symbols already upsampled
  Vec<T3> shape_samples(const Vec<T1> &input);

  //! Clear internal states
  void clear(void);

protected:
  //! The impulse resounse of the pulse shaping filter
  Vec<T2> impulse_response;
  //! The pulse shaping filter
  MA_Filter<T1, T2, T3> shaping_filter;
  //! Length in symbols
  int pulse_length;
  //! Samples per input symbol
  int upsampling_factor;
  //! Ensures that setup is called before any other member function
  bool setup_done;
};

/*!
  \brief Raised Cosine (RC) Pulse Shaper

  Upsamples and shapes symbols as raised cosine pulses with a given
  roll-off factor \f$ \alpha \f$.
  The raised cosine pulse shape is defined as:
  \f[
  p(t) = \frac{\sin(\pi t / T)}{\pi t / T}
  \frac{\cos(\alpha \pi t / T)}{1 - (2 \alpha t / T)^2}
  \f]
  For more details see e.g. Lee & Messerschmitt, p. 190. Observe
  that the shaping is done with a FIR filter where the size is given
  by \a filter_length * \a over_sample_factor + 1.  The first
  samples in the output will therefore be zero or small before the
  memory of the filter is filled.

  What is important, when using RC shaping in a transmission system
  with the AWGN channel, the mean power of the output samples is not
  normalised, so the channel noise variance (or shaped signal)
  should be scaled appropriately.

  The class is templated as follows: \c T1 is the type of the input
  and the output samples.
  An example of usage is:

  \code
  #include "itpp/itcomm.h"

  Raised_Cosine<double> rc(0.5, 6, 8);
  BPSK bpsk;
  vec symbols, samples;

  symbols = bpsk.modulate_bits(randb(20));
  samples = rc.shape_symbols(symbols);
  \endcode
*/
template<class T1>
class Raised_Cosine : public Pulse_Shape<T1, double, T1>
{
public:
  //! Constructor
  Raised_Cosine() {}
  //! Constructor
  Raised_Cosine(double roll_off, int filter_length = 6, int upsampling_factor = 8);
  //! Destructor
  virtual ~Raised_Cosine() {}
  //! Set pulse shape (roll_off_factor between 0 and 1, filter_length even)
  void set_pulse_shape(double roll_off_factor, int filter_length = 6, int upsampling_factor = 8);
  //! Get the roll-off factor
  double get_roll_off(void) const;

protected:
  //! The roll off factor (i.e. alpha)
  double roll_off_factor;
};

/*!
  \brief (Square) Root Raised Cosine (RRC) Pulse Shaper

  Upsamples and shapes symbols as square root raised cosine pulses
  with a given roll-off factor \f$ \alpha \f$ (zero is not allowed
  - use raised cosine instead).
  The Root Raised Cosine pulse shape is defined as:
  \f[
  p(t) = \frac{4 \alpha}{\pi \sqrt{T}} \frac{\cos((1+\alpha)\pi t / T)
  + T \sin((1-\alpha)\pi t / T) / (4 \alpha t) }{1 - (4 \pi t / T)^2}
  \f]

  It is called square root raised cosine since it is defined as the
  square root of the raised cosine pulse in the frequency domain.
  Thus with a transmitter pulse shape of root raised cosine and a
  receiver filter that is root raised cosine, the overall response
  will be a raised cosine.

  For more details see e.g. Lee & Messerschmitt, p. 228. Observe
  that the shaping is done with a FIR filter where the size is given
  by \a filter_length * \a over_sample_factor + 1. The first samples
  in the output will therefore be zero or small before the memory of
  the filter is filled.

  What is important, when using RRC shaping in a transmission system
  with the AWGN channel, the mean power of the output samples is not
  normalised, so the channel noise variance (or shaped signal)
  should be scaled appropriately.

  The class is templated as follows: \c T1 is the type of the input
  and the output samples.
  An example of usage is:

  \code
  #include "itpp/itcomm.h"

  Root_Raised_Cosine<double> rrc(0.5,6,8);
  BPSK bpsk;
  vec symbols, samples;

  symbols = bpsk.modulate_bits(randb(20));
  samples = rrc.shape_symbols(symbols);
  \endcode
*/
template<class T1>
class Root_Raised_Cosine : public Pulse_Shape<T1, double, T1>
{
public:
  //! Constructor
  Root_Raised_Cosine() {}
  //! Constructor
  Root_Raised_Cosine(double roll_off_factor, int filter_length = 6, int upsampling_factor = 8);
  //! Destructor
  virtual ~Root_Raised_Cosine() {}
  //! Set pulse_shape, roll_off_factor between 0 and 1, filter_length even
  void set_pulse_shape(double roll_off_factor, int filter_length = 6, int upsampling_factor = 8);
  //! Get the Roll-off factor
  double get_roll_off(void) const;

protected:
  //! The roll off factor (i.e. alpha)
  double roll_off_factor;
};

//-------------------------------------------------------------------------
// Implementation of templated code starts here
//-------------------------------------------------------------------------

//---------------------------- Pulse_Shape --------------------------------

template<class T1, class T2, class T3>
Pulse_Shape<T1, T2, T3>::Pulse_Shape()
{
  setup_done = false;
  pulse_length = 0;
  upsampling_factor = 0;
}


template<class T1, class T2, class T3>
Pulse_Shape<T1, T2, T3>::Pulse_Shape(const Vec<T2> &impulse_response, int upsampling_factor)
{
  set_pulse_shape(impulse_response, upsampling_factor);
}

template<class T1, class T2, class T3>
void Pulse_Shape<T1, T2, T3>::set_pulse_shape(const Vec<T2> &impulse_response_in, int upsampling_factor_in)
{
  it_error_if(impulse_response_in.size() == 0, "Pulse_Shape: impulse response is zero length");
  it_error_if(upsampling_factor_in < 1, "Pulse_Shape: incorrect upsampling factor");

  pulse_length = (impulse_response_in.size() - 1) / upsampling_factor_in;
  upsampling_factor = upsampling_factor_in;

  impulse_response = impulse_response_in;
  shaping_filter.set_coeffs(impulse_response);
  shaping_filter.clear();
  setup_done = true;
}

template<class T1, class T2, class T3>
Vec<T2> Pulse_Shape<T1, T2, T3>::get_pulse_shape(void) const
{
  return impulse_response;
}

template<class T1, class T2, class T3>
int Pulse_Shape<T1, T2, T3>::get_upsampling_factor(void) const
{
  return upsampling_factor;
}

template<class T1, class T2, class T3>
int Pulse_Shape<T1, T2, T3>::get_pulse_length(void) const
{
  return pulse_length;
}

template<class T1, class T2, class T3>
int Pulse_Shape<T1, T2, T3>::get_filter_length(void) const
{
  return impulse_response.size();
}

template<class T1, class T2, class T3>
void Pulse_Shape<T1, T2, T3>::shape_symbols(const Vec<T1>& input, Vec<T3> &output)
{
  it_assert(setup_done, "Pulse_Shape must be set up before using");
  it_error_if(pulse_length == 0, "Pulse_Shape: impulse response is zero length");
  it_error_if(input.size() == 0, "Pulse_Shape: input is zero length");

  if (upsampling_factor > 1)
    output = shaping_filter(upsample(input, upsampling_factor));
  else
    output = input;
}

template<class T1, class T2, class T3>
Vec<T3> Pulse_Shape<T1, T2, T3>::shape_symbols(const Vec<T1>& input)
{
  it_assert(setup_done, "Pulse_Shape must be set up before using");
  Vec<T3> temp;
  shape_symbols(input, temp);
  return temp;
}

template<class T1, class T2, class T3>
void Pulse_Shape<T1, T2, T3>::shape_samples(const Vec<T1>& input, Vec<T3> &output)
{
  it_assert(setup_done, "Pulse_Shape must be set up before using");
  it_error_if(pulse_length == 0, "Pulse_Shape: impulse response is zero length");
  it_error_if(input.size() == 0, "Pulse_Shape: input is zero length");

  if (upsampling_factor > 1)
    output = shaping_filter(input);
  else
    output = input;
}

template<class T1, class T2, class T3>
Vec<T3> Pulse_Shape<T1, T2, T3>::shape_samples(const Vec<T1>& input)
{
  it_assert(setup_done, "Pulse_Shape must be set up before using");
  Vec<T3> temp;
  shape_samples(input, temp);
  return temp;
}

template<class T1, class T2, class T3>
void Pulse_Shape<T1, T2, T3>::clear(void)
{
  it_assert(setup_done, "Pulse_Shape must be set up before using");
  shaping_filter.clear();
}

//-------------------- Raised_Cosine -----------------------------------

template<class T1>
Raised_Cosine<T1>::Raised_Cosine(double roll_off_factor, int filter_length, int upsampling_factor)
{
  set_pulse_shape(roll_off_factor, filter_length, upsampling_factor);
}

template<class T1>
void Raised_Cosine<T1>::set_pulse_shape(double roll_off_factor_in, int filter_length, int upsampling_factor_in)
{
  it_error_if(roll_off_factor_in < 0 || roll_off_factor_in > 1, "Raised_Cosine: roll-off out of range");
  roll_off_factor = roll_off_factor_in;

  it_assert(is_even(filter_length), "Raised_Cosine: Filter length not even");

  int i;
  double t, den;
  this->upsampling_factor = upsampling_factor_in;
  this->pulse_length = filter_length;
  this->impulse_response.set_size(filter_length * upsampling_factor_in + 1,
                                  false);

  for (i = 0; i < this->impulse_response.size(); i++) {
    // delayed to be casual
    t = (double)(i - filter_length * upsampling_factor_in / 2)
        / upsampling_factor_in;
    den = 1 - sqr(2 * roll_off_factor * t);
    if (den == 0) {
      // exception according to "The Care and feeding of digital,
      // pulse-shaping filters" by Ken Gentile,
      // the limit of raised cosine impulse responce function,
      // as (alpha * t / tau) approaches (+- 0.5) is given as:
      this->impulse_response(i) = sinc(t) * pi / 4;
    }
    else {
      this->impulse_response(i) = std::cos(roll_off_factor * pi * t)
                                  * sinc(t) / den;
    }
  }

  // BUGFIX: Commented out to achieve similar results to Matlab
  // rcosfil function. Now the concatenation of two root-raised
  // cosine filters gives tha same results as a one raised cosine
  // shaping function.
  // this->impulse_response /= std::sqrt(double(this->upsampling_factor));
  this->shaping_filter.set_coeffs(this->impulse_response);
  this->shaping_filter.clear();
  this->setup_done = true;
}

template<class T1>
double Raised_Cosine<T1>::get_roll_off(void) const
{
  it_assert(this->setup_done, "Pulse_Shape must be set up before using");
  return roll_off_factor;
}

//-------------------- Root_Raised_Cosine -----------------------------------

template<class T1>
Root_Raised_Cosine<T1>::Root_Raised_Cosine(double roll_off_factor, int filter_length, int upsampling_factor)
{
  set_pulse_shape(roll_off_factor, filter_length, upsampling_factor);
}

template<class T1>
void Root_Raised_Cosine<T1>::set_pulse_shape(double roll_off_factor_in, int filter_length, int upsampling_factor_in)
{
  it_error_if(roll_off_factor_in <= 0 || roll_off_factor_in > 1,
              "Root_Raised_Cosine: roll-off out of range");
  roll_off_factor = roll_off_factor_in;

  it_assert(is_even(filter_length),
            "Root_Raised_Cosine: Filter length not even");

  int i;
  double t, num, den, tmp_arg;
  this->upsampling_factor = upsampling_factor_in;
  this->pulse_length = filter_length;
  this->impulse_response.set_size(filter_length * upsampling_factor_in + 1,
                                  false);

  for (i = 0; i < this->impulse_response.size(); i++) {
    // delayed to be casual
    t = (double)(i - filter_length * upsampling_factor_in / 2)
        / upsampling_factor_in;
    den = 1 - sqr(4 * roll_off_factor * t);
    if (t == 0) {
      this->impulse_response(i) = 1 + (4 * roll_off_factor / pi)
                                  - roll_off_factor;
    }
    else if (den == 0) {
      tmp_arg = pi / (4 * roll_off_factor);
      this->impulse_response(i) = roll_off_factor / std::sqrt(2.0)
                                  * ((1 + 2 / pi) * std::sin(tmp_arg) + (1 - 2 / pi) * std::cos(tmp_arg));
    }
    else {
      num = std::sin(pi * (1 - roll_off_factor) * t)
            + std::cos(pi * (1 + roll_off_factor) * t) * 4 * roll_off_factor * t;
      this->impulse_response(i) = num / (pi * t * den);
    }
  }

  this->impulse_response /= std::sqrt(double(upsampling_factor_in));
  this->shaping_filter.set_coeffs(this->impulse_response);
  this->shaping_filter.clear();
  this->setup_done = true;
}

template<class T1>
double Root_Raised_Cosine<T1>::get_roll_off(void) const
{
  it_assert(this->setup_done, "Pulse_Shape must be set up before using");
  return roll_off_factor;
}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Pulse_Shape<double, double, double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Pulse_Shape < std::complex<double>, double,
  std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Pulse_Shape < std::complex<double>, std::complex<double>,
  std::complex<double> >;

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Root_Raised_Cosine<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Root_Raised_Cosine<std::complex<double> >;

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Raised_Cosine<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Raised_Cosine<std::complex<double> >;

//! \endcond

} // namespace itpp

#endif // #ifndef PULSE_SHAPE_H
