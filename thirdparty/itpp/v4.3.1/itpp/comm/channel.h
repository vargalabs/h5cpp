/*!
 * \file
 * \brief Communication channels' classes - header file
 * \author Tony Ottosson, Pal Frenger, Adam Piatyszek and Zbigniew Dlugaszewski
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

#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/mat.h>
#include <itpp/base/array.h>
#include <itpp/base/random.h>
#include <itpp/signal/filter.h>
#include <itpp/itexports.h>

/*!
 * \addtogroup channels
 * \brief Communication Channel Models
 * \author Tony Ottosson, Pal Frenger, Adam Piatyszek and Zbigniew Dlugaszewski
 *
 *
 * \section toc Table of Contents
 *
 * - \ref channel_intro
 * - \ref channel_doppler
 * - \ref channel_freq_sel
 * - \ref channel_los
 * - \ref channel_ref
 *
 *
 * \section channel_intro Introduction
 *
 * When simulating a communication link, a model of the channel behaviour
 * is usually needed. Such a model typically consist of three parts: the
 * propagation attenuation, the shadowing (log-normal fading) and the
 * multipath fading (small scale fading). In the following we will focus
 * on the small scale (or multipath) fading.
 *
 * Multipath fading is the process where the received signal is a sum of
 * many reflections, each with different propagation time, phase and
 * attenuation. The sum signal will vary in time if the receiver (or
 * transmitter) moves, or if some of the reflectors move. We usually refer
 * to this process as a fading process and try to model it as a stochastic
 * process. The most common model is the Rayleigh fading model where the
 * process is modeled as a sum of infinitely many (in practise it is
 * enough with only a few) received reflections from all angles (uniformly
 * distributed) around the receiver. Mathematically we write the received
 * signal, \f$ r(t) \f$ as
 *
 * \f[ r(t) = a(t) * s(t), \f]
 *
 * where \f$ s(t) \f$ is the transmitted signal and \f$ a(t) \f$ is the
 * complex channel coefficient (or fading process). If this process is
 * modeled as a Rayleigh fading process then \f$ a(t) \f$ is a complex
 * Gaussian process and the envelope \f$ \|a(t)\| \f$ is Rayleigh
 * distributed.
 *
 *
 * \section channel_doppler Doppler
 *
 * The speed by which the channel changes is decided by the speed of the
 * mobile (transmitter or receiver or both). This movement will cause the
 * channel coefficient \f$ a(t) \f$ to be correlated in time (or
 * equivalently in frequency). Different models of this correlation exist,
 * but the most common is the classical Jakes model where the correlation
 * function is given as
 *
 * \f[ R(\tau) = E[a^*(t) a(t+\tau)] = J_0(2 \pi f_\mathrm{max} \tau), \f]
 *
 * where \f$ f_\mathrm{max} \f$ is the maximum Doppler frequency given by
 *
 * \f[ f_\mathrm{max} = \frac{v}{\lambda} = \frac{v}{c_0} f_c. \f]
 *
 * Here \f$ c_0 \f$ is the speed of light and \f$ f_c \f$ is the carrier
 * frequency. Often the maximum Doppler frequency is given as the
 * normalized Doppler \f$ f_\mathrm{max} T_s \f$ , where \f$ T_s \f$ is
 * the sample duration (often the symbol time) of the simulated system.
 * Instead of specifying the correlation function \f$ R(\tau) \f$ one can
 * specify the Doppler spectrum (the Fourier transform of \f$ R(\tau)
 * \f$).
 *
 *
 * \section channel_freq_sel Frequency-selective Channels
 *
 * Since \f$ a(t) \f$ affects the transmitted signal as a constant scaling
 * factor at a given time, this channel model is often referred to as
 * flat-fading (or frequency non-selective fading). On the other hand, if
 * time arrivals of the reflections are very different (compared to the
 * sampling time), we cannot model the received signal only as a scaled
 * version of the transmitted signal. Instead we model the channel as
 * frequency-selective but time-invariant (or at least wide-sense
 * stationary) with the impulse response
 *
 * \f[ h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k )
 * \delta(t-\tau_k), \f]
 *
 * where \f$ N_\mathrm{taps} \f$ is the number of channel taps, \f$ a_k
 * \f$ is the average amplitude at delay \f$ \tau_k \f$, and \f$ \theta_k
 * \f$ is the channel phase of the \f$ k^{th} \f$ channel tap.
 *
 * The average power and delay profiles are defined as:
 *
 * \f[ \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}] \f]
 *
 * and
 *
 * \f[ \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}],
 * \f]
 *
 * respectively. We assume without loss of generality that \f$ \tau_0 = 0
 * \f$ and \f$ \tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1} \f$.
 * Now the received signal is simply a linear filtering (or convolution)
 * of the transmitted signal, where \f$ h(t) \f$ is the impulse response
 * of the filter.
 *
 * In practice, when simulating a communication link, the impulse response
 * \f$ h(t) \f$ is sampled with a sample period \f$ T_s \f$ that is
 * related to the symbol rate of the system of investigation (often 2-8
 * times higher). Hence, the impulse response of the channel can now be
 * modeled as a time-discrete filter or tapped-delay line (TDL) where the
 * delays are given as \f$ \tau_k = d_k T_s \f$, and \f$ d_k \f$ are
 * positive integers.
 *
 *
 * \section channel_los Line of Sight (LOS) or Rice Fading
 *
 * If there is a line of sight (LOS) between the transmitter and receiver
 * the first component received (or a few first) will have a static
 * component that depends only on the Doppler frequency. In practice the
 * difference in time between the first LOS component(s) and the reflected
 * components is small and hence in a discretized system the first tap is
 * usually modeled as a LOS component and a fading component. Such a
 * process is usually called a Rice fading process.
 *
 * The LOS component can be expressed as:
 *
 * \f[ \rho \exp(2 \pi f_\rho t + \theta_\rho), \f]
 *
 * where \f$ \rho \f$, \f$ f_\rho \f$, and \f$ \theta_\rho \f$ are the
 * amplitude, Doppler frequency and phase of the LOS component,
 * respectively. Instead of stating the amplitude itself, the ratio of the
 * LOS power and total fading process power, called the Rice factor (or
 * relative power), is often used. The Doppler frequency is limited by the
 * maximum Doppler frequency \f$ f_\mathrm{max} \f$ and hence typically
 * the Doppler of the LOS is expressed relative to its maximum (common is
 * \f$ f_\rho = 0.7 f_\mathrm{max} \f$). The phase is usually assumed to
 * be random and can without loss of generality be set to 0 (does not
 * affect the statistics of the process).
 *
 *
 * \section channel_ref References
 *
 * - [Pat02] Matthias Patzold, Mobile fading channels, Wiley, 2002.
 * - [Stu01] Gordon L. Stuber, Principles of mobile communication, 2nd.
 * ed., Kluwer, 2001.
 * - [Rap96] Theodore S. Rappaport, Wireless communications: principles
 * and practise, Prentice Hall, 1996.
 *
 */


namespace itpp
{

//! \addtogroup channels
//@{

//! Predefined channel profiles. Includes LOS and Doppler spectrum settings.
enum CHANNEL_PROFILE {
  ITU_Vehicular_A, ITU_Vehicular_B, ITU_Pedestrian_A, ITU_Pedestrian_B,
  COST207_RA, COST207_RA6,
  COST207_TU, COST207_TU6alt, COST207_TU12, COST207_TU12alt,
  COST207_BU, COST207_BU6alt, COST207_BU12, COST207_BU12alt,
  COST207_HT, COST207_HT6alt, COST207_HT12, COST207_HT12alt,
  COST259_TUx, COST259_RAx, COST259_HTx
};

//! Fading generator type: Independent (default), Static or Correlated.
enum FADING_TYPE { Independent, Static, Correlated };

//! Correlated fading generation methods: Rice_MEDS (default), IFFT or FIR.
enum CORRELATED_METHOD { Rice_MEDS, IFFT, FIR };

//! Rice fading generation methods: MEDS.
enum RICE_METHOD { MEDS };

//! Predefined Doppler spectra
enum DOPPLER_SPECTRUM {
  Jakes = 0, J = 0, Classic = 0, C = 0,
  GaussI = 1, Gauss1 = 1, GI = 1, G1 = 1,
  GaussII = 2, Gauss2 = 2, GII = 2, G2 = 2
};


/*!
 * \brief Fading generator class
 * \author Tony Ottosson, Adam Piatyszek and Zbigniew Dlugaszewski
 *
 * Abstract base class defining the interface of a single tap fading
 * generator. Besides pure interface methods it implements a common
 * \a set_LOS_power() method for setting up the Rice factor to be used in
 * fading generators, which inherit from this class.
 */
class ITPP_EXPORT Fading_Generator
{
public:
  //! Default constructor
  Fading_Generator();
  //! Destructor
  virtual ~Fading_Generator() {}

  //! Set relative LOS power
  void set_LOS_power(double relative_power);
  //! Set relative Doppler of the LOS component (for correlated fading generators)
  virtual void set_LOS_doppler(double relative_doppler);
  //! Set time offset in samples (for correlated fading generators)
  virtual void set_time_offset(int offset);
  //! Set FIR filter length (for FIR fading generator)
  virtual void set_filter_length(int filter_length);
  //! Set normalized Doppler (for correlated fading generators)
  virtual void set_norm_doppler(double norm_doppler);
  //! Set Doppler spectrum (for Rice fading generator)
  virtual void set_doppler_spectrum(DOPPLER_SPECTRUM spectrum);
  //! Set number of sine frequencies (for Rice fading generator)
  virtual void set_no_frequencies(int no_freq);
  //! Set calculation method of Doppler frequencies and amplitudes (for Rice fading generator)
  virtual void set_rice_method(RICE_METHOD method);

  //! Get relative power of LOS component (Rice factor)
  double get_LOS_power() const { return los_power; }
  //! Get relative Doppler of the LOS component (for correlated fading generators)
  virtual double get_LOS_doppler() const;
  //! Get time offset in samples (for correlated fading generators)
  virtual double get_time_offset() const;
  //! Set FIR filter length (for FIR fading generator)
  virtual int get_filter_length() const;
  //! Return normalized Doppler (for correlated fading generators)
  virtual double get_norm_doppler() const;
  //! Return Doppler spectrum (for Rice fading generator)
  virtual DOPPLER_SPECTRUM get_doppler_spectrum() const;
  //! Get number of sine frequencies (for Rice fading generator)
  virtual int get_no_frequencies() const;
  //! Get calculation method of Doppler frequencies and amplitudes (for Rice fading generator)
  virtual RICE_METHOD get_rice_method() const;

  //! Shift generator time offset by a number of samples (for correlated fading generators)
  virtual void shift_time_offset(int no_samples);

  //! Initialize the generator
  virtual void init() = 0;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec &output) = 0;
  //! Generate \a no_samples values from the fading process
  cvec generate(int no_samples);

protected:
  bool init_flag; //!< signals if generator is initialized or not
  double los_power; //!< Relative power of LOS component compared to diffuse component (K factor)
  double los_diffuse; //!< Diffuse component: sqrt(1 / (1 + los_power))
  double los_direct; //!< Direct component: sqrt(los_power / (1 + los_power))
};

//! \cond

#if (!defined(_MSC_VER) || (defined(_MSC_VER) && defined (ITPP_SHARED_LIB)))
//define two input stream operators for data types involved
//otherwise explicit instantiation of Array is impossible
inline std::istream& operator>>(std::istream& is, Fading_Generator* &pfg)
{
  void* p;
  is >> p;
  if(is) pfg = (Fading_Generator*)p;
  return is;
}

inline std::istream& operator>>(std::istream& is, DOPPLER_SPECTRUM& sp)
{
  int val;
  is >> val;
  if(!is) return is;
  if(val > 0 && val <= G2) sp = (DOPPLER_SPECTRUM)val;
  else is.setstate(std::ios_base::failbit);
  return is;
}

//MSVC explicitly instantiate required template while building the shared library
#ifdef _WIN32
template class ITPP_EXPORT Array<DOPPLER_SPECTRUM>;
template class ITPP_EXPORT Array<Fading_Generator*>;
#endif
#endif

//! \endcond

/*!
 * \brief Independent (random) fading generator class
 * \author Adam Piatyszek
 *
 * This class implements the intependent fading generator, which can be
 * used on each tap of the TDL channel model. This generator produces a
 * set of independent Rayleigh or Rice distributed channel coefficients.
 */
class ITPP_EXPORT Independent_Fading_Generator : public Fading_Generator
{
public:
  //! Default constructor
  Independent_Fading_Generator() : Fading_Generator() {}
  //! Destructor
  virtual ~Independent_Fading_Generator() {}

  //! Initialize the generator
  virtual void init() { init_flag = true; }

  using Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec& output);
};


/*!
 * \brief Static fading generator class
 * \author Adam Piatyszek
 *
 * This class implements the static fading generator, which can be
 * used on each tap of the TDL channel model. This generator produces a
 * set of identical (static) Rayleigh or Rice distributed channel
 * coefficients.
 */
class ITPP_EXPORT Static_Fading_Generator : public Fading_Generator
{
public:
  //! Default constructor
  Static_Fading_Generator() : Fading_Generator() {}
  //! Destructor
  virtual ~Static_Fading_Generator() {}

  //! Initialize the generator
  virtual void init();

  using Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec& output);

protected:
  //! Static Rayleigh distributed sample
  double static_sample_re;
  double static_sample_im;
};


/*!
 * \brief Correlated (random) fading generator class
 * \author Adam Piatyszek
 *
 * Correlated fading generator class implements an abstract interface for
 * a set of generators. Parameter that define the correlated fading
 * process are the normalized Doppler and optionally the relative Doppler
 * of a LOS component.
 *
 * References:
 * - [Pat02] Matthias Patzold, Mobile fading channels, Wiley, 2002.
 */
class ITPP_EXPORT Correlated_Fading_Generator : public Fading_Generator
{
public:
  //! Default constructor
  Correlated_Fading_Generator(double norm_doppler);
  //! Destructor
  virtual ~Correlated_Fading_Generator() {}

  //! Set normalized Doppler
  virtual void set_norm_doppler(double norm_doppler);
  //! Set relative Doppler (compared to the maximum Doppler) for the LOS component
  virtual void set_LOS_doppler(double relative_doppler);
  //! Set time offset in samples
  virtual void set_time_offset(int offset);

  //! Return normalized Doppler
  virtual double get_norm_doppler() const { return n_dopp; }
  //! Get relative Doppler (compared to the maximum doppler) for the LOS component
  virtual double get_LOS_doppler() const { return los_dopp; }
  //! Get time offset in samples
  virtual double get_time_offset() const { return time_offset; }

  //! Shift generator time offset by a number of samples
  virtual void shift_time_offset(int no_samples);

  //! Initialize the generator
  virtual void init() = 0;

  using Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec& output) = 0;

protected:
  double n_dopp; //!< Normalized maximum Doppler frequency
  double los_dopp; //!< Relative Doppler on LOS component (0.7 by default)
  double time_offset; //!< Time offset in samples (time state in the generator)

  //! add LOS component to the \a sample with index \a idx
  void add_LOS(int idx, std::complex<double>& sample);
};


/*!
 * \brief Rice type fading generator class
 * \author Tony Ottosson, Adam Piatyszek and Zbigniew Dlugaszewski
 *
 * A Rice generator is a generator of the form:
 *
 * \f[ \tilde \mu_i(t) = \sum_{n=1}^{N_i} c_{i,n} \cos(2\pi f_{i,n} t +
 * \theta_{i,n}) \f]
 *
 * Here \f$ c_{i,n} \f$, \f$ f_{i,n} \f$, and \f$ \theta_{i,n} \f$ are the
 * Doppler coefficients, discrete Doppler frequencies, and Doppler phases,
 * respectively. Rice showed that a generator of this form can perfectly
 * model a Gaussian process when \f$ N_i \rightarrow \infty \f$. When
 * generating a fading pattern we need a complex-valued generator
 *
 * \f[ \tilde \mu(t) = \tilde \mu_1(t) + j \tilde \mu_2(t) \f]
 *
 * Parameters that define the generator are the normalized Doppler and the
 * doppler spectrum. Possible values of the Doppler spectrum are:
 * - <em>Jakes, J, Classic, C</em>: The classical Jakes spectrum shape.
 * This method can be also used for modelling LOS paths, by setting the
 * LOS relative power (Rice factor) and LOS Doppler parameters.
 * - <em>GaussI, GI, Gauss1, G1</em>
 * - <em>GaussII, GII, Gauss2, G2</em>
 *
 * Furthermore also the number of sine waves, \f$ N_i \f$ and method used
 * to calculate the parameters \f$ c_{i,n} \f$, \f$ f_{i,n} \f$, and \f$
 * \theta_{i,n} \f$ can be specified. For now the only method defined for
 * calculating the parameters is the Method of Exact Doppler Spread
 * (MEDS). See [Pat02] for more details.
 *
 * References:
 * - [Pat02] Matthias Patzold, Mobile fading channels, Wiley, 2002.
 */
class ITPP_EXPORT Rice_Fading_Generator : public Correlated_Fading_Generator
{
public:
  //! Default constructor
  Rice_Fading_Generator(double norm_doppler, DOPPLER_SPECTRUM spectrum = Jakes,
                        int no_freq = 16, RICE_METHOD method = MEDS);
  //! Destructor
  virtual ~Rice_Fading_Generator() {}

  //! Set Doppler spectrum
  virtual void set_doppler_spectrum(DOPPLER_SPECTRUM spectrum);
  //! Set number of Doppler frequencies
  virtual void set_no_frequencies(int no_freq);
  //! Set calculation method of Doppler frequencies and amplitudes
  virtual void set_rice_method(RICE_METHOD method);

  //! Return Doppler spectrum
  virtual DOPPLER_SPECTRUM get_doppler_spectrum() const { return dopp_spectrum; }
  //! Get number of Doppler frequencies
  virtual int get_no_frequencies() const { return Ni; }
  //! Get calculation method of Doppler frequencies and amplitudes
  virtual RICE_METHOD get_rice_method() const { return rice_method; }

  //! Initialize the generator
  virtual void init();

  using Correlated_Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec &output);

protected:
  DOPPLER_SPECTRUM dopp_spectrum; //!< Doppler spectrum type (Jakes by default)
  //! Number of sine waves in a Gaussian process
  int Ni;
  //! Rice process generation method
  RICE_METHOD rice_method;
  /*! Doppler frequencies, amplitudes and phases
   * @{ */
  vec f1, f2, c1, c2, th1, th2;
  /*! @} */
  /*! Frequency shift values of the Doppler spectrum in GaussI and GaussII
   * @{ */
  double f01, f02;
  /*! @} */

  //! Init function for MEDS method
  void init_MEDS();
};


/*!
 * \brief FIR type Fading generator class
 * \author Tony Ottosson and Adam Piatyszek
 *
 * A FIR generator is a linear finite impulse response (FIR) filter
 * implementation of a filter method for generation of a fading process.
 * Parameters that define the generator are the normalized Doppler and
 * length of the FIR filter. The default value of filter length is 500. If
 * the normalized Doppler frequency is lower than 0.1 an equivalent
 * process of a higher normalized Doppler is generated and linearly
 * interpolated.
 *
 * References:
 * - [Stu01] Gordon L. Stuber, Principles of mobile communication, 2nd.
 * ed., Kluwer, 2001.
 * - [Rap96] Theodore S. Rappaport, Wireless communications: principles
 * and practise, Prentice Hall, 1996.
 */
class ITPP_EXPORT FIR_Fading_Generator : public Correlated_Fading_Generator
{
public:
  //! Default constructor
  FIR_Fading_Generator(double norm_doppler, int filter_length = 500);
  //! Destructor
  virtual ~FIR_Fading_Generator() {}

  //! Set FIR filter length
  virtual void set_filter_length(int filter_length);
  //! Get filter length
  virtual int get_filter_length() const { return fir_length; }

  //! Initialize the generator
  virtual void init();

  using Correlated_Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec &output);

protected:
  int fir_length; //!< Size of FIR filter
  int upsample_rate; //!< Upsampling rate for linear interpolation
  //! Filter used for fading generation
  MA_Filter<std::complex<double>, double, std::complex<double> > fir_filter;
  cvec left_overs; //!< Left-overs from upsampling

  /*!
   * \brief Jakes spectrum filter
   *
   * Function that generates the taps in the Jakes filter.
   * \param order Number of taps in the filter
   * \param norm_dopp Normalized Doppler frequency, i.e. \f$ f_{norm} =
   * f_{max} T_{s} \f$, where \f$ f_{max} \f$ is the actual Doppler
   * frequency and \f$ T_{s} \f$ is the sampling interval.
   * \return A vector containing the filter taps of the Jakes filter.
   */
  vec Jakes_filter(double norm_dopp, int order = 100);
};


/*!
 * \brief IFFT type Fading generator class
 * \author Tony Ottosson and Adam Piatyszek
 *
 * A IFFT generator is a frequency domain filter implementation of filter
 * method for generation of a fading process. Parameters that define the
 * generator is the normalized Doppler.
 *
 * The method is block-based and consecutive blocks are independent
 * fading. The method is very efficient for large blocks. The size of the
 * FFT, \f$ N_\mathrm{fft} \f$, is given by the nearest higher integer
 * power of two that is larger than \a no_samples. For small blocks the
 * FFT size is increased to keep a good accuracy (at least 10 samples in
 * the representation of the Doppler spectrum). However, to keep the
 * program reasonably efficient the largest upsampling factor is 64.
 * Higher factors will result in a run-time error. If so, please use
 * another method.
 *
 * References:
 * - [Stu01] Gordon L. Stuber, Principles of mobile communication, 2nd.
 * ed., Kluwer, 2001.
 * - [Rap96] Theodore S. Rappaport, Wireless communications: principles
 * and practise, Prentice Hall, 1996.
 */
class ITPP_EXPORT IFFT_Fading_Generator : public Correlated_Fading_Generator
{
public:
  //! Default constructor
  IFFT_Fading_Generator(double norm_doppler) :
    Correlated_Fading_Generator(norm_doppler) {}
  //! Destructor
  virtual ~IFFT_Fading_Generator() {}

  //! Initialize the generator
  virtual void init() { init_flag = true; }

  using Correlated_Fading_Generator::generate;

  //! Generate \a no_samples values from the fading process
  virtual void generate(int no_samples, cvec &output);

protected:
  //! Generator for Jakes spectrum
  void generate_Jakes(int no_samples, cvec &output);
};


/*!
 * \brief General specification of a time-domain multipath channel
 * \author Tony Ottosson and Adam Piatyszek
 *
 * This class does NOT generate any channel values. It is only used to
 * specify the channel model. To generate channel coefficients use the
 * Tapped-Delay Line (TDL) class \a TDL_Channel.
 *
 * A time invariant (or at least wide-sense stationary) channel have an
 * impulse response that can be modeled as:
 *
 * \f[ h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k)
 * \delta(t-\tau_k), \f]
 *
 * where \f$ N_{taps} \f$ is the number of channel taps, \f$ a_k \f$ is
 * the average amplitude at delay \f$ \tau_k \f$, and \f$ \theta_k \f$ is
 * the channel phase of the \f$ k^{th} \f$ channel tap. The average power
 * profile, and the delay profiles are defined as:
 *
 * \f[ \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}] \f]
 *
 * and
 *
 * \f[ \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}],
 * \f]
 *
 * respectively. We assume without loss of generality that \f$ \tau_0 = 0
 * \f$ and \f$ \tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1} \f$.
 *
 * To initialize the class the following parameters should be defined:
 * - \a avg_power_dB Average power profile \f$ \mathbf{a} \f$, given in dB
 * - \a delay_profile Delay profile \f$ \mathbf{\tau} \f$, given in seconds.
 * The delay profile should be sorted (increasing delay) and the first
 * delay has to be 0.
 *
 * Optionally one can define LOS parameters: \a relative_power and
 * \a relative_doppler, and additionally the kind of Doppler spectrum for
 * each tap.
 *
 * It is also possible to specify a predefined channel model. The
 * implemented ones are as follows:
 * - ITU Channel models:
 *   - \a ITU_Vehicular_A, 6-tap channel
 *   - \a ITU_Vehicular_B, 6-tap channel
 *   - \a ITU_Pedestrian_A, 4-tap channel
 *   - \a ITU_Pedestrian_B, 6-tap channel
 * - COST 207 Channel models (see [Pat02] pp. 259-266):
 *   - \a COST207_RA, Rural area, 4-tap channel
 *   - \a COST207_RA6, Rural area, 6-tap channel
 *   - \a COST207_TU, Typical urban, 6-tap channel
 *   - \a COST207_TU6alt, Typical urban, alternative 6-tap channel
 *   - \a COST207_TU12, Typical urban, 12-tap channel
 *   - \a COST207_TU12alt, Typical urban, alternative 12-tap channel
 *   - \a COST207_BU, Bad urban, 6-tap channel
 *   - \a COST207_BU6alt, Bad urban, alternative 6-tap channel
 *   - \a COST207_BU12, Bad urban, 12-tap channel
 *   - \a COST207_BU12alt, Bad urban, alternative 12-tap channel
 *   - \a COST207_HT, Hilly terrain, 6-tap channel
 *   - \a COST207_HT6alt, Hilly terrain, alternative 6-tap channel
 *   - \a COST207_HT12, Hilly terrain, 12-tap channel
 *   - \a COST207_HT12alt, Hilly terrain, alternative 12-tap channel
 * - COST 259 Channel models (see [3GPP_TR_25.943]):
 *   - \a COST259_TUx, Typical urban, 20-tap channel
 *   - \a COST259_RAx, Rural ara, 10-tap channel
 *   - \a COST259_HTx, Hilly terrain, 20-tap channel
 *
 * References:
 * - [Pat02] Matthias Patzold, Mobile fading channels, Wiley, 2002.
 * - [3GPP_TR_25.943] Technical Specification Group Radio Access Networs;
 * Deployment aspects. Version 5.1.0 (2002-06).
 */
class ITPP_EXPORT Channel_Specification
{
public:
  //! Default constructor (power profile in dB, delay profile in seconds)
  Channel_Specification(const vec &avg_power_dB = "0", const vec &delay_prof = "0");
  //! Initialize with predetermined channel profile
  Channel_Specification(const CHANNEL_PROFILE profile);
  //! Destructor
  virtual ~Channel_Specification() {}

  //! Set both average power profile in dB and power delay profile in seconds
  void set_channel_profile(const vec &avg_power_dB, const vec &delay_prof);
  //! Set channel profile to a predetermined profile
  void set_channel_profile(const CHANNEL_PROFILE profile);

  //! Set doppler spectrum for each tap in the channel profile
  void set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum);
  //! Set doppler spectrum for tap \a tap_number in the channel profile.
  void set_doppler_spectrum(int tap_number, DOPPLER_SPECTRUM tap_spectrum);

  //! Set LOS (Rice) components for tap \a tap_number
  void set_LOS(int tap_number, double relative_power, double relative_doppler = 0.7);
  //! Set LOS (Rice) components for all taps
  void set_LOS(const vec& relative_power, const vec& relative_doppler = "");

  //! Get both average power profile in dB and power delay profile in seconds
  void get_channel_profile(vec &avg_power_dB, vec &delay_prof) const;
  //! Return power profile in dB
  vec get_avg_power_dB() const { return a_prof_dB; }
  //! Return delay profile in seconds
  vec get_delay_prof() const { return d_prof; }
  //! Get doppler spectrum for tap \a index
  Array<DOPPLER_SPECTRUM> get_doppler_spectrum() const { return tap_doppler_spectrum; }
  //! Get doppler spectrum for tap \a index
  DOPPLER_SPECTRUM get_doppler_spectrum(int index) const;
  //! Get relative power (Rice factor) for each tap
  vec get_LOS_power() const { return los_power; }
  //! Get relative Doppler for each tap
  vec get_LOS_doppler() const { return los_dopp; }
  //! Get relative power (Rice factor) for tap \a tap_number
  double get_LOS_power(int tap_number) const { return los_power(tap_number); }
  //! Get relative Doppler for tap \a tap_number
  double get_LOS_doppler(int tap_number) const { return los_dopp(tap_number); }

  //! Return the number of channel taps
  int taps() const { return N_taps; }

  //! Calculate mean excess delay in samples
  double calc_mean_excess_delay() const;
  //! Calculate RMS delay spread in samples
  double calc_rms_delay_spread() const;

protected:
  vec a_prof_dB; //!< Power profile in dB
  vec d_prof; //!< Delay profile in seconds
  Array<DOPPLER_SPECTRUM> tap_doppler_spectrum; //!< Doppler spectrum for each tap
  int N_taps; //!< Number of taps
  vec los_power; //!< Relative power for each Rice component
  vec los_dopp; //!< Relative Rice Doppler for each Rice component
};


/*!
 * \brief Tapped Delay Line (TDL) channel model
 * \author Tony Ottosson, Adam Piatyszek and Zbigniew Dlugaszewski
 *
 * A time invariant (or at least wide-sense stationary) channel have an
 * impulse response that can be modeled as:
 *
 * \f[ h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k )
 * \delta(t-\tau_k), \f]
 *
 * where \f$ N_{taps} \f$ is the number of channel taps, \f$ a_k \f$ is
 * the average amplitude at delay \f$ \tau_k \f$, and \f$ \theta_k \f$ is
 * the channel phase of the \f$ k^{th} \f$ channel tap. The average power
 * profile, and delay profile are defined as:
 *
 * \f[ \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}] \f]
 *
 * and
 *
 * \f[ \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}],
 * \f]
 *
 * respectively. We assume without loss of generality that \f$ \tau_0 =
 * 0 \f$ and \f$ \tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1} \f$.
 *
 * To initialize the channel profile the following parameters should be
 * defined:
 * - \a avg_power_dB Average power profile \f$ \mathbf{a} \f$, given in dB
 * - \a delay_profile Delay profile \f$ \mathbf{\tau} \f$, given in samples.
 * The delay profile should be sorted (increasing delay) and the first
 * delay has to be 0.
 *
 * In the case of correlated channel, the correlation in time is decided
 * by the normalized Doppler frequency and Doppler spectrum. The
 * normalized Doppler frequency is calculated as \f$ f_\mathrm{max} T_s
 * \f$, where \f$ f_\mathrm{max} \f$ is the maximum Doppler frequency and
 * \f$ T_s \f$ is the sample duration.
 *
 * Two main types of generation methods exist: the filter method and Rice
 * method. In the filter method the correlated fading process is generated
 * with a filtering of the complex Gaussian process to achieve a given
 * Doppler spectrum (\a Jakes by default). Currently there are two filter
 * implementations:
 * - \a FIR - Finite Impulse Response (FIR) filter. Rather slow and
 * inaccurate. Use with care.
 * - \a IFFT - Blockwise filtering in the frequency-domain. Consecutive
 * blocks have independent fading. Very efficient method for large blocks.
 *
 * The preferred method for generating the correlated fading process is
 * the \a Rice method that approximate the fading process as a sum of
 * sinusoids. Currently there is only one implementation, the Rice Method
 * of Exact Doppler Spread (\a Rice_MEDS), which is also the default choice.
 *
 * To summarize, the currently supported correlated fading generation
 * methods are:
 * - \a FIR - filter FIR method
 * - \a IFFT - filter IFFT method
 * - \a Rice_MEDS - sum of sine waves MEDS method
 *
 * Beside the above described correlated methods, two additional fading
 * generators exists: the \a Independent_Fading_Generator and
 * \a Static_Fading_Generator ones. Their names are self-explanatory. The
 * only optional parameter of these two generators is a \a relative_power
 * of the LOS component.
 *
 * Example: Simulation of WCDMA
 * \code
 * #include <itpp/itcomm.h>
 * using namespace itpp;
 *
 * int main() {
 *   // set sampling time at a half of chip rate (0.5 / 3.84e6)
 *   double Ts = 130.2e-9;
 *   // select the COST259 Rural Area model
 *   Channel_Specification channel_spec(COST259_RAx);
 *   // initialize with the predefined channel profile
 *   TDL_Channel channel(channel_spec, Ts);
 *   // set the normalized Doppler; fading type will be set to Correlated
 *   // and Rice_MEDS method will be used (default settings)
 *   channel.set_norm_doppler(0.01);
 *
 *   cvec transmitted_signal;
 *   // -----------------------------------------------
 *   // Your code here generates the transmitted signal
 *   // -----------------------------------------------
 *
 *   // Channel coefficients are returned in the 'coeff' array of complex values
 *   Array<cvec> coeff;
 *   cvec received_signal = channel(transmitted_signal, coeff);
 * }
 * \endcode
 *
 * References:
 * - [Pat02] Matthias Patzold, Mobile fading channels, Wiley, 2002.
 * - [Stu01] Gordon L. Stuber, Principles of mobile communication, 2nd.
 * ed., Kluwer, 2001.
 * - [Rap96] Theodore S. Rappaport, Wireless communications: principles
 * and practise, Prentice Hall, 1996.
 */
class ITPP_EXPORT TDL_Channel
{
public:
  //! Default constructor
  TDL_Channel(const vec &avg_power_dB = "0", const ivec &delay_prof = "0");
  //! Constructor with defined Channel_Specification. Delay profile will be discretized with \a sampling_time in seconds.
  TDL_Channel(const Channel_Specification &channel_spec, double sampling_time);
  //! Destructor
  virtual ~TDL_Channel();

  //! Set both average power profile in dB and power delay profile in samples
  void set_channel_profile(const vec &avg_power_dB, const ivec &delay_prof);
  //! Set channel profile to uniform with \a no_taps taps
  void set_channel_profile_uniform(int no_taps);
  //! Set channel profile to exponential with \a no_taps taps
  void set_channel_profile_exponential(int no_taps);
  //! Set channel profile using Channel_Specification. Delay profile will be discretized with \a sampling_time in seconds.
  void set_channel_profile(const Channel_Specification &channel_spec, double sampling_time);

  //! Set the fading generation method to \a method
  void set_correlated_method(CORRELATED_METHOD method);
  //! Set fading type to one of \a Independent, \a Static or \a Correlated
  void set_fading_type(FADING_TYPE fading_type);

  //! Set normalized Doppler rate. A \a Correlated fading type will be used.
  void set_norm_doppler(double norm_doppler);

  //! Set LOS parameters for each tap. LOS Doppler will be set to 0.7 by default.
  void set_LOS(const vec& relative_power, const vec& relative_doppler = "");
  //! Set LOS power for each tap. LOS Doppler will be set to 0.7 by default.
  void set_LOS_power(const vec& relative_power);
  //! Set LOS doppler for each tap. A \a Correlated fading type will be used.
  void set_LOS_doppler(const vec& relative_doppler);

  //! Set doppler spectrum for each tap in the channel profile. \a Rice_MEDS method will be used.
  void set_doppler_spectrum(const DOPPLER_SPECTRUM *tap_spectrum);
  //! Set doppler spectrum of tap \a tap_number. \a Rice_MEDS method will be used.
  void set_doppler_spectrum(int tap_number, DOPPLER_SPECTRUM tap_spectrum);
  //! Set number of sine frequencies. \a Rice_MEDS method will be used.
  void set_no_frequencies(int no_freq);

  //! Set fading generators' time offset in samples. A \a Correlated fading type will be used.
  void set_time_offset(int offset);
  //! Shift fading generators' time offset. A \a Correlated fading type will be used.
  void shift_time_offset(int no_samples);

  //! Set fading generator filter length. \a FIR method will be used.
  void set_filter_length(int filter_length);


  //! Return the number of channel taps
  int taps() const { return N_taps; }

  //! Get both average power profile in dB and power delay profile in samples
  void get_channel_profile(vec &avg_power_dB, ivec &delay_prof) const;
  //! Return power profile in dB
  vec get_avg_power_dB() const;
  //! Return delay profile in samples
  ivec get_delay_prof() const { return d_prof; }

  //! Return fading generation method
  CORRELATED_METHOD get_correlated_method() const { return method; }
  //! Return fading type
  FADING_TYPE get_fading_type() const { return fading_type; }

  //! Return normalized doppler rate
  double get_norm_doppler() const { return n_dopp; }

  //! Get relative power (Rice factor) for each tap
  vec get_LOS_power() const { return los_power; }
  //! Get relative Doppler (to the maximum Doppler) for each tap
  vec get_LOS_doppler() const { return los_dopp; }
  //! Get relative power (Rice factor) for tap \a tap_number
  double get_LOS_power(int tap_number) const { return los_power(tap_number); }
  //! Get relative Doppler (to the maximum Doppler) for tap \a tap_number
  double get_LOS_doppler(int tap_number) const { return los_dopp(tap_number); }

  //! Get the minimum number of frequencies used in Rice MEDS fading generator
  int get_no_frequencies() const { return nrof_freq; }

  //! Get fading generators' time ofset
  double get_time_offset() const;

  //! Calculate mean excess delay in samples
  double calc_mean_excess_delay() const;
  //! Calculate RMS delay spread in samples
  double calc_rms_delay_spread() const;

  //! Initialize all fading generators. Automatically invoked in generate() or filter() functions.
  void init();

  //! Generate \a no_samples values of the channel
  void generate(int no_samples, Array<cvec> &channel_coeff);
  //! Generate \a no_samples values of the channel. Returns the matrix with one tap per column.
  void generate(int no_samples, cmat &channel_coeff);

  //! Filter the \a input with the known channel values \a channel_coeff (e.g. from the generate function)
  void filter_known_channel(const cvec &input, cvec &output, const Array<cvec> &channel_coeff);
  //! Filter the \a input with the known channel values \a channel_coeff (e.g. from the generate function)
  void filter_known_channel(const cvec &input, cvec &output, const cmat &channel_coeff);

  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  void filter(const cvec &input, cvec &output, Array<cvec> &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  void filter(const cvec &input, cvec &output, cmat &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  cvec filter(const cvec &input, Array<cvec> &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  cvec filter(const cvec &input, cmat &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Only return the output
  void filter(const cvec &input, cvec &output);
  //! Generate channel coefficients and filter the \a input. Only return output
  cvec filter(const cvec &input);

  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  void operator()(const cvec &input, cvec &output, Array<cvec> &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  void operator()(const cvec &input, cvec &output, cmat &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  cvec operator()(const cvec &input, Array<cvec> &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Return output and channel coefficients
  cvec operator()(const cvec &input, cmat &channel_coeff);
  //! Generate channel coefficients and filter the \a input. Only return output
  cvec operator()(const cvec &input);

  //! Calculate impulse-response on the supplied channel coefficients (produced by the generate() function)
  void calc_impulse_response(const Array<cvec> &channel_coeff, Array<cvec> &impulse_response);

  //! Calculate frequency-response on the supplied channel coefficients (produced by the generate() function)
  void calc_frequency_response(const Array<cvec> &channel_coeff, Array<cvec> &frequency_response, const int fft_size);
  //! Calculate frequency-response on the supplied channel coefficients (produced by the generate() function)
  void calc_frequency_response(const cmat &channel_coeff, cmat &frequency_response, const int fft_size);
  //! Return channel sampling time (used for discretization)
  double get_sampling_time() const { return discrete_Ts; }

protected:
  bool init_flag; //!< Channel ready to produce data
  vec a_prof; //!< Average amplitude of each tap
  ivec d_prof; //!< Delay in samples for each tap
  vec los_power; //!< Relative power for each Rice component
  vec los_dopp; //!< Relative LOS Doppler for each Rice component
  int N_taps; //!< Number of taps
  double n_dopp; //!< Normalized Doppler of the correlated fading
  FADING_TYPE fading_type; //!< Fading type: Independent (default), Static or Correlated
  CORRELATED_METHOD method; //!< Correlated fading generation method: Rice_MEDS (default), IFFT or FIR
  Array<DOPPLER_SPECTRUM> tap_doppler_spectrum; //!< Doppler spectrum for each tap
  Array<Fading_Generator *> fading_gen; //!< Fading generators for each tap
  int filter_length; //!< Filter length of FIR fading generator
  int nrof_freq; //!< Number of sine frequencies in the Rice MEDS fading generator
  double discrete_Ts; //!< Sampling time of discretization

  /*!
   * \brief Discretize the delay profile with \a discrete_Ts (Ts). All
   * taps within ((i-0.5)Ts,(i+0.5)Ts] belong to the ith discrete tap.
   *
   * \param delay_profile Delay profile in seconds.
   */
  void discretize(const vec &delay_profile);
};



/*!
  \brief A Binary Symetric Channel with crossover probability p.

  Input and output are of type \a bvec with 0 and 1.
  Example:
  \code
  #include <itpp/itcomm.h>
  using namespace itpp;

  int main() {
    // Initiate the BSC with cross-over probability 0.1
    BSC bsc(0.1);

    bvec transmitted_bits = randb(100);
    bvec received_bits = bsc(transmitted_bits);
  }
  \endcode
*/
class ITPP_EXPORT BSC
{
public:
  //! Class constructor. Sets the error probability to p
  BSC(double in_p = 0.0) : u(0.0, 1.0) { p = in_p; };
  //! Set crossover (bit error) probability
  void set_prob(double in_p) { p = in_p; };
  //! Get crossover (bit error) probability
  double get_prob() const { return p; };
  //! Feed \a input through the BSC channel
  bvec operator()(const bvec &input);
private:
  Uniform_RNG u;
  double p;
};



/*!
  \brief Ordinary AWGN Channel for cvec or vec inputs and outputs.

  For real signals, the input parameter (\a noisevar) denotes the noise
  variance per real dimension. Therefore, it should be set to \f$N_0/2\f$,
  where \f$N_0\f$ is the noise power spectral density. However, in case of
  complex signals, the input parameter (\a noisevar) represents the
  noise variance per complex dimension, i.e. the sum of the variances in
  the real and imaginary parts, and thus is equal to \f$N_0\f$.

  Example:
  \code
  #include <itpp/itcomm.h>
  using namespace itpp;

  int main() {
    // Initiate the AWGN_Channel class
    double noisevar = 0.1;
    AWGN_Channel awgn_channel(noisevar);

    // Initiate a QPSK-modulator, and generate the transmitted signal
    QPSK qpsk;
    bvec transmitted_bits = randb(20);
    cvec transmitted_signal = qpsk.modulate_bits(transmitted_bits);

    // Usage of the member operator ()
    cvec received_signal = awgn_channel(transmitted_signal);

    // Demodulate the bits
    bvec received_bits = qpsk.demodulate_bits(received_signal);
  }
  \endcode
*/
class ITPP_EXPORT AWGN_Channel
{
public:
  //! Class constructor. Sets the noise variance (for complex-valued channels the sum of real and imaginary parts)
  AWGN_Channel(double noisevar = 0.0): sigma(std::sqrt(noisevar)) {}
  //! Set noise variance (for complex-valued channels the sum of real and imaginary parts)
  void set_noise(double noisevar) { sigma = std::sqrt(noisevar); }
  //! Get noise variance (for complex-valued channels the sum of real and imaginary parts)
  double get_noise() const { return sqr(sigma); }
  //! Feed the complex input \a input through the complex-valued AWGN channel
  cvec operator()(const cvec &input);
  //! Feed the input \a through the real-valued AWGN channel
  vec operator()(const vec &input);
private:
  Complex_Normal_RNG rng_cn;
  Normal_RNG rng_n;
  double sigma;
};

//@}

} // namespace itpp

#endif // #ifndef CHANNEL_H
