/*!
 * \file
 * \brief Communication Channels' classes - source file
 * \author Tony Ottosson, Pal Frenger, Adam Piatyszek and Zbigniew Dlugaszewski
 *
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

#include <itpp/comm/channel.h>
#include <itpp/base/math/error.h>
#include <itpp/base/math/trig_hyp.h>
#include <itpp/base/bessel.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/signal/resampling.h>
#include <itpp/signal/transforms.h>
#include <itpp/signal/window.h>
#include <itpp/base/math/min_max.h>
#include <itpp/stat/misc_stat.h>


namespace itpp
{


// --------------------------------------------------------------------------
// Fading_Generator class
// --------------------------------------------------------------------------

Fading_Generator::Fading_Generator() : init_flag(false)
{
  // no default LOS component
  set_LOS_power(0.0);
}

void Fading_Generator::set_LOS_power(double relative_power)
{
  it_assert(relative_power >= 0.0,
            "Fading_Generator::set_LOS_power(): Relative_power can not be negative");
  los_power = relative_power;
  los_diffuse = std::sqrt(1.0 / (1.0 + los_power));
  los_direct = los_diffuse * std::sqrt(los_power);
}

void Fading_Generator::set_LOS_doppler(double)
{
  it_warning("Fading_Generator::set_LOS_doppler(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_time_offset(int)
{
  it_warning("Fading_Generator::set_time_offset(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_norm_doppler(double)
{
  it_warning("Fading_Generator::set_norm_doppler(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_filter_length(int)
{
  it_warning("Fading_Generator::set_filter_length(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_doppler_spectrum(DOPPLER_SPECTRUM)
{
  it_warning("Fading_Generator::set_doppler_spectrum(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_no_frequencies(int)
{
  it_warning("Fading_Generator::set_no_frequencies(): This function has no effect on this kind of generator");
}

void Fading_Generator::set_rice_method(RICE_METHOD)
{
  it_warning("Fading_Generator::set_rice_method(): This function has no effect on this kind of generator");
}

double Fading_Generator::get_LOS_doppler() const
{
  it_warning("Fading_Generator::get_LOS_doppler(): This function has no effect on this kind of generator");
  return 0;
}

double Fading_Generator::get_time_offset() const
{
  it_warning("Fading_Generator::get_time_offset(): This function has no effect on this kind of generator");
  return 0;
}

int Fading_Generator::get_filter_length() const
{
  it_warning("Fading_Generator::get_filter_length(): This function has no effect on this kind of generator");
  return 0;
}

double Fading_Generator::get_norm_doppler() const
{
  it_warning("Fading_Generator::get_norm_doppler(): This function has no effect on this kind of generator");
  return 0;
}

DOPPLER_SPECTRUM Fading_Generator::get_doppler_spectrum() const
{
  it_warning("Fading_Generator::get_doppler_spectrum(): This function has no effect on this kind of generator");
  return Jakes;
}

int Fading_Generator::get_no_frequencies() const
{
  it_warning("Fading_Generator::get_no_frequencies(): This function has no effect on this kind of generator");
  return 0;
}

RICE_METHOD Fading_Generator::get_rice_method() const
{
  it_warning("Fading_Generator::get_rice_method(): This function has no effect on this kind of generator");
  return MEDS;
}

void Fading_Generator::shift_time_offset(int)
{
  it_warning("Fading_Generator::shift_time_offset(): This function has no effect on this kind of generator");
}

cvec Fading_Generator::generate(int no_samples)
{
  cvec output;
  this->generate(no_samples, output);
  return output;
}


// --------------------------------------------------------------------------
// Independent_Fading_Generator class
// --------------------------------------------------------------------------

void Independent_Fading_Generator::generate(int no_samples, cvec& output)
{
  output.set_size(no_samples, false);
  if (los_power > 0.0) {
    for (int i = 0; i < no_samples; ++i) {
      output(i) = los_diffuse * randn_c() + los_direct;
    }
  }
  else {
    output = randn_c(no_samples);
  }
}


// --------------------------------------------------------------------------
// Static_Fading_Generator class
// --------------------------------------------------------------------------

void Static_Fading_Generator::init()
{
  std::complex<double> static_sample = randn_c();
  if (los_power > 0.0) {
    static_sample *= los_diffuse;
    static_sample += los_direct;
  }
  static_sample_re = static_sample.real();
  static_sample_im = static_sample.imag();
  init_flag = true;
}

void Static_Fading_Generator::generate(int no_samples, cvec& output)
{
  if (init_flag == false)
    init();

  output.set_size(no_samples, false);
  output = std::complex<double>(static_sample_re,static_sample_im);
}


// --------------------------------------------------------------------------
// Correlated_Fading_Generator class
// --------------------------------------------------------------------------

Correlated_Fading_Generator::Correlated_Fading_Generator(double norm_doppler) :
    Fading_Generator(), los_dopp(0.7), time_offset(0.0)
{
  set_norm_doppler(norm_doppler);
}

void Correlated_Fading_Generator::set_norm_doppler(double norm_doppler)
{
  it_assert((norm_doppler > 0) && (norm_doppler <= 1.0),
            "Correlated_Fading_Generator: Normalized Doppler out of range");
  n_dopp = norm_doppler;
  init_flag = false;
}

void Correlated_Fading_Generator::set_LOS_doppler(double relative_doppler)
{
  it_assert((relative_doppler >= 0) && (relative_doppler <= 1.0),
            "Correlated_Fading_Generator::set_LOS_doppler(): Relative Doppler out of range");
  los_dopp = relative_doppler;
}

void Correlated_Fading_Generator::set_time_offset(int offset)
{
  time_offset = static_cast<double>(offset);
}

void Correlated_Fading_Generator::shift_time_offset(int no_samples)
{
  time_offset += static_cast<double>(no_samples);
}

void Correlated_Fading_Generator::add_LOS(int idx, std::complex<double>& sample)
{
  double tmp_arg = m_2pi * los_dopp * n_dopp * (idx + time_offset);
  sample *= los_diffuse;
  sample += los_direct * std::complex<double>(std::cos(tmp_arg),
            std::sin(tmp_arg));
}


// --------------------------------------------------------------------------
// Rice_Fading_Generator class
// --------------------------------------------------------------------------

Rice_Fading_Generator::Rice_Fading_Generator(double norm_doppler,
    DOPPLER_SPECTRUM spectrum,
    int no_freq, RICE_METHOD method) :
    Correlated_Fading_Generator(norm_doppler)
{
  set_doppler_spectrum(spectrum);
  set_no_frequencies(no_freq);
  set_rice_method(method);
}

void Rice_Fading_Generator::set_doppler_spectrum(DOPPLER_SPECTRUM spectrum)
{
  dopp_spectrum = spectrum;
  init_flag = false;
}

void Rice_Fading_Generator::set_no_frequencies(int no_freq)
{
  it_assert(no_freq >= 7,
            "Rice_Fading_Generator::set_no_frequencies(): Too low number of Doppler frequencies");
  Ni = no_freq;
  init_flag = false;
}

void Rice_Fading_Generator::set_rice_method(RICE_METHOD method)
{
  // check if this method works for the given spectrum
  rice_method = method;
  init_flag = false;
}

void Rice_Fading_Generator::init()
{
  switch (rice_method) {
  case MEDS: // Method of Exact Doppler Spread (MEDS)
    init_MEDS();
    break;
  default:
    it_error("Rice_Fading_Generator::init(): Wrong Rice method for this fading generator");
  };

  init_flag = true; // generator ready to use
}

void Rice_Fading_Generator::generate(int no_samples, cvec &output)
{
  if (init_flag == false)
    init();

  output.set_size(no_samples, false);

  switch (dopp_spectrum) {
  case Jakes: {
    double tmp_re, tmp_im;
    if (los_power > 0.0) { // LOS component exists
      for (int i = 0; i < no_samples; i++) {
        tmp_re = sum(elem_mult(c1, cos(m_2pi * f1 * n_dopp * (i + time_offset) + th1)));
        tmp_im = sum(elem_mult(c2, cos(m_2pi * f2 * n_dopp * (i + time_offset) + th2)));
        output(i) = std::complex<double>(tmp_re, tmp_im);
        add_LOS(i, output(i));
      }
    }
    else {
      for (int i = 0; i < no_samples; i++) {
        tmp_re = sum(elem_mult(c1, cos(m_2pi * f1 * n_dopp * (i + time_offset) + th1)));
        tmp_im = sum(elem_mult(c2, cos(m_2pi * f2 * n_dopp * (i + time_offset) + th2)));
        output(i) = std::complex<double>(tmp_re, tmp_im);
      }
    }
    break;
  }
  case GaussI:
  case GaussII: {
    double tmp;
    for (int i = 0; i < no_samples; i++) {
      tmp = m_2pi * n_dopp * (i + time_offset);
      output(i) = sum(elem_mult(c1, cos(f1 * tmp + th1)))
                  * std::complex<double>(std::cos(f01 * tmp), -std::sin(f01 * tmp))
                  + sum(elem_mult(c2, cos(f2 * tmp + th2)))
                  * std::complex<double>(std::cos(f02 * tmp), -std::sin(f02 * tmp));
    }
    break;
  }
  }

  time_offset += no_samples;
}

void Rice_Fading_Generator::init_MEDS()
{
  vec n;
  double sgm_0_2;

  switch (dopp_spectrum) {
  case Jakes:
    n = linspace(1, Ni, Ni);
    f1 = sin(pi / (2 * Ni) * (n - 0.5));
    c1 = std::sqrt(1.0 / Ni) * ones(Ni);
    th1 = randu(Ni) * 2 * pi;
    n = linspace(1, Ni + 1, Ni + 1);
    f2 = sin(pi / (2 * (Ni + 1)) * (n - 0.5));
    c2 = std::sqrt(1.0 / (Ni + 1)) * ones(Ni + 1);
    th2 = randu(Ni + 1) * 2 * pi;
    f01 = f02 = 0;
    break;
  case GaussI:
    n = linspace(1, Ni, Ni);
    sgm_0_2 = 5.0 / 6.0;
    c1 = std::sqrt(sgm_0_2 * 2.0 / Ni) * ones(Ni);
    f1 = std::sqrt(2.0) * 0.05 * erfinv((2 * n - 1) / (2 * Ni));
    th1 = randu(Ni) * 2 * pi;
    sgm_0_2 = 1.0 / 6.0;
    c2 = std::sqrt(sgm_0_2 * 2.0 / Ni) * ones(Ni);
    f2 = std::sqrt(2.0) * 0.1 * erfinv((2 * n - 1) / (2 * Ni));
    th2 = randu(Ni) * 2 * pi;
    f01 = 0.8;
    f02 = -0.4;
    break;
  case GaussII:
    n = linspace(1, Ni, Ni);
    sgm_0_2 = std::sqrt(10.0) / (std::sqrt(10.0) + 0.15);
    c1 = std::sqrt(sgm_0_2 * 2.0 / Ni) * ones(Ni);
    f1 = std::sqrt(2.0) * 0.1 * erfinv((2 * n - 1) / (2 * Ni));
    th1 = randu(Ni) * 2 * pi;
    sgm_0_2 = 0.15 / (std::sqrt(10.0) + 0.15);
    c2 = std::sqrt(sgm_0_2 * 2.0 / Ni) * ones(Ni);
    f2 = std::sqrt(2.0) * 0.15 * erfinv((2 * n - 1) / (2 * Ni));
    th2 = randu(Ni) * 2 * pi;
    f01 = -0.7;
    f02 = 0.4;
    break;
  default:
    it_error("Rice_Fading_Generator::init_MEDS(): Wrong spectrum method for this fading generator");
  };
}


// --------------------------------------------------------------------------
// FIR_Fading_Generator class methods
// --------------------------------------------------------------------------

FIR_Fading_Generator::FIR_Fading_Generator(double norm_doppler,
    int filter_length) :
    Correlated_Fading_Generator(norm_doppler)
{
  set_filter_length(filter_length);
}

void FIR_Fading_Generator::set_filter_length(int filter_length)
{
  it_assert(filter_length >= 50,
            "FIR_Fading_Generator::set_filter_length(): Filter length should be at least 50");
  fir_length = filter_length;
  init_flag = false;
}

void FIR_Fading_Generator::init()
{
  // calculate a reasonable upsample rate so that normalized doppler is > 0.1
  double norm_dopp = n_dopp;
  upsample_rate = 1;
  while (norm_dopp < 0.1) {
    norm_dopp *= 2;
    upsample_rate *= 2;
  }
  fir_filter.set_coeffs(Jakes_filter(norm_dopp, fir_length));

  // fill filter with dummy data
  cvec dummy = fir_filter(randn_c(fir_length));

  left_overs.set_size(0, false);

  init_flag = true; // generator ready to use
}

void FIR_Fading_Generator::generate(int no_samples, cvec &output)
{
  if (init_flag == false)
    init();

  // calculate number of samples before upsampling
  int no_upsamples = ceil_i(static_cast<double>(no_samples - left_overs.size())
                            / upsample_rate) + 1;

  // should make a smarter interpolation here!!!
  lininterp(fir_filter(randn_c(no_upsamples)), upsample_rate, output);
  output = concat(left_overs, output); // add left-overs from previous filtering
  left_overs = output.right(output.size() - no_samples); // save left-overs for next round of filtering
  output.set_size(no_samples, true);

  if (los_power > 0.0) { // LOS component exist
    for (int i = 0; i < no_samples; i++) {
      add_LOS(i, output(i));
    }
  }

  time_offset += no_samples;
}

vec FIR_Fading_Generator::Jakes_filter(double norm_dopp, int order)
{
  int L = order / 2;
  vec x_pos(L), x_neg(L), x(2*L + 1), h(2*L + 1);
  for (int i = 1; i <= L; i++) {
    x_pos(i - 1) = besselj(0.25, m_2pi * norm_dopp * i) / std::pow(i, 0.25);
    // / std::sqrt(std::sqrt(static_cast<double>(i)));
  }
  double x0 = 1.468813 * std::pow(norm_dopp, 0.25); // std::sqrt(std::sqrt(norm_dopp));
  x_neg = reverse(x_pos);
  x = concat(concat(x_neg, x0), x_pos);
  h = elem_mult(hamming(2 * L + 1), x);
  h /= norm(h);
  return h;
}


// --------------------------------------------------------------------------
// IFFT_Fading_Generator class methods
// --------------------------------------------------------------------------

void IFFT_Fading_Generator::generate(int no_samples, cvec &output)
{
  if (init_flag == false)
    init();

  generate_Jakes(no_samples, output);

  if (los_power > 0.0) { // LOS component exist
    for (int i = 0; i < no_samples; i++) {
      add_LOS(i, output(i));
    }
  }

  time_offset += no_samples;
}

void IFFT_Fading_Generator::generate_Jakes(int no_samples, cvec &output)
{
  int Nfft = pow2i(levels2bits(no_samples));
  double df = 1.0 / Nfft;
  int noisesamp = ceil_i(n_dopp / df);
  int no_upsample = 1;

  while (noisesamp <= 10) { // if too few samples, increase the FFT size
    Nfft *= 2;
    no_upsample *= 2;
    df = 1.0 / Nfft;
    noisesamp = ceil_i(n_dopp / df);
    it_assert(no_upsample < 128,
              "IFFT_Fading_Generator::generate_Jakes(): Too low normalized doppler or too small blocks of data. Results in an inefficient algorithm with lots of zero-padding");
  }

  vec Fpos = linspace(0, 0.5, Nfft / 2 + 1);
  vec F = concat(Fpos, reverse(-Fpos(1, Nfft / 2 - 1)));
  vec S = zeros(Nfft);

  for (int i = 0; i < F.size(); i++) {
    if (std::fabs(F(i)) < n_dopp)
      S(i) = std::sqrt(1.5 / (pi * n_dopp * std::sqrt(1 - std::pow(F(i) / n_dopp, 2))));
    else if (std::fabs(F(i)) == n_dopp)
      S(i) = 1000000;
  }

  S /= norm(S, 2);
  S *= Nfft;

  cvec x = zeros_c(Nfft);

  for (int i = 0; i < noisesamp; ++i) {
    x(i) = S(i) * randn_c();
    x(Nfft - 1 - i) = S(Nfft - 1 - i) * randn_c();
  }

  x = ifft(x);

  output = x.mid(0, no_samples);
}


// --------------------------------------------------------------------------
// Channel_Specification class methods
// --------------------------------------------------------------------------

Channel_Specification::Channel_Specification(const vec &avg_power_dB,
    const vec &delay_prof)
{
  set_channel_profile(avg_power_dB, delay_prof);
}

Channel_Specification::Channel_Specification(const CHANNEL_PROFILE profile)
{
  set_channel_profile(profile);
}

void Channel_Specification::set_channel_profile(const vec &avg_power_dB, const vec &delay_prof)
{
  it_assert(min(delay_prof) == 0,
            "Channel_Specification::set_channel_profile(): Minimum relative delay must be 0");
  it_assert(avg_power_dB.size() == delay_prof.size(),
            "Channel_Specification::set_channel_profile(): Power and delay vectors must be of equal length");
  it_assert(delay_prof(0) == 0,
            "Channel_Specification::set_channel_profile(): First tap must be at zero delay");
  for (int i = 1; i < delay_prof.size(); i++) {
    it_assert(delay_prof(i) > delay_prof(i - 1),
              "Channel_Specification::set_channel_profile(): Delays should be sorted and unique");
  }

  N_taps = delay_prof.size();
  a_prof_dB = avg_power_dB;
  d_prof = delay_prof;

  // set doppler spectrum to Jakes per default
  tap_doppler_spectrum.set_size(N_taps, false);
  tap_doppler_spectrum = Jakes;

  // set LOS parameters to zeros per default
  set_LOS(zeros(N_taps));
}

void Channel_Specification::set_channel_profile(const CHANNEL_PROFILE profile)
{
  switch (profile) {
    // -------------- ITU Channel models -----------------
  case ITU_Vehicular_A:
    set_channel_profile(vec("0 -1 -9 -10 -15 -20"),
                        vec("0 310 710 1090 1730 2510") * 1e-9);
    break;

  case ITU_Vehicular_B:
    set_channel_profile(vec("-2.5 0 -12.8 -10 -25.2 -16"),
                        vec("0 300 8900 12900 17100 20000") * 1e-9);
    break;

  case ITU_Pedestrian_A:
    set_channel_profile(vec("0 -9.7 -19.2 -22.8"),
                        vec("0 110 190 410") * 1e-9);
    break;

  case ITU_Pedestrian_B:
    set_channel_profile(vec("0 -0.9 -4.9 -8 -7.8 -23.9"),
                        vec("0 200 800 1200 2300 3700") * 1e-9);
    break;

    // -------------- COST259 Channel models -----------------
  case COST259_TUx:
    set_channel_profile(vec("-5.7 -7.6 -10.1 -10.2 -10.2 -11.5 -13.4 -16.3 -16.9 -17.1 -17.4 -19 -19 -19.8 -21.5 -21.6 -22.1 -22.6 -23.5 -24.3"),
                        vec("0 217 512 514 517 674 882 1230 1287 1311 1349 1533 1535 1622 1818 1836 1884 1943 2048 2140") * 1e-9);
    break;

  case COST259_RAx:
    set_channel_profile(vec("-5.2 -6.4 -8.4 -9.3 -10 -13.1 -15.3 -18.5 -20.4 -22.4"),
                        vec("0 42 101 129 149 245 312 410 469 528") * 1e-9);
    set_LOS(0, sqr(0.91 / 0.41), 0.7);
    break;

  case COST259_HTx:
    set_channel_profile(vec("-3.6 -8.9 -10.2 -11.5 -11.8 -12.7 -13.0 -16.2 -17.3 -17.7 -17.6 -22.7 -24.1 -25.8 -25.8 -26.2 -29 -29.9 -30 -30.7"),
                        vec("0 356 441 528 546 609 625 842 916 941 15000 16172 16492 16876 16882 16978 17615 17827 17849 18016") * 1e-9);
    break;

    // -------------- COST207 Channel models -----------------
  case COST207_RA:
    set_channel_profile(vec("0 -2 -10 -20"),
                        vec("0 200 400 600") * 1e-9);
    set_LOS(0, sqr(0.91 / 0.41), 0.7);
    break;

  case COST207_RA6:
    set_channel_profile(vec("0 -4 -8 -12 -16 -20"),
                        vec("0 100 200 300 400 500") * 1e-9);
    set_LOS(0, sqr(0.91 / 0.41), 0.7);
    break;

  case COST207_TU:
    set_channel_profile(vec("-3 0 -2 -6 -8 -10"),
                        vec("0 200 600 1600 2400 5000") * 1e-9);
    set_doppler_spectrum(2, GaussI);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_TU6alt:
    set_channel_profile(vec("-3 0 -2 -6 -8 -10"),
                        vec("0 200 500 1600 2300 5000") * 1e-9);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_TU12:
    set_channel_profile(vec("-4 -3 0 -2 -3 -5 -7 -5 -6 -9 -11 -10"),
                        vec("0 200 400 600 800 1200 1400 1800 2400 3000 3200 5000") * 1e-9);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussI);
    set_doppler_spectrum(6, GaussI);
    set_doppler_spectrum(7, GaussI);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;

  case COST207_TU12alt:
    set_channel_profile(vec("-4 -3 0 -2.6 -3 -5 -7 -5 -6.5 -8.6 -11 -10"),
                        vec("0 200 400 600 800 1200 1400 1800 2400 3000 3200 5000") * 1e-9);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussI);
    set_doppler_spectrum(6, GaussI);
    set_doppler_spectrum(7, GaussI);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;

  case COST207_BU:
    set_channel_profile(vec("-3 0 -3 -5 -2 -4"),
                        vec("0 400 1000 1600 5000 6600") * 1e-9);
    set_doppler_spectrum(2, GaussI);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_BU6alt:
    set_channel_profile(vec("-2.5 0 -3 -5 -2 -4"),
                        vec("0 300 1000 1600 5000 6600") * 1e-9);
    set_doppler_spectrum(2, GaussI);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_BU12:
    set_channel_profile(vec("-7 -3 -1 0 -2 -6 -7 -1 -2 -7 -10 -15"),
                        vec("0 200 400 800 1600 2200 3200 5000 6000 7200 8200 10000") * 1e-9);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussII);
    set_doppler_spectrum(6, GaussII);
    set_doppler_spectrum(7, GaussII);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;

  case COST207_BU12alt:
    set_channel_profile(vec("-7.7 -3.4 -1.3 0 -2.3 -5.6 -7.4 -1.4 -1.6 -6.7 -9.8 -15.1"),
                        vec("0 100 300 700 1600 2200 3100 5000 6000 7200 8100 10000") * 1e-9);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussII);
    set_doppler_spectrum(6, GaussII);
    set_doppler_spectrum(7, GaussII);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;


  case COST207_HT:
    set_channel_profile(vec("0 -2 -4 -7 -6 -12"),
                        vec("0 200 400 600 15000 17200") * 1e-9);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_HT6alt:
    set_channel_profile(vec("0 -1.5 -4.5 -7.5 -8 -17.7"),
                        vec("0 100 300 500 15000 17200") * 1e-9);
    set_doppler_spectrum(4, GaussII);
    set_doppler_spectrum(5, GaussII);
    break;

  case COST207_HT12:
    set_channel_profile(vec("-10 -8 -6 -4 0 0 -4 -8 -9 -10 -12 -14"),
                        vec("0 200 400 600 800 2000 2400 15000 15200 15800 17200 20000") * 1e-9);
    set_doppler_spectrum(3, GaussI);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussI);
    set_doppler_spectrum(6, GaussII);
    set_doppler_spectrum(7, GaussII);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;

  case COST207_HT12alt:
    set_channel_profile(vec("-10 -8 -6 -4 0 0 -4 -8 -9 -10 -12 -14"),
                        vec("0 100 300 500 700 1000 1300 15000 15200 15700 17200 20000") * 1e-9);
    set_doppler_spectrum(4, GaussI);
    set_doppler_spectrum(5, GaussI);
    set_doppler_spectrum(6, GaussI);
    set_doppler_spectrum(7, GaussII);
    set_doppler_spectrum(8, GaussII);
    set_doppler_spectrum(9, GaussII);
    set_doppler_spectrum(10, GaussII);
    set_doppler_spectrum(11, GaussII);
    break;
  };
}


void Channel_Specification::set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum)
{
  for (int i = 0; i < N_taps; i++)
    tap_doppler_spectrum(i) = tap_spectrum[i];
}

void Channel_Specification::set_doppler_spectrum(int tap_number, DOPPLER_SPECTRUM tap_spectrum)
{
  tap_doppler_spectrum(tap_number) = tap_spectrum;
}

void Channel_Specification::set_LOS(int tap_number, double relative_power,
                                    double relative_doppler)
{
  it_assert(N_taps >= 1,
            "Channel_Specification::set_LOS(): Cannot set LOS component if not set channel profile");
  it_assert((tap_number >= 0) && (tap_number < N_taps),
            "Channel_Specification::set_LOS(): Tap number out of range");
  it_assert((relative_doppler >= 0) && (relative_doppler <= 1.0),
            "Channel_Specification::set_LOS(): Normalized Doppler out of range");
  it_assert(relative_power >= 0.0,
            "Channel_Specification::set_LOS(): Rice factor out of range");

  los_power.set_size(N_taps, true);
  los_dopp.set_size(N_taps, true);
  los_power(tap_number) = relative_power;
  los_dopp(tap_number) = relative_doppler;
}

void Channel_Specification::set_LOS(const vec& relative_power,
                                    const vec& relative_doppler)
{
  it_assert((relative_power.size() == N_taps),
            "Channel_Specification::set_LOS(): Improper size of input vectors");

  if (relative_doppler.size() == 0) {
    los_power.set_size(relative_power.size());
    los_dopp.set_size(relative_power.size());
    for (int i = 0; i < relative_power.size(); i++) {
      it_assert(relative_power(i) >= 0.0,
                "Channel_Specification::set_LOS(): Rice factor out of range");
      los_power(i) = relative_power(i);
      los_dopp(i) = 0.7;
    }
  }
  else {
    it_assert(relative_doppler.size() == N_taps,
              "Channel_Specification::set_LOS(): Improper size of input vectors");
    los_power.set_size(relative_power.size());
    los_dopp.set_size(relative_power.size());
    for (int i = 0; i < relative_power.size(); i++) {
      it_assert((relative_doppler(i) >= 0) && (relative_doppler(i) <= 1.0),
                "Channel_Specification::set_LOS(): Normalized Doppler out of range");
      it_assert(relative_power(i) >= 0.0,
                "Channel_Specification::set_LOS(): Rice factor out of range");
      los_power(i) = relative_power(i);
      los_dopp(i) = relative_doppler(i);
    }
  }
}

void Channel_Specification::get_channel_profile(vec &avg_power_dB,
    vec &delay_prof) const
{
  avg_power_dB = a_prof_dB;
  delay_prof = d_prof;
}

DOPPLER_SPECTRUM Channel_Specification::get_doppler_spectrum(int index) const
{
  it_assert((index >= 0) && (index < N_taps),
            "Channel_Specification::get_doppler_spectrum(): Index of of range");
  return tap_doppler_spectrum(index);
}

double Channel_Specification::calc_mean_excess_delay() const
{
  vec a_prof = inv_dB(a_prof_dB);
  return (a_prof * d_prof / sum(a_prof));
}

double Channel_Specification::calc_rms_delay_spread() const
{
  vec a_prof = inv_dB(a_prof_dB);
  double a = a_prof * d_prof / sum(a_prof);
  double b = a_prof * sqr(d_prof) / sum(a_prof);

  return std::sqrt(b - a*a);
}


// --------------------------------------------------------------------------
// TDL_Channel class methods
// --------------------------------------------------------------------------

TDL_Channel::TDL_Channel(const vec &avg_power_dB, const ivec &delay_prof):
    init_flag(false), n_dopp(0.0), fading_type(Independent), method(Rice_MEDS),
    filter_length(0), nrof_freq(16), discrete_Ts(0.0)
{
  set_channel_profile(avg_power_dB, delay_prof);

  // initialize LOS parameters to all zeros
  set_LOS(zeros(delay_prof.size()));

  // initialize Doppler spectra
  tap_doppler_spectrum.set_size(delay_prof.size());
  tap_doppler_spectrum = Jakes;
}

TDL_Channel::TDL_Channel(const Channel_Specification &channel_spec, double sampling_time):
    init_flag(false), n_dopp(0.0), fading_type(Independent), method(Rice_MEDS),
    filter_length(0), nrof_freq(16), discrete_Ts(sampling_time)
{
  set_channel_profile(channel_spec, sampling_time);

  // set Doppler spectrum
  tap_doppler_spectrum = channel_spec.get_doppler_spectrum();
}

TDL_Channel::~TDL_Channel()
{
  if (fading_gen.size() > 0) { // delete all old generators
    for (int i = 0; i < fading_gen.size(); i++) {
      if (fading_gen(i) != NULL) {
        delete fading_gen(i);
        fading_gen(i) = NULL;
      }
    }
  }
}

void TDL_Channel::set_channel_profile(const vec &avg_power_dB,
                                      const ivec &delay_prof)
{
  it_assert(min(delay_prof) == 0,
            "TDL_Channel::set_channel_profile(): Minimum relative delay must be 0.");
  it_assert(avg_power_dB.size() == delay_prof.size(),
            "TDL_Channel::set_channel_profile(): Power and delay vectors must be of equal length!");
  it_assert(delay_prof(0) == 0,
            "TDL_Channel::set_channel_profile(): First tap must be at zero delay");
  for (int i = 1; i < delay_prof.size(); i++) {
    it_assert(delay_prof(i) > delay_prof(i - 1),
              "TDL_Channel::set_channel_profile(): Delays should be sorted and unique");
  }

  N_taps = delay_prof.size();
  a_prof = pow(10.0, avg_power_dB / 20.0); // Convert power profile to amplitude profile
  a_prof /= norm(a_prof); // Normalize amplitude profile
  d_prof = delay_prof;

  // initialize Doppler spectra
  tap_doppler_spectrum.set_size(d_prof.size());
  tap_doppler_spectrum = Jakes;

  // set size of Rice parameters according to the new channel profile
  set_LOS(zeros(N_taps));

  // changes in PDP require initialisation
  init_flag = false;
}

void TDL_Channel::set_channel_profile_uniform(int no_taps)
{
  it_assert(no_taps >= 1, "TDL_Channel::set_channel_profile_uniform(): Minimum number of taps is 1.");

  vec avg_power_dB = zeros(no_taps);
  ivec delay_prof(no_taps);
  for (int i = 0; i < no_taps; i++)
    delay_prof(i) = i;

  set_channel_profile(avg_power_dB, delay_prof);
}

void TDL_Channel::set_channel_profile_exponential(int no_taps)
{
  it_assert(no_taps >= 1, "TDL_Channel::set_channel_profile_exponential(): Minimum number of taps is 1.");

  vec avg_power_dB(no_taps);
  ivec delay_prof(no_taps);
  for (int i = 0; i < no_taps; i++) {
    delay_prof(i) = i;
    // p(i*ts) = exp(-i*ts),    k = 0...no_taps-1
    avg_power_dB(i) = dB(std::exp(static_cast<double>(-i)));
  }

  set_channel_profile(avg_power_dB, delay_prof);
}

void TDL_Channel::set_channel_profile(const Channel_Specification &channel_spec, double sampling_time)
{
  vec avg_power_dB;
  vec delay_profile;

  // set power profile and delays
  channel_spec.get_channel_profile(avg_power_dB, delay_profile);
  discrete_Ts = sampling_time;
  N_taps = avg_power_dB.size();
  a_prof = pow(10.0, avg_power_dB / 20.0); // Convert power profile to amplitude profile
  a_prof /= norm(a_prof); // Normalize amplitude profile

  // set size of Rice parameters according to the new channel profile
  set_LOS(channel_spec.get_LOS_power(), channel_spec.get_LOS_doppler());

  // set Doppler spectrum
  tap_doppler_spectrum = channel_spec.get_doppler_spectrum();

  // sets discretized delay profile
  discretize(delay_profile);

  init_flag = false;
}


void TDL_Channel::set_correlated_method(CORRELATED_METHOD correlated_method)
{
  fading_type = Correlated;
  method = correlated_method;
  init_flag = false;
}

void TDL_Channel::set_fading_type(FADING_TYPE fading_type_in)
{
  fading_type = fading_type_in;
  init_flag = false;
}


void TDL_Channel::set_norm_doppler(double norm_doppler)
{
  it_assert((norm_doppler > 0) && (norm_doppler <= 1.0),
            "TDL_Channel::set_norm_doppler(): Normalized Doppler out of range");
  n_dopp = norm_doppler;
  // if normalized Doppler is set, we have correlated fading
  fading_type = Correlated;
  init_flag = false;
}


void TDL_Channel::set_LOS(const vec& relative_power, const vec& relative_doppler)
{
  it_assert((relative_power.size() == N_taps),
            "TDL_Channel::set_LOS(): Improper size of input vectors");

  if (relative_doppler.size() == 0) {
    los_power.set_size(relative_power.size());
    los_dopp.set_size(relative_power.size());
    for (int i = 0; i < relative_power.size(); i++) {
      it_assert(relative_power(i) >= 0.0,
                "TDL_Channel::set_LOS(): Rice factor out of range");
      los_power(i) = relative_power(i);
      los_dopp(i) = (relative_power(i) > 0) ? 0.7 : 0.0;
    }
  }
  else {
    it_assert(relative_doppler.size() == N_taps,
              "TDL_Channel::set_LOS(): Improper size of input vectors");
    los_power.set_size(relative_power.size());
    los_dopp.set_size(relative_power.size());
    for (int i = 0; i < relative_power.size(); i++) {
      it_assert((relative_doppler(i) >= 0) && (relative_doppler(i) <= 1.0),
                "TDL_Channel::set_LOS(): Normalized Doppler out of range");
      it_assert(relative_power(i) >= 0.0,
                "TDL_Channel::set_LOS(): Rice factor out of range");
      los_power(i) = relative_power(i);
      los_dopp(i) = relative_doppler(i);
    }
  }
}

void TDL_Channel::set_LOS_power(const vec& relative_power)
{
  it_assert(relative_power.size() == N_taps,
            "TDL_Channel::set_LOS_power(): Improper size of input vector");

  los_power.set_size(relative_power.size());
  los_dopp.set_size(relative_power.size());
  for (int i = 0; i < los_power.size(); ++i) {
    los_power(i) = relative_power(i);
    los_dopp(i) = (relative_power(i) > 0) ? 0.7 : 0.0;
  }
  init_flag = false;
}

void TDL_Channel::set_LOS_doppler(const vec& relative_doppler)
{
  it_assert(relative_doppler.size() == los_power.size(),
            "TDL_Channel::set_LOS_doppler(): Improper size of input vector");

  it_assert(n_dopp > 0, "TDL_Channel::set_LOS_doppler(): Normalized Doppler needs to be non zero to set the LOS Doppler in a Correlated fading generator");

  los_dopp.set_size(relative_doppler.size());
  for (int i = 0; i < relative_doppler.size(); ++i) {
    it_assert((relative_doppler(i) >= 0) && (relative_doppler(i) <= 1.0),
              "TDL_Channel::set_LOS_doppler(): Normalized Doppler out of range");
    los_dopp(i) = relative_doppler(i);
  }

  init_flag = false;
}


void TDL_Channel::set_doppler_spectrum(const DOPPLER_SPECTRUM *tap_spectrum)
{
  it_assert(N_taps > 0, "TDL_Channel::set_doppler_spectrum(): Channel profile not defined yet");

  it_assert(n_dopp > 0, "TDL_Channel::set_doppler_spectrum(): Normalized Doppler needs to be non zero to set the Doppler spectrum in the Correlated Rice MEDS fading generator");

  if (method != Rice_MEDS)
    method = Rice_MEDS;

  tap_doppler_spectrum.set_size(N_taps, false);
  for (int i = 0; i < N_taps; i++)
    tap_doppler_spectrum(i) = tap_spectrum[i];

  init_flag = false;
}

void TDL_Channel::set_doppler_spectrum(int tap_number, DOPPLER_SPECTRUM tap_spectrum)
{
  it_assert((tap_number >= 0) && (tap_number < N_taps),
            "TDL_Channel::set_doppler_spectrum(): Improper tap number");

  it_assert(n_dopp > 0, "TDL_Channel::set_doppler_spectrum(): Normalized Doppler needs to be non zero to set the Doppler spectrum in the Correlated Rice MEDS fading generator");

  if (method != Rice_MEDS)
    method = Rice_MEDS;

  tap_doppler_spectrum.set_size(N_taps, true);
  tap_doppler_spectrum(tap_number) = tap_spectrum;

  init_flag = false;
}

void TDL_Channel::set_no_frequencies(int no_freq)
{
  it_assert(n_dopp > 0, "TDL_Channel::set_no_frequencies(): Normalized Doppler needs to be non zero to set the number of frequencies in the Correlated Rice MEDS fading generator");
  nrof_freq = no_freq;
  if (method != Rice_MEDS)
    method = Rice_MEDS;

  init_flag = false;
}


void TDL_Channel::set_filter_length(int fir_length)
{
  it_assert(n_dopp > 0, "TDL_Channel::set_filter_length(): Normalized Doppler needs to be non zero to use the Correlated FIR fading generator");

  filter_length = fir_length;
  if (method != FIR)
    method = FIR;

  init_flag = false;
}


void TDL_Channel::set_time_offset(int offset)
{
  it_assert(n_dopp > 0, "TDL_Channel::set_time_offset(): Normalized Doppler needs to be non zero to set time offset in a Correlated fading generator");

  if (init_flag == false)
    init();

  for (int i = 0; i < N_taps; i++) {
    fading_gen(i)->set_time_offset(offset);
  }
}


void TDL_Channel::shift_time_offset(int no_samples)
{
  it_assert(n_dopp > 0, "TDL_Channel::shift_time_offset(): Normalized Doppler needs to be non zero to shift time offset in a Correlated fading generator");

  if (init_flag == false)
    init();

  for (int i = 0; i < N_taps; i++) {
    fading_gen(i)->shift_time_offset(no_samples);
  }
}


void TDL_Channel::get_channel_profile(vec &avg_power_dB,
                                      ivec &delay_prof) const
{
  avg_power_dB = 20 * log10(a_prof);
  delay_prof = d_prof;
}

vec TDL_Channel::get_avg_power_dB() const
{
  return (20 * log10(a_prof));
}

double TDL_Channel::get_time_offset() const
{
  if (fading_gen(0) != NULL)
    return fading_gen(0)->get_time_offset();
  else
    return -1.0;
}

double TDL_Channel::calc_mean_excess_delay() const
{
  return (sqr(a_prof)*d_prof / sum_sqr(a_prof));
}

double TDL_Channel::calc_rms_delay_spread() const
{
  double a = (sqr(a_prof) * d_prof / sum_sqr(a_prof));
  double b = (sqr(a_prof) * sqr(to_vec(d_prof)) / sum_sqr(a_prof));

  return (std::sqrt(b - a*a));
}

void TDL_Channel::init()
{
  it_assert(N_taps > 0, "TDL_Channel::init(): Channel profile not defined yet");
  it_assert(N_taps == los_power.size(),
            "TDL_Channel::init(): LOS profile does not mach the channel profile");

  if (fading_gen.size() > 0) { // delete all old generators
    for (int i = 0; i < fading_gen.size(); i++) {
      if (fading_gen(i) != NULL) {
        delete fading_gen(i);
        fading_gen(i) = NULL;
      }
    }
  }

  // create all generators and set the parameters
  fading_gen.set_size(N_taps, false);
  switch (fading_type) {

  case Independent:
    for (int i = 0; i < N_taps; ++i) {
      fading_gen(i) = new Independent_Fading_Generator();
      if (los_power(i) > 0)
        fading_gen(i)->set_LOS_power(los_power(i));
      fading_gen(i)->init();
    }
    break;

  case Static:
    for (int i = 0; i < N_taps; ++i) {
      fading_gen(i) = new Static_Fading_Generator();
      if (los_power(i) > 0)
        fading_gen(i)->set_LOS_power(los_power(i));
      fading_gen(i)->init();
    }
    break;

  case Correlated:
    it_assert(n_dopp > 0,
              "TDL_Channel::init(): Correlated fading requires non zero normalized Doppler");

    switch (method) {
    case Rice_MEDS:
      // The third parameter is the number of sine waveforms that create the
      // fading process. It is increased by 2 for each tap to make taps
      // uncorrelated. Minimum number of waveforms set explicitly to 16.
      for (int i = 0; i < N_taps; ++i) {
        fading_gen(i) = new Rice_Fading_Generator(n_dopp, tap_doppler_spectrum(i),
            nrof_freq + 2*i, MEDS);
        if (los_power(i) > 0) {
          fading_gen(i)->set_LOS_power(los_power(i));
          fading_gen(i)->set_LOS_doppler(los_dopp(i));
        }
        fading_gen(i)->init();
      }
      break;

    case FIR:
      for (int i = 0; i < N_taps; ++i) {
        it_assert(tap_doppler_spectrum(i) == Jakes,
                  "TDL_Channel::init(): FIR fading generator can be used with Jakes spectrum only");
        fading_gen(i) = new FIR_Fading_Generator(n_dopp);
        if (los_power(i) > 0) {
          fading_gen(i)->set_LOS_power(los_power(i));
          fading_gen(i)->set_LOS_doppler(los_dopp(i));
        }
        if (filter_length > 0)
          fading_gen(i)->set_filter_length(filter_length);
        fading_gen(i)->init();
      }
      break;

    case IFFT:
      for (int i = 0; i < N_taps; ++i) {
        it_assert(tap_doppler_spectrum(i) == Jakes,
                  "TDL_Channel::init(): IFFT fading generator can be used with Jakes spectrum only");
        fading_gen(i) = new IFFT_Fading_Generator(n_dopp);
        if (los_power(i) > 0) {
          fading_gen(i)->set_LOS_power(los_power(i));
          fading_gen(i)->set_LOS_doppler(los_dopp(i));
        }
        fading_gen(i)->init();
      }
      break;

    default:
      it_error("TDL_Channel::init(): No such fading generation method");
    }
    break;

  default:
    it_error("TDL_Channel::init(): No such fading type");
  };

  init_flag = true;
}

void TDL_Channel::generate(int no_samples, Array<cvec> &channel_coeff)
{
  if (init_flag == false)
    init();

  channel_coeff.set_size(N_taps, false);
  for (int i = 0; i < N_taps; i++)
    channel_coeff(i) = a_prof(i) * fading_gen(i)->generate(no_samples);
}

void TDL_Channel::generate(int no_samples, cmat &channel_coeff)
{
  if (init_flag == false)
    init();

  channel_coeff.set_size(no_samples, N_taps, false);
  for (int i = 0; i < N_taps; i++)
    channel_coeff.set_col(i, a_prof(i) * fading_gen(i)->generate(no_samples));
}


void TDL_Channel::filter_known_channel(const cvec &input, cvec &output, const Array<cvec> &channel_coeff)
{
  int maxdelay = max(d_prof);

  output.set_size(input.size() + maxdelay, false);
  output.zeros();

  for (int i = 0; i < N_taps; i++)
    output += concat(zeros_c(d_prof(i)), elem_mult(input, channel_coeff(i)), zeros_c(maxdelay - d_prof(i)));
}

void TDL_Channel::filter_known_channel(const cvec &input, cvec &output, const cmat &channel_coeff)
{
  int maxdelay = max(d_prof);

  output.set_size(input.size() + maxdelay, false);
  output.zeros();

  for (int i = 0; i < N_taps; i++)
    output += concat(zeros_c(d_prof(i)), elem_mult(input, channel_coeff.get_col(i)), zeros_c(maxdelay - d_prof(i)));
}

void TDL_Channel::filter(const cvec &input, cvec &output, Array<cvec> &channel_coeff)
{
  generate(input.size(), channel_coeff);
  filter_known_channel(input, output, channel_coeff);
}

void TDL_Channel::filter(const cvec &input, cvec &output, cmat &channel_coeff)
{
  generate(input.size(), channel_coeff);
  filter_known_channel(input, output, channel_coeff);
}

cvec TDL_Channel::filter(const cvec &input, Array<cvec> &channel_coeff)
{
  cvec output;
  filter(input, output, channel_coeff);
  return output;
}

cvec TDL_Channel::filter(const cvec &input, cmat &channel_coeff)
{
  cvec output;
  filter(input, output, channel_coeff);
  return output;
}

void TDL_Channel::filter(const cvec &input, cvec &output)
{
  Array<cvec> channel_coeff;
  filter(input, output, channel_coeff);
}

cvec TDL_Channel::filter(const cvec &input)
{
  cvec output;
  filter(input, output);
  return output;
}


void TDL_Channel::operator()(const cvec &input, cvec &output, Array<cvec> &channel_coeff)
{
  filter(input, output, channel_coeff);
}

void TDL_Channel::operator()(const cvec &input, cvec &output, cmat &channel_coeff)
{
  filter(input, output, channel_coeff);
}


cvec TDL_Channel::operator()(const cvec &input, Array<cvec> &channel_coeff)
{
  return filter(input, channel_coeff);
}

cvec TDL_Channel::operator()(const cvec &input, cmat &channel_coeff)
{
  return filter(input, channel_coeff);
}

cvec TDL_Channel::operator()(const cvec &input)
{
  return filter(input);
}


void TDL_Channel::calc_impulse_response(const Array<cvec> &channel_coeff, Array<cvec> &impulse_response)
{
  it_assert(init_flag == true, "calc_impulse_response: TDL_Channel is not initialized");
  it_assert(N_taps == channel_coeff.size(), "calc_impulse_response: number of channel taps do not match");

  int no_samples = channel_coeff(0).size();
  it_assert(no_samples > 0, "calc_impulse_response: channel_coeff must contain samples");

  impulse_response.set_size(no_samples);

  for (int i = 0; i < no_samples; i++) {
    impulse_response(i).set_size(d_prof(N_taps - 1) + 1, false);
    impulse_response(i).zeros();

    for (int j = 0; j < N_taps; j++)
      impulse_response(i)(d_prof(j)) = channel_coeff(j)(i);

  }
}

void TDL_Channel::calc_frequency_response(const Array<cvec> &channel_coeff, Array<cvec> &frequency_response, const int fft_size)
{
  it_assert(init_flag == true, "calc_frequency_response: TDL_Channel is not initialized");
  it_assert(N_taps == channel_coeff.size(), "calc_frequency_response: number of channel taps do not match");

  int no_samples = channel_coeff(0).size();
  it_assert(no_samples > 0, "calc_frequency_response: channel_coeff must contain samples");

  frequency_response.set_size(no_samples);

  it_assert(fft_size > d_prof(N_taps - 1), "calc_frequency_response: fft_size must be larger than the maximum delay in samples");
  cvec impulse_response(fft_size);

  for (int i = 0; i < no_samples; i++) {
    impulse_response.zeros();

    for (int j = 0; j < N_taps; j++)
      impulse_response(d_prof(j)) = channel_coeff(j)(i);

    fft(impulse_response, frequency_response(i));

  }
}

void TDL_Channel::calc_frequency_response(const cmat &channel_coeff, cmat &frequency_response, const int fft_size)
{
  it_assert(init_flag == true, "calc_frequency_response: TDL_Channel is not initialized");
  it_assert(N_taps == channel_coeff.cols(), "calc_frequency_response: number of channel taps do not match");

  int no_samples = channel_coeff.rows();
  it_assert(no_samples > 0, "calc_frequency_response: channel_coeff must contain samples");

  frequency_response.set_size(fft_size, no_samples, false);

  it_assert(fft_size > d_prof(N_taps - 1), "calc_frequency_response: fft_size must be larger than the maximum delay in samples");
  cvec impulse_response(fft_size);
  cvec freq;

  for (int i = 0; i < no_samples; i++) {
    impulse_response.zeros();

    for (int j = 0; j < N_taps; j++)
      impulse_response(d_prof(j)) = channel_coeff(i, j);

    fft(impulse_response, freq);
    frequency_response.set_col(i, freq);
  }
}

void TDL_Channel::discretize(const vec &delay_profile)
{
  it_assert(N_taps > 0, "TDL_Channel::discretize(): No channel profile specified");
  it_assert(delay_profile(0) == 0, "TDL_Channel::discretize(): First tap should be at zero delay");
  it_assert(discrete_Ts > 0, "TDL_Channel::discretize(): Incorrect sampling time");
  it_assert((a_prof.size() == N_taps) && (delay_profile.size() == N_taps)
            && (los_power.size() == N_taps) && (tap_doppler_spectrum.size() == N_taps),
            "TDL_Channel:: discretize(): Channel profile lenghts must be equal to the number of taps");

  vec p_prof = sqr(a_prof); // Power profile
  ivec delay_prof(N_taps);
  vec power(N_taps);
  double spower;
  vec scattered(N_taps), direct(N_taps);
  vec los_doppler(N_taps);
  Array <DOPPLER_SPECTRUM> tap_spectrum(N_taps);

  delay_prof(0) = round_i(delay_profile(0) / discrete_Ts);  // should be at zero delay anyway
  power(0) = p_prof(0);
  spower = p_prof(0) / (1 + los_power(0));
  scattered(0) = spower;
  direct(0) = los_power(0) * spower;
  los_doppler(0) = los_dopp(0);
  tap_spectrum(0) = tap_doppler_spectrum(0);

  // taps within ((j-0.5)Ts,(j+0.5)Ts] are included in the j-th tap
  int j = 0, j_delay = 0;
  for (int i = 1; i < N_taps; i++) {
    if (delay_profile(i) > (j_delay + 0.5)*discrete_Ts) {
      // first skip empty taps
      while (delay_profile(i) > (j_delay + 0.5)*discrete_Ts) { j_delay++; }
      // create a new tap at (j+1)Ts
      j++;
      delay_prof(j) = j_delay;
      power(j) = p_prof(i);
      spower = p_prof(i) / (1 + los_power(i));
      scattered(j) = spower;
      direct(j) = los_power(i) * spower;
      los_doppler(j) = los_dopp(i);
      tap_spectrum(j) = tap_doppler_spectrum(i);
    }
    else {
      // add to the previously created tap
      power(j) += p_prof(i);
      spower = p_prof(i) / (1 + los_power(i));
      scattered(j) += spower;
      direct(j) += los_power(i) * spower;
      it_assert(tap_spectrum(j) == tap_doppler_spectrum(i),
                "TDL_Channel::discretize(): Sampling frequency too low. Can not discretize the channel with different Doppler spectra on merged taps.");
      it_warning("TDL_Channel::discretize(): Sampling frequency too low. Merging original tap " << i << " with new tap " << j << ".");
      if (los_doppler(j) != los_dopp(i)) {
        los_doppler(j) = 0.7;
        it_warning("TDL_Channel::discretize(): LOS Doppler value reset to 0.7 for tap " << j << " due to the merging process.");
      }
    }
  }

  int no_taps = j + 1; // number of taps found
  if (no_taps < N_taps) {
    delay_prof.set_size(no_taps, true);
    power.set_size(no_taps, true);
    direct.set_size(no_taps, true);
    scattered.set_size(no_taps, true);
    los_doppler.set_size(no_taps, true);
    tap_spectrum.set_size(no_taps, true);

    // write over the existing channel profile with its new version
    N_taps = no_taps;
    a_prof = sqrt(power);
    los_power = elem_div(direct, scattered);
    los_dopp = los_doppler;
    tap_doppler_spectrum = tap_spectrum;
  }
  // new discretized path's delays
  d_prof = delay_prof; // in samples
}


// --------------------------------------------------------------------------
// Binary Symetric Channel class methods
// --------------------------------------------------------------------------

bvec BSC::operator()(const bvec &input)
{
  int i, length = input.length();
  bvec output(length);

  for (i = 0; i < length; i++) {
    if (u() <= p) {
      output(i) = input(i) + bin(1);
    }
    else {
      output(i) = input(i);
    }
  }
  return output;
}


// --------------------------------------------------------------------------
// AWGN_Channel class methods
// --------------------------------------------------------------------------

cvec AWGN_Channel::operator()(const cvec &input)
{
  int n = input.size();
  cvec noise(n);
  rng_cn.sample_vector(n, noise);
  noise *= sigma;
  noise += input;
  return noise;
}

vec AWGN_Channel::operator()(const vec &input)
{
  int n = input.size();
  vec noise(n);
  rng_n.sample_vector(n, noise);
  noise *= sigma;
  noise += input;
  return noise;
}


} // namespace itpp
