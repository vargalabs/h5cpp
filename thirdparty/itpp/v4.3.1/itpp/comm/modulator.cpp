/*!
 * \file
 * \brief One- and two-dimensional modulators - source file
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

#include <itpp/comm/modulator.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/specmat.h>


namespace itpp
{

//! \cond
//MSVC explicitely instantiate required template while building the shared library
template class ITPP_EXPORT Modulator<double>;
template class ITPP_EXPORT Modulator<std::complex<double> >;
//! \endcond

// ----------------------------------------------------------------------
// QAM
// ----------------------------------------------------------------------

void QAM::set_M(int Mary)
{
  k = levels2bits(Mary);
  M = Mary;
  it_assert((pow2i(k) == M) && (is_even(k)),
            "QAM::set_M(): M = " << M << " is not an even power of 2");
  L = round_i(std::sqrt(static_cast<double>(M)));

  double average_energy = (M - 1) * 2.0 / 3.0;
  scaling_factor = std::sqrt(average_energy);

  symbols.set_size(M);
  bitmap.set_size(M, k);
  bits2symbols.set_size(M);

  bmat gray_code = graycode(levels2bits(L));

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      symbols(i*L + j) = std::complex<double>(((L - 1) - j * 2) / scaling_factor,
                                              ((L - 1) - i * 2) / scaling_factor);
      bitmap.set_row(i*L + j, concat(gray_code.get_row(i),
                                     gray_code.get_row(j)));
      bits2symbols(bin2dec(bitmap.get_row(i*L + j))) = i * L + j;
    }
  }

  calculate_softbit_matrices();

  setup_done = true;
}


void QAM::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "QAM::demodulate_bits(): Modulator not ready.");
  out.set_size(k*signal.size(), false);

  int temp_real, temp_imag;

  for (int i = 0; i < signal.size(); i++) {
    temp_real = round_i((L - 1) - (std::real(signal(i) * scaling_factor)
                                   + (L - 1)) / 2.0);
    temp_imag = round_i((L - 1) - (std::imag(signal(i) * scaling_factor)
                                   + (L - 1)) / 2.0);
    if (temp_real < 0)
      temp_real = 0;
    else if (temp_real > (L - 1))
      temp_real = (L - 1);
    if (temp_imag < 0)
      temp_imag = 0;
    else if (temp_imag > (L - 1))
      temp_imag = (L - 1);
    out.replace_mid(k*i, bitmap.get_row(temp_imag * L + temp_real));
  }
}

bvec QAM::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}


// ----------------------------------------------------------------------
// PSK
// ----------------------------------------------------------------------

void PSK::set_M(int Mary)
{
  k = levels2bits(Mary);
  M = Mary;
  it_assert(pow2i(k) == M, "PSK::set_M(): M is not a power of 2");

  symbols.set_size(M);
  bitmap = graycode(k);
  bits2symbols.set_size(M);

  double delta = m_2pi / M;
  double epsilon = delta / 10000.0;
  std::complex<double> symb;
  for (int i = 0; i < M; i++) {
    symb = std::complex<double>(std::polar(1.0, delta * i));
    if (std::fabs(std::real(symb)) < epsilon) {
      symbols(i) = std::complex<double>(0.0, std::imag(symb));
    }
    else if (std::fabs(std::imag(symb)) < epsilon) {
      symbols(i) = std::complex<double>(std::real(symb), 0.0);
    }
    else {
      symbols(i) = symb;
    }

    bits2symbols(bin2dec(bitmap.get_row(i))) = i;
  }

  calculate_softbit_matrices();

  setup_done = true;
}


void PSK::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "PSK::demodulate_bits(): Modulator not ready.");
  int est_symbol;
  double ang, temp;

  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
    ang = std::arg(signal(i));
    temp = (ang < 0) ? (m_2pi + ang) : ang;
    est_symbol = round_i(temp * (M >> 1) / pi) % M;
    out.replace_mid(i*k, bitmap.get_row(est_symbol));
  }
}

bvec PSK::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}


// ----------------------------------------------------------------------
// QPSK
// ----------------------------------------------------------------------

void QPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                vec &soft_bits, Soft_Method) const
{
  soft_bits.set_size(k * rx_symbols.size());
  std::complex<double> temp;
  double factor = 2 * std::sqrt(2.0) / N0;
  std::complex<double> exp_pi4 = std::complex<double>(std::cos(pi / 4),
                                 std::sin(pi / 4));
  for (int i = 0; i < rx_symbols.size(); i++) {
    temp = rx_symbols(i) * exp_pi4;
    soft_bits((i << 1) + 1) = std::real(temp) * factor;
    soft_bits(i << 1) = std::imag(temp) * factor;
  }
}

vec QPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                               Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}


void QPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
                                double N0, vec &soft_bits,
                                Soft_Method) const
{
  soft_bits.set_size(2*rx_symbols.size(), false);
  std::complex<double> temp;
  double factor = 2 * std::sqrt(2.0) / N0;
  std::complex<double> exp_pi4 = std::complex<double>(std::cos(pi / 4),
                                 std::sin(pi / 4));
  for (int i = 0; i < rx_symbols.size(); i++) {
    temp = rx_symbols(i) * std::conj(channel(i)) * exp_pi4;
    soft_bits((i << 1) + 1) = std::real(temp) * factor;
    soft_bits(i << 1) = std::imag(temp) * factor;
  }
}

vec QPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
                               double N0, Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, channel, N0, out, method);
  return out;
}


// ----------------------------------------------------------------------
// BPSK_c
// ----------------------------------------------------------------------

void BPSK_c::modulate_bits(const bvec &bits, cvec &out) const
{
  out.set_size(bits.size(), false);
  for (int i = 0; i < bits.size(); i++) {
    out(i) = (bits(i) == 0 ? 1.0 : -1.0);
  }
}

cvec BPSK_c::modulate_bits(const bvec &bits) const
{
  cvec out(bits.size());
  modulate_bits(bits, out);
  return out;
}


void BPSK_c::demodulate_bits(const cvec &signal, bvec &out) const
{
  out.set_size(signal.size(), false);
  for (int i = 0; i < signal.length(); i++) {
    out(i) = (std::real(signal(i)) > 0) ? bin(0) : bin(1);
  }
}

bvec BPSK_c::demodulate_bits(const cvec &signal) const
{
  bvec out(signal.size());
  demodulate_bits(signal, out);
  return out;
}


void BPSK_c::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                  vec &soft_bits, Soft_Method) const
{
  double factor = 4 / N0;
  soft_bits.set_size(rx_symbols.size(), false);

  for (int i = 0; i < rx_symbols.size(); i++) {
    soft_bits(i) = factor * std::real(rx_symbols(i));
  }
}

vec BPSK_c::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                 Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}


void BPSK_c::demodulate_soft_bits(const cvec &rx_symbols,
                                  const cvec &channel,
                                  double N0, vec &soft_bits,
                                  Soft_Method) const
{
  double factor = 4 / N0;
  soft_bits.set_size(rx_symbols.size(), false);

  for (int i = 0; i < rx_symbols.size(); i++) {
    soft_bits(i) = factor * std::real(rx_symbols(i) * std::conj(channel(i)));
  }
}

vec BPSK_c::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
                                 double N0, Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, channel, N0, out, method);
  return out;
}


// ----------------------------------------------------------------------
// BPSK
// ----------------------------------------------------------------------

void BPSK::modulate_bits(const bvec &bits, vec &out) const
{
  out.set_size(bits.size(), false);
  for (int i = 0; i < bits.size(); i++) {
    out(i) = (bits(i) == 0 ? 1.0 : -1.0);
  }
}

vec BPSK::modulate_bits(const bvec &bits) const
{
  vec out(bits.size());
  modulate_bits(bits, out);
  return out;
}


void BPSK::demodulate_bits(const vec &signal, bvec &out) const
{
  out.set_size(signal.size(), false);
  for (int i = 0; i < signal.length(); i++) {
    out(i) = (signal(i) > 0) ? bin(0) : bin(1);
  }
}

bvec BPSK::demodulate_bits(const vec &signal) const
{
  bvec out(signal.size());
  demodulate_bits(signal, out);
  return out;
}


void BPSK::demodulate_soft_bits(const vec &rx_symbols, double N0,
                                vec &soft_bits, Soft_Method) const
{
  double factor = 4 / N0;
  soft_bits.set_size(rx_symbols.size(), false);

  for (int i = 0; i < rx_symbols.size(); i++) {
    soft_bits(i) = factor * rx_symbols(i);
  }
}

vec BPSK::demodulate_soft_bits(const vec &rx_symbols, double N0,
                               Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}


void BPSK::demodulate_soft_bits(const vec &rx_symbols, const vec &channel,
                                double N0, vec &soft_bits,
                                Soft_Method) const
{
  double factor = 4 / N0;
  soft_bits.set_size(rx_symbols.size(), false);

  for (int i = 0; i < rx_symbols.size(); i++) {
    soft_bits(i) = factor * (rx_symbols(i) * channel(i));
  }
}

vec BPSK::demodulate_soft_bits(const vec &rx_symbols, const vec &channel,
                               double N0, Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, channel, N0, out, method);
  return out;
}


// ----------------------------------------------------------------------
// PAM_c
// ----------------------------------------------------------------------

void PAM_c::set_M(int Mary)
{
  M = Mary;
  k = levels2bits(M);
  it_assert(pow2i(k) == M, "PAM_c::set_M(): M is not a power of 2");

  symbols.set_size(M, false);
  bits2symbols.set_size(M, false);
  bitmap = graycode(k);
  double average_energy = (sqr(M) - 1) / 3.0;
  scaling_factor = std::sqrt(average_energy);

  for (int i = 0; i < M; i++) {
    symbols(i) = ((M - 1) - i * 2) / scaling_factor;
    bits2symbols(bin2dec(bitmap.get_row(i))) = i;
  }

  calculate_softbit_matrices();

  setup_done = true;
}


void PAM_c::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "PAM_c::demodulate_bits(): Modulator not ready.");
  int est_symbol;
  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
    est_symbol = round_i((M - 1) - (std::real(signal(i)) * scaling_factor
                                    + (M - 1)) / 2);
    if (est_symbol < 0)
      est_symbol = 0;
    else if (est_symbol > (M - 1))
      est_symbol = M - 1;
    out.replace_mid(i*k, bitmap.get_row(est_symbol));
  }
}

bvec PAM_c::demodulate_bits(const cvec &signal) const
{
  bvec temp(signal.size());
  demodulate_bits(signal, temp);
  return temp;
}


void PAM_c::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                 vec &soft_bits, Soft_Method method) const
{
  it_assert_debug(setup_done, "PAM_c::demodulate_soft_bits(): Modulator not ready.");
  double P0, P1, d0min, d1min, temp;
  vec metric(M);

  soft_bits.set_size(k * rx_symbols.size());

  if (method == LOGMAP) {
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = std::exp(-sqr(std::real(rx_symbols(l) - symbols(j)))
                             / N0);
      }
      for (int i = 0; i < k; i++) {
        P0 = P1 = 0;
        for (int j = 0; j < (M >> 1); j++) {
          P0 += metric(S0(i, j));
          P1 += metric(S1(i, j));
        }
        soft_bits(l*k + i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }
  else { // method == APPROX
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = sqr(std::real(rx_symbols(l) - symbols(j)));
      }
      for (int i = 0; i < k; i++) {
        d0min = d1min = std::numeric_limits<double>::max();
        for (int j = 0; j < (M >> 1); j++) {
          temp = metric(S0(i, j));
          if (temp < d0min) { d0min = temp; }
          temp = metric(S1(i, j));
          if (temp < d1min) { d1min = temp; }
        }
        soft_bits(l*k + i) = (-d0min + d1min) / N0;
      }
    }
  }
}

vec PAM_c::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}


void PAM_c::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
                                 double N0, vec &soft_bits,
                                 Soft_Method method) const
{
  it_assert_debug(setup_done, "PAM_c::demodulate_soft_bits(): Modulator not ready.");
  double P0, P1, d0min, d1min, temp;
  vec metric(M);

  soft_bits.set_size(k * rx_symbols.size());

  if (method == LOGMAP) {
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = std::exp(-sqr(std::real(rx_symbols(l)
                                            - channel(l) * symbols(j))) / N0);
      }
      for (int i = 0; i < k; i++) {
        P0 = P1 = 0;
        for (int j = 0; j < (M >> 1); j++) {
          P0 += metric(S0(i, j));
          P1 += metric(S1(i, j));
        }
        soft_bits(l*k + i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }
  else { // method == APPROX
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = sqr(std::real(rx_symbols(l) - channel(l) * symbols(j)));
      }
      for (int i = 0; i < k; i++) {
        d0min = d1min = std::numeric_limits<double>::max();
        for (int j = 0; j < (M >> 1); j++) {
          temp = metric(S0(i, j));
          if (temp < d0min) { d0min = temp; }
          temp = metric(S1(i, j));
          if (temp < d1min) { d1min = temp; }
        }
        soft_bits(l*k + i) = (-d0min + d1min) / N0;
      }
    }
  }
}

vec PAM_c::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
                                double N0, Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, channel, N0, out, method);
  return out;
}


// ----------------------------------------------------------------------
// PAM
// ----------------------------------------------------------------------

void PAM::set_M(int Mary)
{
  M = Mary;
  k = levels2bits(M);
  it_assert(pow2i(k) == M, "PAM::set_M(): M is not a power of 2");

  symbols.set_size(M, false);
  bits2symbols.set_size(M, false);
  bitmap = graycode(k);
  double average_energy = (sqr(M) - 1) / 3.0;
  scaling_factor = std::sqrt(average_energy);

  for (int i = 0; i < M; i++) {
    symbols(i) = ((M - 1) - i * 2) / scaling_factor;
    bits2symbols(bin2dec(bitmap.get_row(i))) = i;
  }

  calculate_softbit_matrices();

  setup_done = true;
}


void PAM::demodulate_bits(const vec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "PAM::demodulate_bits(): Modulator not ready.");
  int est_symbol;
  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
    est_symbol = round_i((M - 1) - (signal(i) * scaling_factor + (M - 1)) / 2);
    if (est_symbol < 0)
      est_symbol = 0;
    else if (est_symbol > (M - 1))
      est_symbol = M - 1;
    out.replace_mid(i*k, bitmap.get_row(est_symbol));
  }
}

bvec PAM::demodulate_bits(const vec &signal) const
{
  bvec temp(signal.size());
  demodulate_bits(signal, temp);
  return temp;
}

} // namespace itpp
