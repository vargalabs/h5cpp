/*!
 * \file
 * \brief Filter design functions
 * \author Tony Ottosson
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

#include <itpp/signal/filter_design.h>
#include <itpp/signal/poly.h>
#include <itpp/signal/filter.h>
#include <itpp/signal/transforms.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/algebra/ls_solve.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/math/trig_hyp.h>
#include <itpp/base/converters.h>


namespace itpp
{


void polystab(const vec &a, vec &out)
{
  cvec r;
  roots(a, r);

  for (int i = 0; i < r.size(); i++) {
    if (abs(r(i)) > 1)
      r(i) = std::complex<double>(1.0) / conj(r(i));
  }
  out = real(std::complex<double>(a(0)) * poly(r));
}

void polystab(const cvec &a, cvec &out)
{
  cvec r;
  roots(a, r);

  for (int i = 0; i < r.size(); i++) {
    if (abs(r(i)) > 1)
      r(i) = std::complex<double>(1.0) / conj(r(i));
  }
  out = a(0) * poly(r);
}


// ----------------------- freqz() ---------------------------------------------------------
void freqz(const cvec &b, const cvec &a, const int N, cvec &h, vec &w)
{
  w = pi * linspace(0, N - 1, N) / double(N);

  cvec ha, hb;
  hb = fft(b, 2 * N);
  hb = hb(0, N - 1);

  ha = fft(a, 2 * N);
  ha = ha(0, N - 1);

  h = elem_div(hb, ha);
}

cvec freqz(const cvec &b, const cvec &a, const int N)
{
  cvec h;
  vec w;

  freqz(b, a, N, h, w);

  return h;
}


cvec freqz(const cvec &b, const cvec &a, const vec &w)
{
  int la = a.size(), lb = b.size(), k = std::max(la, lb);

  cvec h, ha, hb;

  // Evaluate the nominator and denominator at the given frequencies
  hb = polyval(zero_pad(b, k), to_cvec(cos(w), sin(w)));
  ha = polyval(zero_pad(a, k), to_cvec(cos(w), sin(w)));

  h = elem_div(hb, ha);

  return h;
}

void freqz(const vec &b, const vec &a, const int N, cvec &h, vec &w)
{
  w = pi * linspace(0, N - 1, N) / double(N);

  cvec ha, hb;
  hb = fft_real(b, 2 * N);
  hb = hb(0, N - 1);

  ha = fft_real(a, 2 * N);
  ha = ha(0, N - 1);

  h = elem_div(hb, ha);
}

cvec freqz(const vec &b, const vec &a, const int N)
{
  cvec h;
  vec w;

  freqz(b, a, N, h, w);

  return h;
}


cvec freqz(const vec &b, const vec &a, const vec &w)
{
  int la = a.size(), lb = b.size(), k = std::max(la, lb);

  cvec h, ha, hb;

  // Evaluate the nominator and denominator at the given frequencies
  hb = polyval(zero_pad(b, k), to_cvec(cos(w), sin(w)));
  ha = polyval(zero_pad(a, k), to_cvec(cos(w), sin(w)));

  h = elem_div(hb, ha);

  return h;
}



void filter_design_autocorrelation(const int N, const vec &f, const vec &m, vec &R)
{
  it_assert(f.size() == m.size(), "filter_design_autocorrelation: size of f and m vectors does not agree");
  int N_f = f.size();

  it_assert(f(0) == 0.0, "filter_design_autocorrelation: first frequency must be 0.0");
  it_assert(f(N_f - 1) == 1.0, "filter_design_autocorrelation: last frequency must be 1.0");

  // interpolate frequency-response
  int N_fft = 512;
  vec m_interp(N_fft + 1);
  // unused variable:
  // double df_interp = 1.0/double(N_fft);

  m_interp(0) = m(0);
  double inc;

  int jstart = 0, jstop;

  for (int i = 0; i < N_f - 1; i++) {
    // calculate number of points to the next frequency
    jstop = floor_i(f(i + 1) * (N_fft + 1)) - 1;
    //std::cout << "jstart=" << jstart << "jstop=" << jstop << std::endl;

    for (int j = jstart; j <= jstop; j++) {
      inc = double(j - jstart) / double(jstop - jstart);
      m_interp(j) = m(i) * (1 - inc) + m(i + 1) * inc;
    }
    jstart = jstop + 1;
  }

  vec S = sqr(concat(m_interp, reverse(m_interp(2, N_fft))));  // create a complete frequency response with also negative frequencies

  R = ifft_real(to_cvec(S)); // calculate correlation

  R = R.left(N);
}


// Calculate the AR coefficients of order \c n of the ARMA-process defined by the autocorrelation R
// using the deternined modified Yule-Walker method
// maxlag determines the size of the system to solve N>= n.
// If N>m then the system is overdetermined and a least squares solution is used.
// as a rule of thumb use N = 4*n
void modified_yule_walker(const int m, const int n, const int N, const vec &R, vec &a)
{
  it_assert(m > 0, "modified_yule_walker: m must be > 0");
  it_assert(n > 0, "modified_yule_walker: n must be > 0");
  it_assert(N <= R.size(), "modified_yule_walker: autocorrelation function too short");

  // create the modified Yule-Walker equations Rm * a = - rh
  // see eq. (3.7.1) in Stoica and Moses, Introduction to spectral analysis
  int M = N - m - 1;

  mat Rm;
  if (m - n + 1 < 0)
    Rm = toeplitz(R(m, m + M - 1), reverse(concat(R(1, std::abs(m - n + 1)), R(0, m))));
  else
    Rm = toeplitz(R(m, m + M - 1), reverse(R(m - n + 1, m)));


  vec rh = - R(m + 1, m + M);

  // solve overdetermined system
  a = backslash(Rm, rh);

  // prepend a_0 = 1
  a = concat(1.0, a);

  // stabilize polynomial
  a = polystab(a);
}



void arma_estimator(const int m, const int n, const vec &R, vec &b, vec &a)
{
  it_assert(m > 0, "arma_estimator: m must be > 0");
  it_assert(n > 0, "arma_estimator: n must be > 0");
  it_assert(2*(m + n) <= R.size(), "arma_estimator: autocorrelation function too short");


  // windowing the autocorrelation
  int N = 2 * (m + n);
  vec Rwindow = elem_mult(R.left(N), 0.54 + 0.46 * cos(pi * linspace(0.0, double(N - 1), N) / double(N - 1))); // Hamming windowing

  // calculate the AR part using the overdetmined Yule-Walker equations
  modified_yule_walker(m, n, N, Rwindow, a);

  // --------------- Calculate MA part --------------------------------------
  // use method in ref [2] section VII.
  vec r_causal = Rwindow;
  r_causal(0) *= 0.5;

  vec h_inv_a = filter(1, a, concat(1.0, zeros(N - 1))); // see eq (50) of [2]
  mat H_inv_a = toeplitz(h_inv_a, concat(1.0, zeros(m)));

  vec b_causal = backslash(H_inv_a, r_causal);

  // calculate the double-sided spectrum
  int N_fft = 256;
  vec H = 2.0 * real(elem_div(fft_real(b_causal, N_fft), fft_real(a, N_fft))); // calculate spectrum

  // Do weighting and windowing in cepstrum domain
  cvec cepstrum = log(to_cvec(H));
  cvec q = ifft(cepstrum);

  // keep only causal part of spectrum (windowing)
  q.set_subvector(N_fft / 2, zeros_c(N_fft / 2));
  q(0) *= 0.5;

  cvec h = ifft(exp(fft(q))); // convert back to frequency domain, from cepstrum and do inverse transform to calculate impulse response
  b = real(backslash(to_cmat(H_inv_a), h(0, N - 1))); // use Shank's method to calculate b coefficients
}


void yulewalk(const int N, const vec &f, const vec &m, vec &b, vec &a)
{
  it_assert(f.size() == m.size(), "yulewalk: size of f and m vectors does not agree");
  int N_f = f.size();

  it_assert(f(0) == 0.0, "yulewalk: first frequency must be 0.0");
  it_assert(f(N_f - 1) == 1.0, "yulewalk: last frequency must be 1.0");


  vec R;
  filter_design_autocorrelation(4*N, f, m, R);

  arma_estimator(N, N, R, b, a);
}


} // namespace itpp
