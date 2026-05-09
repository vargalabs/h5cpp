/*!
 * \file
 * \brief Implementation of signal processing functions
 * \author Tony Ottosson, Thomas Eriksson, Pal Frenger, and Tobias Ringstrom
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

#include <itpp/signal/sigfun.h>
#include <itpp/signal/transforms.h>
#include <itpp/signal/window.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/itcompat.h>
#include <itpp/stat/misc_stat.h>


namespace itpp
{

vec xcorr_old(const vec &x, const int max_lag, const std::string scaleopt)
{
  vec out;
  xcorr_old(x, x, out, max_lag, scaleopt);
  return out;
}

vec xcorr(const vec &x, const int max_lag, const std::string scaleopt)
{
  cvec out(2*x.length() - 1); //Initial size does ont matter, it will get adjusted
  xcorr(to_cvec(x), to_cvec(x), out, max_lag, scaleopt, true);

  return real(out);
}

cvec xcorr(const cvec &x, const int max_lag, const std::string scaleopt)
{
  cvec out(2*x.length() - 1); //Initial size does ont matter, it will get adjusted
  xcorr(x, x, out, max_lag, scaleopt, true);

  return out;
}

vec xcorr(const vec &x, const vec &y, const int max_lag, const std::string scaleopt)
{
  cvec out(2*x.length() - 1); //Initial size does ont matter, it will get adjusted
  xcorr(to_cvec(x), to_cvec(y), out, max_lag, scaleopt, false);

  return real(out);
}

cvec xcorr(const cvec &x, const cvec &y, const int max_lag, const std::string scaleopt)
{
  cvec out(2*x.length() - 1); //Initial size does ont matter, it will get adjusted
  xcorr(x, y, out, max_lag, scaleopt, false);

  return out;
}

void xcorr(const vec &x, const vec &y, vec &out, const int max_lag, const std::string scaleopt)
{
  cvec xx = to_cvec(x);
  cvec yy = to_cvec(y);
  cvec oo = to_cvec(out);
  xcorr(xx, yy, oo, max_lag, scaleopt, false);

  out = real(oo);
}

void xcorr_old(const vec &x, const vec &y, vec &out, const int max_lag, const std::string scaleopt)
{
  int m, n;
  double s_plus, s_minus, M_double, coeff_scale = 0.0;
  int M, N;

  M = x.size();
  M = std::max(x.size(), y.size());
  M_double = double(M);

  if (max_lag == -1) {
    N = std::max(x.size(), y.size());
  }
  else {
    N = max_lag + 1;
  }

  out.set_size(2*N - 1, false);

  it_assert(N <= std::max(x.size(), y.size()), "max_lag cannot be as large as, or larger than, the maximum length of x and y.");

  if (scaleopt == "coeff") {
    coeff_scale = std::sqrt(energy(x)) * std::sqrt(energy(y));
  }

  for (m = 0; m < N; m++) {
    s_plus = 0;
    s_minus = 0;

    for (n = 0;n < M - m;n++) {
      s_minus += index_zero_pad(x, n) * index_zero_pad(y, n + m);
      s_plus += index_zero_pad(x, n + m) * index_zero_pad(y, n);
    }

    if (scaleopt == "none") {
      out(N + m - 1) = s_plus;
      out(N - m - 1) = s_minus;
    }
    else if (scaleopt == "biased") {
      out(N + m - 1) = s_plus / M_double;
      out(N - m - 1) = s_minus / M_double;
    }
    else if (scaleopt == "unbiased") {
      out(N + m - 1) = s_plus / double(M - m);
      out(N - m - 1) = s_minus / double(M - m);
    }
    else if (scaleopt == "coeff") {
      out(N + m - 1) = s_plus / coeff_scale;
      out(N - m - 1) = s_minus / coeff_scale;
    }
    else
      it_error("Incorrect scaleopt specified.");
  }
}


vec xcorr_old(const vec &x, const vec &y, const int max_lag, const std::string scaleopt)
{
  vec out;
  xcorr_old(x, y, out, max_lag, scaleopt);
  return out;
}

//Correlation
void xcorr(const cvec &x, const cvec &y, cvec &out, const int max_lag, const std::string scaleopt, bool autoflag)
{
  int N = std::max(x.length(), y.length());

  //Compute the FFT size as the "next power of 2" of the input vector's length (max)
  int b = ceil_i(::log2(2.0 * N - 1));
  int fftsize = pow2i(b);

  int end = fftsize - 1;

  cvec temp2;
  if (autoflag == true) {
    //Take FFT of input vector
    cvec X = fft(zero_pad(x, fftsize));

    //Compute the abs(X).^2 and take the inverse FFT.
    temp2 = ifft(elem_mult(X, conj(X)));
  }
  else {
    //Take FFT of input vectors
    cvec X = fft(zero_pad(x, fftsize));
    cvec Y = fft(zero_pad(y, fftsize));

    //Compute the crosscorrelation
    temp2 = ifft(elem_mult(X, conj(Y)));
  }

  // Compute the total number of lags to keep. We truncate the maximum number of lags to N-1.
  int maxlag;
  if ((max_lag == -1) || (max_lag >= N))
    maxlag = N - 1;
  else
    maxlag = max_lag;


  //Move negative lags to the beginning of the vector. Drop extra values from the FFT/IFFt
  if (maxlag == 0) {
    out.set_size(1, false);
    out = temp2(0);
  }
  else
    out = concat(temp2(end - maxlag + 1, end), temp2(0, maxlag));


  //Scale data
  if (scaleopt == "biased")
    //out = out / static_cast<double_complex>(N);
    out = out / static_cast<std::complex<double> >(N);
  else if (scaleopt == "unbiased") {
    //Total lag vector
    vec lags = linspace(-maxlag, maxlag, 2 * maxlag + 1);
    cvec scale = to_cvec(static_cast<double>(N) - abs(lags));
    out /= scale;
  }
  else if (scaleopt == "coeff") {
    if (autoflag == true) // Normalize by Rxx(0)
      out /= out(maxlag);
    else { //Normalize by sqrt(Rxx(0)*Ryy(0))
      double rxx0 = sum(abs(elem_mult(x, x)));
      double ryy0 = sum(abs(elem_mult(y, y)));
      out /= std::sqrt(rxx0 * ryy0);
    }
  }
  else if (scaleopt == "none") {}
  else
    it_warning("Unknow scaling option in XCORR, defaulting to <none> ");

}


mat cov(const mat &X, bool is_zero_mean)
{
  int d = X.cols(), n = X.rows();
  mat R(d, d), m2(n, d);
  vec tmp;

  R = 0.0;

  if (!is_zero_mean) {
    // Compute and remove mean
    for (int i = 0; i < d; i++) {
      tmp = X.get_col(i);
      m2.set_col(i, tmp - mean(tmp));
    }

    // Calc corr matrix
    for (int i = 0; i < d; i++) {
      for (int j = 0; j <= i; j++) {
        for (int k = 0; k < n; k++) {
          R(i, j) += m2(k, i) * m2(k, j);
        }
        R(j, i) = R(i, j); // When i=j this is unnecassary work
      }
    }
  }
  else {
    // Calc corr matrix
    for (int i = 0; i < d; i++) {
      for (int j = 0; j <= i; j++) {
        for (int k = 0; k < n; k++) {
          R(i, j) += X(k, i) * X(k, j);
        }
        R(j, i) = R(i, j); // When i=j this is unnecassary work
      }
    }
  }
  R /= n;

  return R;
}

vec spectrum(const vec &v, int nfft, int noverlap)
{
  it_assert_debug(pow2i(levels2bits(nfft)) == nfft,
                  "nfft must be a power of two in spectrum()!");

  vec P(nfft / 2 + 1), w(nfft), wd(nfft);

  P = 0.0;
  w = hanning(nfft);
  double w_energy = nfft == 1 ? 1 : (nfft + 1) * .375; // Hanning energy

  if (nfft > v.size()) {
    P = sqr(abs(fft(to_cvec(elem_mult(zero_pad(v, nfft), w)))(0, nfft / 2)));
    P /= w_energy;
  }
  else {
    int k = (v.size() - noverlap) / (nfft - noverlap), idx = 0;
    for (int i = 0; i < k; i++) {
      wd = elem_mult(v(idx, idx + nfft - 1), w);
      P += sqr(abs(fft(to_cvec(wd))(0, nfft / 2)));
      idx += nfft - noverlap;
    }
    P /= k * w_energy;
  }

  P.set_size(nfft / 2 + 1, true);
  return P;
}

vec spectrum(const vec &v, const vec &w, int noverlap)
{
  int nfft = w.size();
  it_assert_debug(pow2i(levels2bits(nfft)) == nfft,
                  "The window size must be a power of two in spectrum()!");

  vec P(nfft / 2 + 1), wd(nfft);

  P = 0.0;
  double w_energy = energy(w);

  if (nfft > v.size()) {
    P = sqr(abs(fft(to_cvec(elem_mult(zero_pad(v, nfft), w)))(0, nfft / 2)));
    P /= w_energy;
  }
  else {
    int k = (v.size() - noverlap) / (nfft - noverlap), idx = 0;
    for (int i = 0; i < k; i++) {
      wd = elem_mult(v(idx, idx + nfft - 1), w);
      P += sqr(abs(fft(to_cvec(wd))(0, nfft / 2)));
      idx += nfft - noverlap;
    }
    P /= k * w_energy;
  }

  P.set_size(nfft / 2 + 1, true);
  return P;
}

vec filter_spectrum(const vec &a, int nfft)
{
  vec s = sqr(abs(fft(to_cvec(a), nfft)));
  s.set_size(nfft / 2 + 1, true);
  return s;
}

vec filter_spectrum(const vec &a, const vec &b, int nfft)
{
  vec s = sqr(abs(elem_div(fft(to_cvec(a), nfft), fft(to_cvec(b), nfft))));
  s.set_size(nfft / 2 + 1, true);
  return s;
}

} // namespace itpp




