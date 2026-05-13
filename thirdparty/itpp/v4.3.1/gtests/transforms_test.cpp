/*!
 * \file
 * \brief Transforms test program
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood, Adam Piatyszek, Andy Panov and Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#include <vector>
#include <itpp/itsignal.h>
#include "gtest/gtest.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace itpp;

//set test tolerance (measure of relative and absolute error)
const double max_rel_error = 1e-6;
const double max_abs_error = 1e-6;

//transform results tester
template<typename T>
inline void test_result(const Vec<T>& in, const Vec<T>& ref)
{
  int n = size(in);
  it_assert(n == ref.size(), "compute_rel_error(): input and reference sizes must be equal.");
  for(int i = 0; i < n ; ++i) {
    if(abs(in(i) - ref(i)) < max_abs_error) continue; //handle numbers with absolute value close to zero (relative error can be huge for them)
    double rel_error = abs(in(i) - ref(i)) / abs(in(i));
    ASSERT_LE(rel_error, max_rel_error);
  }
}

//tester for multiple results collected in std::vector
template<typename T>
inline void test_result(const std::vector<Vec<T> >& in, const Vec<T>& ref)
{
  typename std::vector<Vec<T> >::size_type i;
  for(i = 0; i < in.size(); ++i) test_result(in[i], ref);
}


//Reference transform implementations. These functions are not intended to be fast.
//They just strictly follow the transform definitions

//reference DFT implementation
template<typename T>
inline cvec ref_dft(const Vec<T>& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dft(): zero-sized input detected.");
  cvec ret(n);
  for(int i = 0; i < n; ++i) {
    std::complex<double> res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += std::complex<double>(cos(2 * pi * i * j / n), -sin(2 * pi * i * j / n)) * in(j);
    }
    ret(i) = res;
  }
  return ret;
}

//reference IDFT implementation
inline cvec ref_idft(const cvec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_idft(): zero-sized input detected.");
  cvec ret(n);
  for(int i = 0; i < n; ++i) {
    std::complex<double> res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += std::complex<double>(cos(2 * pi * i * j / n), sin(2 * pi * i * j / n)) * in(j);
    }
    ret(i) = res;
  }
  ret *= 1.0 / n;
  return ret;
}

//Type-II DCT reference implementation
inline vec ref_dct(const vec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dct(): zero-sized input detected.");
  vec ret(n);
  for(int i = 0; i < n; ++i) {
    double res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += cos(pi * (j + 0.5) * i / n) * in(j);
    }
    ret(i) = 2 * res;
  }

  // Scale to matlab definition format
  ret /= std::sqrt(2.0 * n);
  ret(0) /= std::sqrt(2.0);

  return ret;
}

//Type-III DCT (IDCT) reference implementation
inline vec ref_idct(const vec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dct(): zero-sized input detected.");
  vec tmp = in;
  tmp(0) *= std::sqrt(2.0);
  tmp /= std::sqrt(2.0 * n);

  vec ret(n);
  for(int i = 0; i < n; ++i) {
    double res = 0.0;
    for(int j = 1; j < n; ++j) {
      res += cos(pi * (i + 0.5) * j / n) * tmp(j);
    }
    ret(i) = 2 * res + tmp(0);
  }

  return ret;
}

//Transform testers run the same transform function in single or multiple threads (if OMP is enabled).
template<typename InType, typename OutType>
std::vector<OutType> run_transform_test(OutType(*transform_function)(const InType&), const InType& test_input)
{
  //select number of results. Multiple threads run the same code if OMP is enabled
#ifdef _OPENMP
  static const int threads_cnt = omp_get_max_threads();
  omp_set_num_threads(threads_cnt);
#else
  static const int threads_cnt = 1;
#endif
  std::vector<OutType> ret(threads_cnt);
  #pragma omp parallel
  {
    //parallel region start. Spawn the threads.
    #pragma omp for
    for(int j = 0; j < threads_cnt; ++j) {
      ret[j] = transform_function(test_input);
    }
    //parallel region end. Join the threads.
  }
  return ret;
}

template<typename InType, typename OutType>
std::vector<OutType> run_transform_test(OutType(*transform_function)(const InType&, int), const InType& test_input, int N)
{
  //select number of results. Multiple threads run the same code if OMP is enabled
#ifdef _OPENMP
  static const int threads_cnt = omp_get_max_threads();
  omp_set_num_threads(threads_cnt);
#else
  static const int threads_cnt = 1;
#endif
  std::vector<OutType> ret(threads_cnt);
  #pragma omp parallel
  {
    //parallel region start. Spawn the threads.
    #pragma omp for
    for(int j = 0; j < threads_cnt; ++j) {
      ret[j] = transform_function(test_input, N);
    }
    //parallel region end. Join the threads.
  }
  return ret;
}

//----------------------------------------------
//Gtest test cases
//----------------------------------------------

TEST(Transforms, FFTReal)
{
  int N = 16;
  vec x = randn(N);

  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  //vector processing test
  {
    SCOPED_TRACE("y = fft_real(x) test");
    //cvec y = fft_real(x);
    std::vector<cvec> y = run_transform_test(fft_real, x);
    test_result(y, ref_dft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = fft_real(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    //cvec y = fft_real(x, N_sub);
    std::vector<cvec> y = run_transform_test(fft_real, x, N_sub);
    test_result(y, ref_dft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = fft_real(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    //cvec y = fft_real(x, N + N_zp);
    std::vector<cvec> y = run_transform_test(fft_real, x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dft(x));
  }
}

TEST(Transforms, IFFTReal)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  //vector processing test
  {
    SCOPED_TRACE("y = ifft_real(x) test");
    int N = 16;
    cvec t = randn_c(N - 1);
    cvec x(N);
    //generate test complex sequence with real spectra
    x.set_subvector(1, 0.5 * (t + conj(reverse(t))));
    x(0) = randn();
    //run transform & test results
    //vec y = ifft_real(x);
    std::vector<vec> y = run_transform_test(ifft_real, x);
    test_result(y, real(ref_idft(x)));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = ifft_real(x, N) test, N < length(x)");
    int N = 16, N_sub = 11; //define odd subvector length to test odd-length transform
    cvec t = randn_c(N - 1);
    cvec x(N);
    //fill subvector samples with Hermitian sequence
    x.set_subvector(1, 0.5 * (t(1, N_sub - 1) + conj(reverse(t(1, N_sub - 1)))));
    x(0) = randn();
    //fill the rest of x with random data (these data should be ignored by IFFT implementation)
    x.set_subvector(N_sub, randn_c(N - N_sub));
    //run transform & test results
    //vec y = ifft_real(x, N_sub);
    std::vector<vec> y = run_transform_test(ifft_real, x, N_sub);
    test_result(y, real(ref_idft(x(0, N_sub - 1))));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = ifft_real(x, N) test, N > length(x)");
    int N_data = 32, N_zp = 8; //define data and zero-padding length
    cvec t = randn_c(N_data);
    cvec x(N_data + N_zp + 1);
    //generate the test data. sequence posesses Hermitian symmetry after zero-padding with N_zp zeros.
    x(0) = randn();
    x.set_subvector(1, N_zp, std::complex<double>(0));
    x.set_subvector(N_zp + 1, 0.5 * (t + conj(reverse(t))));
    //run transform & test results
    //vec y = ifft_real(x, N_data + 2*N_zp + 1);
    std::vector<vec> y = run_transform_test(ifft_real, x, N_data + 2 * N_zp + 1);
    x.set_size(N_data + 2 * N_zp + 1, true);
    test_result(y, real(ref_idft(x)));
  }
}

TEST(Transforms, FFTCplx)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  int N = 16;
  cvec x = randn_c(N);

  //vector processing test
  {
    SCOPED_TRACE("y = fft(x) test");
    //cvec y = fft(x);
    std::vector<cvec> y = run_transform_test(fft, x);
    test_result(y, ref_dft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = fft(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    //cvec y = fft(x, N_sub);
    std::vector<cvec> y = run_transform_test(fft, x, N_sub);
    test_result(y, ref_dft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = fft(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    //cvec y = fft(x, N + N_zp);
    std::vector<cvec> y = run_transform_test(fft, x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dft(x));
  }
}

TEST(Transforms, IFFTCplx)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  int N = 16;
  cvec x = randn_c(N), y;

  //vector processing test
  {
    SCOPED_TRACE("y = ifft(x) test");
    //cvec y = ifft(x);
    std::vector<cvec> y = run_transform_test(ifft, x);
    test_result(y, ref_idft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = ifft(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    //cvec y = ifft(x, N_sub);
    std::vector<cvec> y = run_transform_test(ifft, x, N_sub);
    test_result(y, ref_idft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = ifft(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    //cvec y = ifft(x, N + N_zp);
    std::vector<cvec> y = run_transform_test(ifft, x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_idft(x));
  }
}

TEST(Transforms, DCT)
{
  if(!have_cosine_transforms()) FAIL() << "Cosine Transforms are not supported with this library build.";

  int N = 16;
  vec x = randn(N);

  //vector processing test
  {
    SCOPED_TRACE("y = dct(x) test");
    //vec y = dct(x);
    std::vector<vec> y = run_transform_test(dct, x);
    test_result(y, ref_dct(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = dct(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    //vec y = dct(x, N_sub);
    std::vector<vec> y = run_transform_test(dct, x, N_sub);
    test_result(y, ref_dct(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = dct(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    //vec y = dct(x, N + N_zp);
    std::vector<vec> y = run_transform_test(dct, x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dct(x));
  }
}

TEST(Transforms, IDCT)
{
  if(!have_cosine_transforms()) FAIL() << "Cosine Transforms are not supported with this library build.";

  int N = 16;
  vec x = randn(N);

  //vector processing test
  {
    SCOPED_TRACE("y = idct(x) test");
    //vec y = idct(x);
    std::vector<vec> y = run_transform_test(idct, x);
    test_result(y, ref_idct(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = idct(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    //vec y = idct(x, N_sub);
    std::vector<vec> y = run_transform_test(idct, x, N_sub);
    test_result(y, ref_idct(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = idct(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    //vec y = idct(x, N + N_zp);
    std::vector<vec> y = run_transform_test(idct, x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_idct(x));
  }
}

//Run several transforms sequentially. This test runs three FFT transforms sequentially in each thread.
//The test is intended to verify FFT operation on larger data sets and test for possible clashes on shared
//data in multithreaded environment.

cvec seq_transforms_test(const cvec& test_input)
{
//run 128,256,512-point FFT on test dataset and store results in output vector
  cvec output(128 + 256 + 512);
  output.set_subvector(0, fft(test_input, 128));
  output.set_subvector(128, fft(test_input, 256));
  output.set_subvector(128 + 256, fft(test_input));
  return output;
}

cvec seq_transforms_ref(const cvec& test_input)
{
//compute reference transforms to verify test results
  cvec output(128 + 256 + 512);
  output.set_subvector(0, ref_dft(test_input(0, 127)));
  output.set_subvector(128, ref_dft(test_input(0, 255)));
  output.set_subvector(128 + 256, ref_dft(test_input));
  return output;

}

TEST(Transforms, FFT_128_256_512)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  cvec test_input = randn_c(512);

  std::vector<cvec> y = run_transform_test(seq_transforms_test, test_input);
  test_result(y, seq_transforms_ref(test_input));
}
