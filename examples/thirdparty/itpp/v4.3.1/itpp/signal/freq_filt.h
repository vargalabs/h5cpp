/*!
 * \file
 * \brief Definition of frequency domain filter class
 * \author Simon Wood
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

#ifndef FREQ_FILT_H
#define FREQ_FILT_H

#include <itpp/base/vec.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/math/min_max.h>
#include <itpp/signal/transforms.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \brief Freq_Filt Frequency domain filtering using the overlap-add technique
  \ingroup filters

  The Freq_Filt class implements an FFT based filter using the overlap-add
  technique. The data is filtered by first transforming the input sequence
  into the frequency domain with an efficient FFT implementation (i.e. FFTW)
  and then multiplied with a Fourier transformed version of the impulse
  response. The resulting data is then inversed Fourier transformed to return
  a filtered time domain signal.

  Freq_Filt is a templated class. The template paramter \c Num_T defines the
  data type for the impulse response \c b, input data \c x and output data
  \c y.

  The class constructor chooses an optimal FFT length and data block
  size that minimizes execution time.

  For example,
  \code
  vec b = "1 2 3 4";
  Freq_Filt<double> FF(b,8000);
  \endcode

  where 8000 corresponds to the expected vector length of the data
  to be filtered. The actual number of elements does not have to be
  exact, but it should be close.

  Here is a complete example:
  \code
  vec b = "1 2 3 4";
  vec x(20);
  x(0) = 1;

  // Define a filter object for doubles
  Freq_Filt<double> FF(b,x.length());

  // Filter the data
  vec y = FF.filter(x);

  // Check the FFT and block sizes that were used
  int fftsize = FF.getFFTSize();
  int blk = FF.getBlkSize();
  \endcode

  To facilitate large data sets the Freq_Filt class has a streaming
  option. In this mode of operation data history is internally
  stored. This allows the class to be used for large data sets that
  are read from disk or streamed in real-time.

  \code
  bool stream = true;
  vec b = "1 2 3 4";
  Freq_Filt<double> FF(b,1000);

  vec x,y;
  while(!EOF) {
  x = "read buffer of data";
  y = FF.filter(x,stream);
  cout << << endl;
  }
  \endcode
*/


template<class Num_T>
class Freq_Filt
{
public:
  /*!
    \brief Constructor

    The empty constructor makes it possible to have other container objects
    of the Freq_Filt class
  */
  Freq_Filt() {}

  /*!
    \brief Constructor with initialization of the impulse response \b b.

    Create a filter object with impulse response \b b. The FFT size and
    data block size are also initialized.
    \code
    vec b = "1 2 3 4";
    vec x(20);
    Freq_Filt FF(b,x.length());
    \endcode
  */
  Freq_Filt(const Vec<Num_T> &b, const int xlength) {init(b, xlength);}

  /*!
    \brief Filter data in the input vector \b x

    Filters data in batch mode or streaming mode
    \code
    FF.filter(x); // Filters data in batch mode
    FF.filter(x,1); // Filters data in streaming mode
    \endcode
  */
  Vec<Num_T> filter(const Vec<Num_T> &x, const int strm = 0);

  //! Return FFT size
  int get_fft_size() { return fftsize; }

  //! Return the data block size
  int get_blk_size() { return blksize; }

  //! Destructor
  ~Freq_Filt() {}

private:
  int fftsize, blksize;
  cvec B; // FFT of impulse vector
  Vec<Num_T> impulse;
  Vec<Num_T> old_data;
  cvec zfinal;

  void init(const Vec<Num_T> &b, const int xlength);
  vec overlap_add(const vec &x);
  svec overlap_add(const svec &x);
  ivec overlap_add(const ivec &x);
  cvec overlap_add(const cvec &x);
  void overlap_add(const cvec &x, cvec &y);
};


// Initialize impulse rsponse, FFT size and block size
template <class Num_T>
void Freq_Filt<Num_T>::init(const Vec<Num_T> &b, const int xlength)
{
  // Set the impulse response
  impulse = b;

  // Get the length of the impulse response
  int num_elements = impulse.length();

  // Initizlize old data
  old_data.set_size(0, false);

  // Initialize the final state
  zfinal.set_size(num_elements - 1, false);
  zfinal.zeros();

  vec Lvec;

  /*
   * Compute the FFT size and the data block size to use.
   * The FFT size is N = L + M -1, where L corresponds to
   * the block size (to be determined) and M corresponds
   * to the impulse length. The method used here is designed
   * to minimize the (number of blocks) * (number of flops per FFT)
   * Use the FFTW flop computation of 5*N*log2(N).
   */
  vec n = linspace(1, 20, 20);
  n = pow(2.0, n);
  ivec fftflops = to_ivec(elem_mult(5.0 * n, log2(n)));

  // Find where the FFT sizes are > (num_elements-1)
  //ivec index = find(n,">", static_cast<double>(num_elements-1));
  ivec index(n.length());
  int cnt = 0;
  for (int ii = 0; ii < n.length(); ii++) {
    if (n(ii) > (num_elements - 1)) {
      index(cnt) = ii;
      cnt += 1;
    }
  }
  index.set_size(cnt, true);

  fftflops = fftflops(index);
  n = n(index);

  // Minimize number of blocks * number of flops per FFT
  Lvec = n - (double)(num_elements - 1);
  int min_ind = min_index(elem_mult(ceil((double)xlength / Lvec), to_vec(fftflops)));
  fftsize = static_cast<int>(n(min_ind));
  blksize = static_cast<int>(Lvec(min_ind));

  // Finally, compute the FFT of the impulse response
  B = fft(to_cvec(impulse), fftsize);
}

// Filter data
template <class Num_T>
Vec<Num_T> Freq_Filt<Num_T>::filter(const Vec<Num_T> &input, const int strm)
{
  Vec<Num_T> x, tempv;
  Vec<Num_T> output;

  /*
   * If we are not in streaming mode we want to process all of the data.
   * However, if we are in streaming mode we want to process the data in
   * blocks that are commensurate with the designed 'blksize' parameter.
   * So, we do a little book keeping to divide the iinput vector into the
   * correct size.
   */
  if (!strm) {
    x = input;
    zfinal.zeros();
    old_data.set_size(0, false);
  }
  else { // we are in streaming mode
    tempv = concat(old_data, input);
    if (tempv.length() <= blksize) {
      x = tempv;
      old_data.set_size(0, false);
    }
    else {
      int end = tempv.length();
      int numblks = end / blksize;
      if ((end % blksize)) {
        x = tempv(0, blksize * numblks - 1);
        old_data = tempv(blksize * numblks, end - 1);
      }
      else {
        x = tempv(0, blksize * numblks - 1);
        old_data.set_size(0, false);
      }
    }
  }
  output = overlap_add(x);

  return output;
}

// Overlap-add routine
template<class Num_T>
void Freq_Filt<Num_T>::overlap_add(const cvec&x, cvec &y)
{
  int nb = impulse.length();
  int nx = x.length();

  y.set_size(nx, false);
  y.zeros();
  cvec X, Y;
  int istart = 0;
  int L = blksize;
  while (istart < nx) {
    int iend = std::min(istart + L - 1, nx - 1);

    X = fft(x(istart, iend), fftsize);
    Y = ifft(elem_mult(X, B));
    Y.set_subvector(0, Y(0, nb - 2) + zfinal);
    int yend = std::min(nx - 1, istart + fftsize - 1);
    y.set_subvector(istart, Y(0, yend - istart));
    zfinal = Y(fftsize - (nb - 1), fftsize - 1);
    istart += L;
  }
}

template<class Num_T>
vec Freq_Filt<Num_T>::overlap_add(const vec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return real(y);
}

template<class Num_T>
svec Freq_Filt<Num_T>::overlap_add(const svec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return to_svec(real(y));
}

template<class Num_T>
ivec Freq_Filt<Num_T>::overlap_add(const ivec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return to_ivec(real(y));
}

template<class Num_T>
cvec Freq_Filt<Num_T>::overlap_add(const cvec &x)
{
  cvec y; // Size of y is set later
  overlap_add(x, y);
  return y;
}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Freq_Filt<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Freq_Filt<std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Freq_Filt<short>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Freq_Filt<int>;

//! \endcond

} // namespace itpp

#endif // #ifndef FREQ_FILT_H
