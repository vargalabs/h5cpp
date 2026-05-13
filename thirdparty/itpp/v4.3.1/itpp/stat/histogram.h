/*!
 * \file
 * \brief Histogram class - header file
 * \author Andy Panov and Adam Piatyszek
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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <itpp/base/mat.h>


namespace itpp
{

//! \addtogroup histogram
//!@{

/*!
  \brief Histogram computation class
  \author Andy Panov

  The Histogram class counts the number of observations of arbitrary
  numerical types that fall into specified bins. Centers of the leftmost
  and rightmost bin along with a total number of bins are passed to a
  histogram object as constructor parameters. Histogram counters are
  updated when calling update() method. It is possible to access bin
  counters and bin interval parameters for all bins at once or separately
  for each bin.

  Example:
  \code
  // Create histogram with 100 bins spanning from 0 to 99 (leftmost bin is
  // centered at 0, rightmost bin is centerd at 99).
  Histogram<double> hist(0, 99, 100);

  // Compute histogram of 100 random variables taken from normal distribution
  hist.update(randn(100));

  // Get position of bin number 5
  double bin5_center = hist.get_bin_center(5);
  // Get corresponding bin counter
  int bin5_counter = hist.get_bin(5);
  // Get bin 5 left boundary:
  double bin5_left = hist.get_bin_left(5);

  // compute PDF & CDF of experimental data
  vec my_data_pdf = hist.get_pdf();
  vec my_data_cdf = hist.get_cdf();
  \endcode
 */
template<typename Num_T>
class Histogram
{
public:
  //! Default constructor. Constructs histogram with 100 bins spanning
  //! values from 0 to 99 by default.
  Histogram(Num_T from = Num_T(0), Num_T to = Num_T(99), int n_bins = 100);
  //! Default destructor
  ~Histogram() {};

  //! Histogram setup
  void setup(Num_T from, Num_T to, int n_bins);

  //! Histogram update
  void update(Num_T value);
  //! Histogram update
  void update(Vec<Num_T> values);
  //! Histogram update
  void update(Mat<Num_T> values);

  //! Bins reset, so accumulation can be restarted
  void reset() { trials_cnt = 0; bins.zeros(); };
  //! Access to single bin counter
  int get_bin(int ix) const { return bins(ix); };
  //! Access to histogram as a vector
  ivec get_bins() const { return bins; };
  //! Access to bin center values (all bins)
  Vec<Num_T> get_bin_centers() const { return center_vals; };
  //! Access to bin center (single bin)
  Num_T get_bin_center(int ix) const { return center_vals(ix); };
  //! Access to left boundary of bin intervals (all bins)
  Vec<Num_T> get_bin_lefts() const { return lo_vals; };
  //! Access to left boundary of single bin
  Num_T get_bin_left(int ix) const { return lo_vals(ix); };
  //! Access to right boundary of bin intervals (all bins)
  Vec<Num_T> get_bin_rights() const { return hi_vals; };
  //! Access to right boundary of single bin
  Num_T get_bin_right(int ix) const { return hi_vals(ix); };

  //! Experimental Probability Density Function (PDF) computation
  vec get_pdf() const;
  //! Experimental Cumulative Density Function (CDF) computation
  vec get_cdf() const;

  //! Current number of bins
  int bins_num() const { return num_bins; };
  //! Current trials counter
  int trials_num() const {return trials_cnt;};

private:
  //! Number of bins
  int num_bins;
  //! Step between bins
  Num_T step;
  //! Low boundaries of histogram bins
  Vec<Num_T> lo_vals;
  //! Upper boundaries of histogram bins
  Vec<Num_T> hi_vals;
  //! Bin centers
  Vec<Num_T> center_vals;
  //! Bins storage
  ivec bins;
  //! Number of processed samples
  int trials_cnt;
};

template<class Num_T>
inline Histogram<Num_T>::Histogram(Num_T from, Num_T to, int n_bins)

{
  setup(from, to, n_bins);
}

template<class Num_T>
inline void Histogram<Num_T>::setup(Num_T from, Num_T to, int n_bins)
{
  num_bins = n_bins;
  lo_vals.set_size(n_bins);
  hi_vals.set_size(n_bins);
  center_vals.set_size(n_bins);
  bins.set_size(n_bins);
  trials_cnt = 0;
  step = (to - from) / (num_bins - 1);
  center_vals = linspace(from, to, num_bins);
  lo_vals = center_vals - step / 2;
  hi_vals = center_vals + step / 2;
  reset();
}

template<class Num_T>
inline void Histogram<Num_T>::update(Num_T value)
{
  // search for the corresponding bin using dichotomy approach
  int start = 0;
  int end = num_bins - 1;
  int test = (start + end) / 2;

  while (start < end) {
    if (value < lo_vals(test))
      end = test - 1;
    else if (value >= hi_vals(test))
      start = test + 1;
    else
      break;
    test = (start + end) / 2;
  };

  bins(test)++;
  trials_cnt++;
}

template<class Num_T>
inline void Histogram<Num_T>::update(Vec<Num_T> values)
{
  for (int i = 0; i < values.length(); i++)
    update(values(i));
}

template<class Num_T>
inline void Histogram<Num_T>::update(Mat<Num_T> values)
{
  for (int i = 0; i < values.rows(); i++)
    for (int j = 0; j < values.cols(); j++)
      update(values(i, j));
}

template<class Num_T>
inline vec Histogram<Num_T>::get_pdf() const
{
  vec pdf(num_bins);
  for (int j = 0; j < num_bins; j++)
    pdf(j) = static_cast<double>(bins(j)) / trials_cnt;
  return pdf;
}

template<class Num_T>
inline vec Histogram<Num_T>::get_cdf() const
{
  ivec tmp = cumsum(bins);
  vec cdf(num_bins);
  for (int j = 0; j < num_bins; j++)
    cdf(j) = static_cast<double>(tmp(j)) / trials_cnt;
  return cdf;
}

//!@}

} // namespace itpp

#endif // #ifndef HISTOGRAM_H
