/*!
 * \file
 * \brief Definitions of Filter classes and functions
 * \author Hakan Eriksson, Thomas Eriksson, Tony Ottosson and Adam Piatyszek
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

#ifndef FILTER_H
#define FILTER_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup filters
*/

/*!
  \brief Virtual Filter Base Class.
  \ingroup filters

  The class is templated as follows:
  <ul>
  <li> \c T1 is the type of the input samples</li>
  <li> \c T2 is the type of the filter coefficients</li>
  <li> \c T3 is the type of the output samples</li>
  </ul>
*/
template <class T1, class T2, class T3>
class Filter
{
public:
  //! Default constructor
  Filter() {}
  //! Filter a single sample.
  virtual T3 operator()(const T1 Sample) { return filter(Sample); }
  //! Filter a vector.
  virtual Vec<T3> operator()(const Vec<T1> &v);
  //! Virtual destructor.
  virtual ~Filter() {}
protected:
  /*!
    \brief Pure virtual filter function. This is where the real
    filtering is done. Implement this function to create a new filter.
  */
  virtual T3 filter(const T1 Sample) = 0;
};

/*!
  \brief Moving Average Filter Base Class.
  \ingroup filters

  This class implements a moving average (MA) filter
  according to
  \f[
  y(n) = b(0)*x(n) + b(1)*x(n-1) + ... + b(N)*x(n-N)
  \f]
  where \a b is the filter coefficients, \a x is the input and \a y
  is the output.

  When filtering a vector, the length of the output vector equals
  the length of the input vector.  Internal states are kept in a
  filter memory. The first time the filter is used the internal
  states have been set to zero.

  The class is templated as follows:
  <ul>
  <li> \c T1 is the type of the input samples</li>
  <li> \c T2 is the type of the filter coefficients</li>
  <li> \c T3 is the type of the output samples</li>
  </ul>
*/
template <class T1, class T2, class T3>
class MA_Filter : public Filter<T1, T2, T3>
{
public:
  //! Class default constructor
  explicit MA_Filter();
  //! Class constructor setting the coefficients in the filter
  explicit MA_Filter(const Vec<T2> &b);
  //! Class destructor
  virtual ~MA_Filter() { }
  //! Filter coefficient access function
  Vec<T2> get_coeffs() const { return coeffs; }
  //! Set the filter coefficients
  void set_coeffs(const Vec<T2> &b);
  //! Clears the filter memory
  void clear() { mem.clear(); }
  //! Get state of filter
  Vec<T3> get_state() const;
  //! Set state of filter
  void set_state(const Vec<T3> &state);

private:
  virtual T3 filter(const T1 Sample);

  Vec<T3> mem;
  Vec<T2> coeffs;
  int inptr;
  bool init;
};

/*!
  \brief Autoregressive (AR) Filter Base Class.
  \ingroup filters

  This class implements a autoregressive (AR) filter
  according to
  \f[
  a(0)*y(n) = x(n) - a(1)*y(n-1) - ... - a(N)*y(n-N)
  \f]
  where \a a is the filter coefficients, \a x is the input and \a y
  is the output.

  When filtering a vector, the length of the output vector equals
  the length of the input vector.  Internal states are kept in a
  filter memory. The first time the filter is used the internal
  states have been set to zero.

  The class is templated as follows:
  <ul>
  <li> \c T1 is the type of the input samples</li>
  <li> \c T2 is the type of the filter coefficients</li>
  <li> \c T3 is the type of the output samples</li>
  </ul>
*/
template <class T1, class T2, class T3>
class AR_Filter : public Filter<T1, T2, T3>
{
public:
  //! Class constructor
  explicit AR_Filter();
  //! Class constructor setting the coefficients in the filter
  explicit AR_Filter(const Vec<T2> &a);
  //! Class destructor
  virtual ~AR_Filter() { }
  //! Filter coefficient access function
  Vec<T2> get_coeffs() const { return coeffs; }
  //! Set the filter coefficients (and order)
  void set_coeffs(const Vec<T2> &a);
  //! Clears the filter memory
  void clear() { mem.clear(); }
  //! Get state of filter
  Vec<T3> get_state() const;
  //! Set state of filter
  void set_state(const Vec<T3> &state);

private:
  virtual T3 filter(const T1 Sample);

  Vec<T3> mem;
  Vec<T2> coeffs;
  Vec<T2> a0;
  int inptr;
  bool init;
};


/*!
  \brief Autoregressive Moving Average (ARMA) Filter Base Class.
  \ingroup filters

  This class implements a autoregressive moving average (ARMA) filter
  according to
  \f[
  a(0)*y(n) = b(0)*x(n) + b(1)*x(n-1) + \ldots + b(N_b)*x(n-N_b)
    - a(1)*y(n-1) - \ldots - a(N_a)*y(n-N_a)
  \f]

  where \a a and \a b are the filter coefficients, \a x is the input
  and \a y is the output.

  When filtering a vector, the length of the output vector equals
  the length of the input vector.  Internal states are kept in a
  filter memory. The first time the filter is used the internal
  states have been set to zero.

  The class is templated as follows:
  <ul>
  <li> \c T1 is the type of the input samples</li>
  <li> \c T2 is the type of the filter coefficients</li>
  <li> \c T3 is the type of the output samples</li>
  </ul>
*/
template <class T1, class T2, class T3>
class ARMA_Filter : public Filter<T1, T2, T3>
{
public:
  //! Class constructor
  explicit ARMA_Filter();
  //! Class constructor setting the coefficients in the filter
  explicit ARMA_Filter(const Vec<T2> &b, const Vec<T2> &a);
  //! Class destructor
  virtual ~ARMA_Filter() { }
  //! Filter \a a coefficient access function
  Vec<T2> get_coeffs_a() const { return acoeffs; }
  //! Filter \a b coefficient access function
  Vec<T2> get_coeffs_b() const { return bcoeffs; }
  //! Filter coefficient access function
  void get_coeffs(Vec<T2> &b, Vec<T2> &a) const { b = bcoeffs; a = acoeffs; }
  //! Set the filter coefficients (and order)
  void set_coeffs(const Vec<T2> &b, const Vec<T2> &a);
  //! Clears the filter memory
  void clear() { mem.clear(); }
  //! Get state of filter
  Vec<T3> get_state() const;
  //! Set state of filter
  void set_state(const Vec<T3> &state);

private:
  virtual T3 filter(const T1 Sample);

  Vec<T3> mem;
  Vec<T2> acoeffs, bcoeffs;
  int inptr;
  bool init;
};



/*!
  \brief ARMA filter function
  \ingroup filters

  These functions implements a autoregressive moving average (ARMA) filter
  according to
  \f[
  a(0)*y(n) = b(0)*x(n) + b(1)*x(n-1) + \ldots + b(N_b)*x(n-N_b)
    - a(1)*y(n-1) - \ldots - a(N_a)*y(n-N_a)
  \f]

  where \a a and \a b are the filter coefficients, \a x is the input
  and \a y is the output.

Setting a=1 gives a MA filter and b=1 gives a AR filter. The
  length of the output vector equals the length of the input
  vector. The state vectors \a state_in and \a state_out is of
  length \f$max(N_a, n_b) - 1\f$.

  If no start state \a state_in is given it is set to zero.

  @{
*/
ITPP_EXPORT vec filter(const vec &b, const vec &a, const vec &input);
ITPP_EXPORT cvec filter(const vec &b, const vec &a, const cvec &input);
ITPP_EXPORT cvec filter(const cvec &b, const cvec &a, const cvec &input);
ITPP_EXPORT cvec filter(const cvec &b, const cvec &a, const vec &input);

ITPP_EXPORT vec filter(const vec &b, const int one, const vec &input);
ITPP_EXPORT cvec filter(const vec &b, const int one, const cvec &input);
ITPP_EXPORT cvec filter(const cvec &b, const int one, const cvec &input);
ITPP_EXPORT cvec filter(const cvec &b, const int one, const vec &input);

ITPP_EXPORT vec filter(const int one, const vec &a, const vec &input);
ITPP_EXPORT cvec filter(const int one, const vec &a, const cvec &input);
ITPP_EXPORT cvec filter(const int one, const cvec &a, const cvec &input);
ITPP_EXPORT cvec filter(const int one, const cvec &a, const vec &input);


ITPP_EXPORT vec filter(const vec &b, const vec &a, const vec &input, const vec &state_in, vec &state_out);
ITPP_EXPORT cvec filter(const vec &b, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const cvec &b, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const cvec &b, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out);

ITPP_EXPORT vec filter(const vec &b, const int one, const vec &input, const vec &state_in, vec &state_out);
ITPP_EXPORT cvec filter(const vec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const cvec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const cvec &b, const int one, const vec &input, const cvec &state_in, cvec &state_out);

ITPP_EXPORT vec filter(const int one, const vec &a, const vec &input, const vec &state_in, vec &state_out);
ITPP_EXPORT cvec filter(const int one, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const int one, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out);
ITPP_EXPORT cvec filter(const int one, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out);
/*! @} */


/*!
  \brief Design a Nth order FIR filter with cut-off frequency \c
  cutoff using the window method.
  \ingroup filters
*/
ITPP_EXPORT vec fir1(int N, double cutoff);

//----------------------------------------------------------------------------
// Implementation of templated functions starts here
//----------------------------------------------------------------------------

//---------------------- class Filter ----------------------------

template <class T1, class T2, class T3>
Vec<T3> Filter<T1, T2, T3>::operator()(const Vec<T1> &x)
{
  Vec<T3> y(x.length());

  for (int i = 0; i < x.length(); i++) {
    y[i] = filter(x[i]);
  }

  return y;
}

//-------------------------- class MA_Filter ---------------------------------

template <class T1, class T2, class T3>
MA_Filter<T1, T2, T3>::MA_Filter() : Filter<T1, T2, T3>()
{
  inptr = 0;
  init = false;
}

template <class T1, class T2, class T3>
MA_Filter<T1, T2, T3>::MA_Filter(const Vec<T2> &b) : Filter<T1, T2, T3>()
{
  set_coeffs(b);
}


template <class T1, class T2, class T3>
void MA_Filter<T1, T2, T3>::set_coeffs(const Vec<T2> &b)
{
  it_assert(b.size() > 0, "MA_Filter: size of filter is 0!");

  coeffs = b;
  mem.set_size(coeffs.size(), false);
  mem.clear();
  inptr = 0;
  init = true;
}

template <class T1, class T2, class T3>
Vec<T3> MA_Filter<T1, T2, T3>::get_state() const
{
  it_assert(init == true, "MA_Filter: filter coefficients are not set!");

  int offset = inptr;
  Vec<T3> state(mem.size());

  for (int n = 0; n < mem.size(); n++) {
    state(n) = mem(offset);
    offset = (offset + 1) % mem.size();
  }

  return state;
}

template <class T1, class T2, class T3>
void MA_Filter<T1, T2, T3>::set_state(const Vec<T3> &state)
{
  it_assert(init == true, "MA_Filter: filter coefficients are not set!");
  it_assert(state.size() == mem.size(), "MA_Filter: Invalid state vector!");

  mem = state;
  inptr = 0;
}

template <class T1, class T2, class T3>
T3 MA_Filter<T1, T2, T3>::filter(const T1 Sample)
{
  it_assert(init == true, "MA_Filter: Filter coefficients are not set!");
  T3 s = 0;

  mem(inptr) = Sample;
  int L = mem.length() - inptr;

  for (int i = 0; i < L; i++) {
    s += coeffs(i) * mem(inptr + i);
  }
  for (int i = 0; i < inptr; i++) {
    s += coeffs(L + i) * mem(i);
  }

  inptr--;
  if (inptr < 0)
    inptr += mem.length();

  return s;
}

//---------------------- class AR_Filter ----------------------------------

template <class T1, class T2, class T3>
AR_Filter<T1, T2, T3>::AR_Filter() : Filter<T1, T2, T3>()
{
  inptr = 0;
  init = false;
}

template <class T1, class T2, class T3>
AR_Filter<T1, T2, T3>::AR_Filter(const Vec<T2> &a) : Filter<T1, T2, T3>()
{
  set_coeffs(a);
}

template <class T1, class T2, class T3>
void AR_Filter<T1, T2, T3>::set_coeffs(const Vec<T2> &a)
{
  it_assert(a.size() > 0, "AR_Filter: size of filter is 0!");
  it_assert(a(0) != T2(0), "AR_Filter: a(0) cannot be 0!");

  a0.set_size(1);//needed to keep the first coefficient for future reuse
  a0(0) = a(0);
  coeffs = a / a0(0);

  mem.set_size(coeffs.size() - 1, false);
  mem.clear();
  inptr = 0;
  init = true;
}


template <class T1, class T2, class T3>
Vec<T3> AR_Filter<T1, T2, T3>::get_state() const
{
  it_assert(init == true, "AR_Filter: filter coefficients are not set!");

  int offset = inptr;
  Vec<T3> state(mem.size());

  for (int n = 0; n < mem.size(); n++) {
    state(n) = mem(offset);
    offset = (offset + 1) % mem.size();
  }

  return state;
}

template <class T1, class T2, class T3>
void AR_Filter<T1, T2, T3>::set_state(const Vec<T3> &state)
{
  it_assert(init == true, "AR_Filter: filter coefficients are not set!");
  it_assert(state.size() == mem.size(), "AR_Filter: Invalid state vector!");

  mem = state;
  inptr = 0;
}

template <class T1, class T2, class T3>
T3 AR_Filter<T1, T2, T3>::filter(const T1 Sample)
{
  it_assert(init == true, "AR_Filter: Filter coefficients are not set!");
  T3 s = Sample;

  if (mem.size() == 0)
    return (s / a0(0));

  int L = mem.size() - inptr;
  for (int i = 0; i < L; i++) {
    s -= mem(i + inptr) * coeffs(i + 1); // All coeffs except a(0)
  }
  for (int i = 0; i < inptr; i++) {
    s -= mem(i) * coeffs(L + i + 1); // All coeffs except a(0)
  }

  inptr--;
  if (inptr < 0)
    inptr += mem.size();
  mem(inptr) = s;

  return (s / a0(0));
}


//---------------------- class ARMA_Filter ----------------------------------
template <class T1, class T2, class T3>
ARMA_Filter<T1, T2, T3>::ARMA_Filter() : Filter<T1, T2, T3>()
{
  inptr = 0;
  init = false;
}

template <class T1, class T2, class T3>
ARMA_Filter<T1, T2, T3>::ARMA_Filter(const Vec<T2> &b, const Vec<T2> &a) : Filter<T1, T2, T3>()
{
  set_coeffs(b, a);
}

template <class T1, class T2, class T3>
void ARMA_Filter<T1, T2, T3>::set_coeffs(const Vec<T2> &b, const Vec<T2> &a)
{
  it_assert(a.size() > 0 && b.size() > 0, "ARMA_Filter: size of filter is 0!");
  it_assert(a(0) != T2(0), "ARMA_Filter: a(0) cannot be 0!");

  acoeffs = a / a(0);
  bcoeffs = b / a(0);

  mem.set_size(std::max(a.size(), b.size()) - 1, false);
  mem.clear();
  inptr = 0;
  init = true;
}

template <class T1, class T2, class T3>
Vec<T3> ARMA_Filter<T1, T2, T3>::get_state() const
{
  it_assert(init == true, "ARMA_Filter: filter coefficients are not set!");

  int offset = inptr;
  Vec<T3> state(mem.size());

  for (int n = 0; n < mem.size(); n++) {
    state(n) = mem(offset);
    offset = (offset + 1) % mem.size();
  }

  return state;
}

template <class T1, class T2, class T3>
void ARMA_Filter<T1, T2, T3>::set_state(const Vec<T3> &state)
{
  it_assert(init == true, "ARMA_Filter: filter coefficients are not set!");
  it_assert(state.size() == mem.size(), "ARMA_Filter: Invalid state vector!");

  mem = state;
  inptr = 0;
}

template <class T1, class T2, class T3>
T3 ARMA_Filter<T1, T2, T3>::filter(const T1 Sample)
{
  it_assert(init == true, "ARMA_Filter: Filter coefficients are not set!");
  T3 z = Sample;
  T3 s;

  for (int i = 0; i < acoeffs.size() - 1; i++) { // All AR-coeff except a(0).
    z -= mem((i + inptr) % mem.size()) * acoeffs(i + 1);
  }
  s = z * bcoeffs(0);

  for (int i = 0; i < bcoeffs.size() - 1; i++) { // All MA-coeff except b(0).
    s += mem((i + inptr) % mem.size()) * bcoeffs(i + 1);
  }

  inptr--;
  if (inptr < 0)
    inptr += mem.size();
  mem(inptr) = z;

  mem(inptr) = z; // Store in the internal state.

  return s;
}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT MA_Filter<double, double, double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT MA_Filter< double, std::complex<double>,
  std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT MA_Filter< std::complex<double>, double,
  std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT MA_Filter< std::complex<double>, std::complex<double>,
  std::complex<double> >;

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT AR_Filter<double, double, double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT AR_Filter< double, std::complex<double>,
  std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT AR_Filter< std::complex<double>,
  double, std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT AR_Filter< std::complex<double>, std::complex<double>,
  std::complex<double> >;

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT ARMA_Filter<double, double, double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT ARMA_Filter< double, std::complex<double>,
  std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT ARMA_Filter< std::complex<double>,
  double, std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT ARMA_Filter< std::complex<double>, std::complex<double>,
  std::complex<double> >;

//! \endcond

} // namespace itpp

#endif // #ifndef FILTER_H
