/*!
 * \file
 * \brief Implementation of Filter classes and functions
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

#include <itpp/signal/filter.h>
#include <itpp/signal/window.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/math/trig_hyp.h>


namespace itpp
{


vec filter(const vec &b, const vec &a, const vec &input)
{
  ARMA_Filter<double, double, double> f(b, a);
  return f(input);
}

cvec filter(const vec &b, const vec &a, const cvec &input)
{
  ARMA_Filter<std::complex<double>, double, std::complex<double> > f(b, a);
  return f(input);
}

cvec filter(const cvec &b, const cvec &a, const cvec &input)
{
  ARMA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(b, a);
  return f(input);
}

cvec filter(const cvec &b, const cvec &a, const vec &input)
{
  ARMA_Filter<double, std::complex<double>, std::complex<double> > f(b, a);
  return f(input);
}


vec filter(const vec &b, const int one, const vec &input)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<double, double, double> f(b);
  return f(input);
}

cvec filter(const vec &b, const int one, const cvec &input)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<std::complex<double>, double, std::complex<double> > f(b);
  return f(input);
}

cvec filter(const cvec &b, const int one, const cvec &input)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(b);
  return f(input);
}

cvec filter(const cvec &b, const int one, const vec &input)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<double, std::complex<double>, std::complex<double> > f(b);
  return f(input);
}


vec filter(const int one, const vec &a, const vec &input)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<double, double, double> f(a);
  return f(input);
}

cvec filter(const int one, const vec &a, const cvec &input)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<std::complex<double>, double, std::complex<double> > f(a);
  return f(input);
}

cvec filter(const int one, const cvec &a, const cvec &input)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(a);
  return f(input);
}

cvec filter(const int one, const cvec &a, const vec &input)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<double, std::complex<double>, std::complex<double> > f(a);
  return f(input);
}





vec filter(const vec &b, const vec &a, const vec &input, const vec &state_in, vec &state_out)
{
  ARMA_Filter<double, double, double> f(b, a);
  f.set_state(state_in);
  vec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const vec &b, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out)
{
  ARMA_Filter<std::complex<double>, double, std::complex<double> > f(b, a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const cvec &b, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out)
{
  ARMA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(b, a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const cvec &b, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out)
{
  ARMA_Filter<double, std::complex<double>, std::complex<double> > f(b, a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}


vec filter(const vec &b, const int one, const vec &input, const vec &state_in, vec &state_out)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<double, double, double> f(b);
  f.set_state(state_in);
  vec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const vec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<std::complex<double>, double, std::complex<double> > f(b);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const cvec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(b);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const cvec &b, const int one, const vec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a MA filter a=1");
  MA_Filter<double, std::complex<double>, std::complex<double> > f(b);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}


vec filter(const int one, const vec &a, const vec &input, const vec &state_in, vec &state_out)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<double, double, double> f(a);
  f.set_state(state_in);
  vec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const int one, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<std::complex<double>, double, std::complex<double> > f(a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const int one, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<std::complex<double>, std::complex<double>, std::complex<double> > f(a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

cvec filter(const int one, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out)
{
  it_assert(one == 1, "filter(): in a AR filter b=1");
  AR_Filter<double, std::complex<double>, std::complex<double> > f(a);
  f.set_state(state_in);
  cvec output = f(input);
  state_out = f.get_state();
  return output;
}

vec fir1(int N, double cutoff)
{
  vec a(N + 1), h = hamming(N + 1);

  for (int i = 0;i < length(a);i++) {
    a[i] = h[i] * sinc(cutoff * (i - N / 2.0));
  }
  a /= sum(a);
  return a;
}

template class ITPP_EXPORT MA_Filter<double, double, double>;
template class ITPP_EXPORT MA_Filter < double, std::complex<double>,
std::complex<double> >;
template class ITPP_EXPORT MA_Filter<std::complex<double>, double, std::complex<double> >;
template class ITPP_EXPORT MA_Filter < std::complex<double>, std::complex<double>,
std::complex<double> >;

template class ITPP_EXPORT AR_Filter<double, double, double>;
template class ITPP_EXPORT AR_Filter < double, std::complex<double>,
std::complex<double> >;
template class ITPP_EXPORT AR_Filter < std::complex<double>, double,
std::complex<double> >;
template class ITPP_EXPORT AR_Filter < std::complex<double>, std::complex<double>,
std::complex<double> >;

template class ITPP_EXPORT ARMA_Filter<double, double, double>;
template class ITPP_EXPORT ARMA_Filter < double, std::complex<double>,
std::complex<double> >;
template class ITPP_EXPORT ARMA_Filter < std::complex<double>, double,
std::complex<double> >;
template class ITPP_EXPORT ARMA_Filter < std::complex<double>, std::complex<double>,
std::complex<double> >;

} // namespace itpp
