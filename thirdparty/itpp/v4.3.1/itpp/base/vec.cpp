/*!
 * \file
 * \brief Templated Vector Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined (HAVE_BLAS)
#  include <itpp/base/blas.h>
#endif

#include <itpp/base/vec.h>
#include <itpp/base/converters.h>
#include <cstdio>
#include <limits>
//! \cond

namespace itpp
{

template<>
int Vec<int>::parse_token(const std::string &s) const
{
  int out;
  std::istringstream buffer(s);
  if (s.find('x', 1) != std::string::npos) {
    buffer >> std::hex >> out;
  }
  else if (((s[0] == '0')
            || (((s[0] == '-') || (s[0] == '+')) && (s[1] == '0')))
           && (s.find('8', 1) == std::string::npos)
           && (s.find('9', 1) == std::string::npos)) {
    buffer >> std::oct >> out;
  }
  else {
    buffer >> std::dec >> out;
  }
  return out;
}

template<>
double Vec<double>::parse_token(const std::string &s) const
{
  double out;
  if ((s == "NaN") || (s == "nan") || (s == "NAN")) {
    if (std::numeric_limits<double>::has_quiet_NaN)
      out = std::numeric_limits<double>::quiet_NaN();
    else if (std::numeric_limits<double>::has_signaling_NaN)
      out = std::numeric_limits<double>::signaling_NaN();
    else
      it_error("Vec<double>::set(): NaN not supported");
  }
  else if ((s =="-Inf") || (s =="-inf") || (s =="-INF")) {
    out = -std::numeric_limits<double>::infinity();
  }
  else if ((s =="Inf") || (s =="inf") || (s =="INF") ||
           (s =="+Inf") || (s =="+inf") || (s =="+INF")) {
    out = std::numeric_limits<double>::infinity();
  }
  else {
    std::istringstream buffer(s);
    buffer >> out;
    it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation failed "
              "(buffer >> out)");
  }
  return out;
}


template<>
void Vec<bin>::set(const std::string &str)
{
  bool abc_format;
  std::vector<std::string> tokens = tokenize(str, abc_format);
  it_assert(!abc_format, "Vec<bin>::set(): \"a:b:c\" format string not "
            "supported for binary vectors");
  set_size(int(tokens.size()));
  for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i) {
    std::istringstream buffer(tokens[i]);
    buffer >> data[i];
    it_assert(!buffer.fail(), "Vec<bin>::set(): Stream operation failed "
              "(buffer >> data)");
  }
}

template<>
void Vec<short int>::set(const std::string &str)
{
  // parser for "short int" is the same as for "int", so reuse it here
  ivec iv(str);
  this->operator=(to_svec(iv));
}

template<>
void Vec<int>::set(const std::string &str)
{
  bool abc_format;
  std::vector<std::string> tokens = tokenize(str, abc_format);
  // no "a:b:c" tokens, so the size of output vector is known
  if (!abc_format) {
    set_size(int(tokens.size()));
    for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i)
      data[i] = parse_token(tokens[i]);
  }
  else {
    int pos = 0;
    set_size((tokens.size() > 20) ? int(tokens.size()) : 20);
    for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i) {
      // check if the current token is in "a:b[:c]" format
      if (tokens[i].find(':', 1) == std::string::npos) {
        if (++pos > datasize) {
          set_size(2 * datasize, true);
        }
        data[pos-1] = parse_token(tokens[i]);
      }
      else {
        int a, b, c;
        parse_abc_token(tokens[i], a, b, c);
        if (++pos > datasize) {
          set_size(2 * datasize, true);
        }
        data[pos-1] = a;

        if ((b > 0) && (c >= a)) {
          while ((data[pos-1] + b) <= c) {
            if (++pos > datasize) {
              set_size(2 * datasize, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if ((b < 0) && (c <= a)) {
          while ((data[pos-1] + b) >= c) {
            if (++pos > datasize) {
              set_size(2 * datasize, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b == 0 && c == a) {
          break;
        }
        else {
          it_error("Vec<int>::set(): Improper data string (a:b:c)");
        }
      }
    }
    set_size(pos, true);
  } // if (!abc_format)
}


template<>
void Vec<double>::set(const std::string &str)
{
  bool abc_format;
  std::vector<std::string> tokens = tokenize(str, abc_format);
  // no "a:b:c" tokens, so the size of output vector is known
  if (!abc_format) {
    set_size(int(tokens.size()));
    for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i)
      data[i] = parse_token(tokens[i]);
  }
  else {
    int pos = 0;
    set_size((tokens.size() > 20) ? int(tokens.size()) : 20);
    for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i) {
      // check if the current token is in "a:b[:c]" format
      if (tokens[i].find(':', 1) == std::string::npos) {
        if (++pos > datasize) {
          set_size(2 * datasize, true);
        }
        data[pos-1] = parse_token(tokens[i]);
      }
      else {
        double a, b, c;
        parse_abc_token(tokens[i], a, b, c);
        if (++pos > datasize) {
          set_size(2 * datasize, true);
        }
        data[pos-1] = a;
        // Adding this margin fixes precision problems in e.g. "0:0.2:3",
        // where the last value was 2.8 instead of 3.
        double eps_margin = std::fabs((c - a) / b) * eps;
        if ((b > 0) && (c >= a)) {
          while ((data[pos-1] + b) <= (c + eps_margin)) {
            if (++pos > datasize) {
              set_size(2 * datasize, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if ((b < 0) && (c <= a)) {
          while ((data[pos-1] + b) >= (c - eps_margin)) {
            if (++pos > datasize) {
              set_size(2 * datasize, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b == 0 && c == a) {
          break;
        }
        else {
          it_error("Vec<double>::set(): Improper data string (a:b:c)");
        }
      }
    }
    set_size(pos, true);
  } // if (!abc_format)
}

template<>
void Vec<std::complex<double> >::set(const std::string &str)
{
  bool abc_format;
  std::vector<std::string> tokens = tokenize(str, abc_format);
  it_assert(!abc_format, "Vec<bin>::set(): \"a:b:c\" format string not "
            "supported for binary vectors");
  set_size(int(tokens.size()));
  for (std::vector<std::string>::size_type i = 0; i < tokens.size(); ++i) {
    std::istringstream buffer(tokens[i]);
    buffer >> data[i];
    it_assert(!buffer.fail(), "Vec<complex>::set(): Stream operation failed "
              "(buffer >> data)");
  }
}

#if defined(HAVE_BLAS)
template<>
double dot(const vec &v1, const vec &v2)
{
  it_assert_debug(v1.length() == v2.length(), "vec::dot(): Wrong sizes");
  int incr = 1; int v1_l = v1.length();
  return blas::ddot_(&v1_l, v1._data(), &incr, v2._data(), &incr);
}
#else
template<>
double dot(const vec &v1, const vec &v2)
{
  it_assert_debug(v1.length() == v2.length(), "vec::dot(): Wrong sizes");
  double r = 0.0;
  for (int i = 0; i < v1.length(); ++i)
    r += v1._data()[i] * v2._data()[i];
  return r;
}
#endif // HAVE_BLAS

#if defined(HAVE_BLAS)
template<>
mat outer_product(const vec &v1, const vec &v2, bool)
{
  it_assert_debug((v1.length() > 0) && (v2.length() > 0),
                  "Vec::outer_product:: Input vector of zero size");
  int v1_l = v1.length(); int v2_l = v2.length();
  mat out(v1_l, v2_l);
  out.zeros();
  double alpha = 1.0;
  int incr = 1;
  blas::dger_(&v1_l, &v2_l, &alpha, v1._data(), &incr,
              v2._data(), &incr, out._data(), &v1_l);
  return out;
}

template<>
cmat outer_product(const cvec &v1, const cvec &v2, bool hermitian)
{
  it_assert_debug((v1.length() > 0) && (v2.length() > 0),
                  "Vec::outer_product:: Input vector of zero size");
  int v1_l = v1.length(); int v2_l = v2.length();
  cmat out(v1_l, v2_l);
  out.zeros();
  std::complex<double> alpha(1.0);
  int incr = 1;
  if (hermitian) {
    blas::zgerc_(&v1_l, &v2_l, &alpha, v1._data(), &incr,
                 v2._data(), &incr, out._data(), &v1_l);
  }
  else {
    blas::zgeru_(&v1_l, &v2_l, &alpha, v1._data(), &incr,
                 v2._data(), &incr, out._data(), &v1_l);
  }
  return out;
}
#else
template<>
mat outer_product(const vec &v1, const vec &v2, bool)
{
  it_assert_debug((v1.length() > 0) && (v2.length() > 0),
                  "Vec::outer_product:: Input vector of zero size");
  int v1_l = v1.length(); int v2_l = v2.length();
  mat out(v1_l, v2_l);
  for (int i = 0; i < v1_l; ++i) {
    for (int j = 0; j < v2_l; ++j) {
      out(i, j) = v1._data()[i] * v2._data()[j];
    }
  }
  return out;
}

template<>
cmat outer_product(const cvec &v1, const cvec &v2, bool hermitian)
{
  it_assert_debug((v1.length() > 0) && (v2.length() > 0),
                  "Vec::outer_product:: Input vector of zero size");
  int v1_l = v1.length(); int v2_l = v2.length();
  cmat out(v1_l, v2_l);
  if (hermitian) {
    for (int i = 0; i < v1_l; ++i) {
      for (int j = 0; j < v2_l; ++j) {
        out(i, j) = v1._data()[i] * conj(v2._data()[j]);
      }
    }
  }
  else {
    for (int i = 0; i < v1_l; ++i) {
      for (int j = 0; j < v2_l; ++j) {
        out(i, j) = v1._data()[i] * v2._data()[j];
      }
    }
  }
  return out;
}
#endif // HAVE_BLAS


template<>
bvec Vec<std::complex<double> >::operator<=(std::complex<double>) const
{
  it_error("operator<=: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator>(std::complex<double>) const
{
  it_error("operator>: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator<(std::complex<double>) const
{
  it_error("operator<: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator>=(std::complex<double>) const
{
  it_error("operator>=: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
Mat<std::complex<double> > Vec<std::complex<double> >::hermitian_transpose() const
{
  Mat<std::complex<double> > temp(1, datasize);
  for (int i = 0; i < datasize; i++)
    temp(i) = std::conj(data[i]);

  return temp;
}


//---------------------------------------------------------------------
// Instantiations
//---------------------------------------------------------------------

template class ITPP_EXPORT Vec<double>;
template class ITPP_EXPORT Vec<int>;
template class ITPP_EXPORT Vec<short int>;
template class ITPP_EXPORT Vec<std::complex<double> >;
template class ITPP_EXPORT Vec<bin>;

// addition operator

template ITPP_EXPORT vec operator+(const vec &v1, const vec &v2);
template ITPP_EXPORT cvec operator+(const cvec &v1, const cvec &v2);
template ITPP_EXPORT ivec operator+(const ivec &v1, const ivec &v2);
template ITPP_EXPORT svec operator+(const svec &v1, const svec &v2);
template ITPP_EXPORT bvec operator+(const bvec &v1, const bvec &v2);

template ITPP_EXPORT vec operator+(const vec &v1, double t);
template ITPP_EXPORT cvec operator+(const cvec &v1, std::complex<double> t);
template ITPP_EXPORT ivec operator+(const ivec &v1, int t);
template ITPP_EXPORT svec operator+(const svec &v1, short t);
template ITPP_EXPORT bvec operator+(const bvec &v1, bin t);

template ITPP_EXPORT vec operator+(double t, const vec &v1);
template ITPP_EXPORT cvec operator+(std::complex<double> t, const cvec &v1);
template ITPP_EXPORT ivec operator+(int t, const ivec &v1);
template ITPP_EXPORT svec operator+(short t, const svec &v1);
template ITPP_EXPORT bvec operator+(bin t, const bvec &v1);

// subraction operator

template ITPP_EXPORT vec operator-(const vec &v1, const vec &v2);
template ITPP_EXPORT cvec operator-(const cvec &v1, const cvec &v2);
template ITPP_EXPORT ivec operator-(const ivec &v1, const ivec &v2);
template ITPP_EXPORT svec operator-(const svec &v1, const svec &v2);
template ITPP_EXPORT bvec operator-(const bvec &v1, const bvec &v2);

template ITPP_EXPORT vec operator-(const vec &v, double t);
template ITPP_EXPORT cvec operator-(const cvec &v, std::complex<double> t);
template ITPP_EXPORT ivec operator-(const ivec &v, int t);
template ITPP_EXPORT svec operator-(const svec &v, short t);
template ITPP_EXPORT bvec operator-(const bvec &v, bin t);

template ITPP_EXPORT vec operator-(double t, const vec &v);
template ITPP_EXPORT cvec operator-(std::complex<double> t, const cvec &v);
template ITPP_EXPORT ivec operator-(int t, const ivec &v);
template ITPP_EXPORT svec operator-(short t, const svec &v);
template ITPP_EXPORT bvec operator-(bin t, const bvec &v);

// unary minus

template ITPP_EXPORT vec operator-(const vec &v);
template ITPP_EXPORT cvec operator-(const cvec &v);
template ITPP_EXPORT ivec operator-(const ivec &v);
template ITPP_EXPORT svec operator-(const svec &v);
template ITPP_EXPORT bvec operator-(const bvec &v);

// multiplication operator
template ITPP_EXPORT std::complex<double> dot(const cvec &v1, const cvec &v2);
template ITPP_EXPORT int dot(const ivec &v1, const ivec &v2);
template ITPP_EXPORT short dot(const svec &v1, const svec &v2);
template ITPP_EXPORT bin dot(const bvec &v1, const bvec &v2);

template ITPP_EXPORT double operator*(const vec &v1, const vec &v2);
template ITPP_EXPORT std::complex<double> operator*(const cvec &v1, const cvec &v2);
template ITPP_EXPORT int operator*(const ivec &v1, const ivec &v2);
template ITPP_EXPORT short operator*(const svec &v1, const svec &v2);
template ITPP_EXPORT bin operator*(const bvec &v1, const bvec &v2);

template ITPP_EXPORT imat outer_product(const ivec &v1, const ivec &v2, bool hermitian);
template ITPP_EXPORT smat outer_product(const svec &v1, const svec &v2, bool hermitian);
template ITPP_EXPORT bmat outer_product(const bvec &v1, const bvec &v2, bool hermitian);

template ITPP_EXPORT vec operator*(const vec &v, double t);
template ITPP_EXPORT cvec operator*(const cvec &v, std::complex<double> t);
template ITPP_EXPORT ivec operator*(const ivec &v, int t);
template ITPP_EXPORT svec operator*(const svec &v, short t);
template ITPP_EXPORT bvec operator*(const bvec &v, bin t);

template ITPP_EXPORT vec operator*(double t, const vec &v);
template ITPP_EXPORT cvec operator*(std::complex<double> t, const cvec &v);
template ITPP_EXPORT ivec operator*(int t, const ivec &v);
template ITPP_EXPORT svec operator*(short t, const svec &v);
template ITPP_EXPORT bvec operator*(bin t, const bvec &v);

// elementwise multiplication

template ITPP_EXPORT vec elem_mult(const vec &a, const vec &b);
template ITPP_EXPORT cvec elem_mult(const cvec &a, const cvec &b);
template ITPP_EXPORT ivec elem_mult(const ivec &a, const ivec &b);
template ITPP_EXPORT svec elem_mult(const svec &a, const svec &b);
template ITPP_EXPORT bvec elem_mult(const bvec &a, const bvec &b);

template ITPP_EXPORT void elem_mult_out(const vec &a, const vec &b, vec &out);
template ITPP_EXPORT void elem_mult_out(const cvec &a, const cvec &b, cvec &out);
template ITPP_EXPORT void elem_mult_out(const ivec &a, const ivec &b, ivec &out);
template ITPP_EXPORT void elem_mult_out(const svec &a, const svec &b, svec &out);
template ITPP_EXPORT void elem_mult_out(const bvec &a, const bvec &b, bvec &out);

template ITPP_EXPORT vec elem_mult(const vec &a, const vec &b, const vec &c);
template ITPP_EXPORT cvec elem_mult(const cvec &a, const cvec &b, const cvec &c);
template ITPP_EXPORT ivec elem_mult(const ivec &a, const ivec &b, const ivec &c);
template ITPP_EXPORT svec elem_mult(const svec &a, const svec &b, const svec &c);
template ITPP_EXPORT bvec elem_mult(const bvec &a, const bvec &b, const bvec &c);

template ITPP_EXPORT void elem_mult_out(const vec &a, const vec &b, const vec &c,
                            vec &out);
template ITPP_EXPORT void elem_mult_out(const cvec &a, const cvec &b, const cvec &c,
                            cvec &out);
template ITPP_EXPORT void elem_mult_out(const ivec &a, const ivec &b, const ivec &c,
                            ivec &out);
template ITPP_EXPORT void elem_mult_out(const svec &a, const svec &b, const svec &c,
                            svec &out);
template ITPP_EXPORT void elem_mult_out(const bvec &a, const bvec &b, const bvec &c,
                            bvec &out);

template ITPP_EXPORT vec elem_mult(const vec &a, const vec &b, const vec &c,
                       const vec &d);
template ITPP_EXPORT cvec elem_mult(const cvec &a, const cvec &b, const cvec &c,
                        const cvec &d);
template ITPP_EXPORT ivec elem_mult(const ivec &a, const ivec &b, const ivec &c,
                        const ivec &d);
template ITPP_EXPORT svec elem_mult(const svec &a, const svec &b, const svec &c,
                        const svec &d);
template ITPP_EXPORT bvec elem_mult(const bvec &a, const bvec &b, const bvec &c,
                        const bvec &d);

template ITPP_EXPORT void elem_mult_out(const vec &a, const vec &b, const vec &c,
                            const vec &d, vec &out);
template ITPP_EXPORT void elem_mult_out(const cvec &a, const cvec &b, const cvec &c,
                            const cvec &d, cvec &out);
template ITPP_EXPORT void elem_mult_out(const ivec &a, const ivec &b, const ivec &c,
                            const ivec &d, ivec &out);
template ITPP_EXPORT void elem_mult_out(const svec &a, const svec &b, const svec &c,
                            const svec &d, svec &out);
template ITPP_EXPORT void elem_mult_out(const bvec &a, const bvec &b, const bvec &c,
                            const bvec &d, bvec &out);

// in-place elementwise multiplication

template ITPP_EXPORT void elem_mult_inplace(const vec &a, vec &b);
template ITPP_EXPORT void elem_mult_inplace(const cvec &a, cvec &b);
template ITPP_EXPORT void elem_mult_inplace(const ivec &a, ivec &b);
template ITPP_EXPORT void elem_mult_inplace(const svec &a, svec &b);
template ITPP_EXPORT void elem_mult_inplace(const bvec &a, bvec &b);

// elementwise multiplication followed by summation

template ITPP_EXPORT double elem_mult_sum(const vec &a, const vec &b);
template ITPP_EXPORT std::complex<double> elem_mult_sum(const cvec &a, const cvec &b);
template ITPP_EXPORT int elem_mult_sum(const ivec &a, const ivec &b);
template ITPP_EXPORT short elem_mult_sum(const svec &a, const svec &b);
template ITPP_EXPORT bin elem_mult_sum(const bvec &a, const bvec &b);

// division operator

template ITPP_EXPORT vec operator/(const vec &v, double t);
template ITPP_EXPORT cvec operator/(const cvec &v, std::complex<double> t);
template ITPP_EXPORT ivec operator/(const ivec &v, int t);
template ITPP_EXPORT svec operator/(const svec &v, short t);
template ITPP_EXPORT bvec operator/(const bvec &v, bin t);

template ITPP_EXPORT vec operator/(double t, const vec &v);
template ITPP_EXPORT cvec operator/(std::complex<double> t, const cvec &v);
template ITPP_EXPORT ivec operator/(int t, const ivec &v);
template ITPP_EXPORT svec operator/(short t, const svec &v);
template ITPP_EXPORT bvec operator/(bin t, const bvec &v);

// elementwise division operator

template ITPP_EXPORT vec elem_div(const vec &a, const vec &b);
template ITPP_EXPORT cvec elem_div(const cvec &a, const cvec &b);
template ITPP_EXPORT ivec elem_div(const ivec &a, const ivec &b);
template ITPP_EXPORT svec elem_div(const svec &a, const svec &b);
template ITPP_EXPORT bvec elem_div(const bvec &a, const bvec &b);

template ITPP_EXPORT vec elem_div(double t, const vec &v);
template ITPP_EXPORT cvec elem_div(std::complex<double> t, const cvec &v);
template ITPP_EXPORT ivec elem_div(int t, const ivec &v);
template ITPP_EXPORT svec elem_div(short t, const svec &v);
template ITPP_EXPORT bvec elem_div(bin t, const bvec &v);

template ITPP_EXPORT void elem_div_out(const vec &a, const vec &b, vec &out);
template ITPP_EXPORT void elem_div_out(const cvec &a, const cvec &b, cvec &out);
template ITPP_EXPORT void elem_div_out(const ivec &a, const ivec &b, ivec &out);
template ITPP_EXPORT void elem_div_out(const svec &a, const svec &b, svec &out);
template ITPP_EXPORT void elem_div_out(const bvec &a, const bvec &b, bvec &out);

// elementwise division followed by summation

template ITPP_EXPORT double elem_div_sum(const vec &a, const vec &b);
template ITPP_EXPORT std::complex<double> elem_div_sum(const cvec &a, const cvec &b);
template ITPP_EXPORT int elem_div_sum(const ivec &a, const ivec &b);
template ITPP_EXPORT short elem_div_sum(const svec &a, const svec &b);
template ITPP_EXPORT bin elem_div_sum(const bvec &a, const bvec &b);

// concat operator

template ITPP_EXPORT vec concat(const vec &v, double a);
template ITPP_EXPORT cvec concat(const cvec &v, std::complex<double> a);
template ITPP_EXPORT ivec concat(const ivec &v, int a);
template ITPP_EXPORT svec concat(const svec &v, short a);
template ITPP_EXPORT bvec concat(const bvec &v, bin a);

template ITPP_EXPORT vec concat(double a, const vec &v);
template ITPP_EXPORT cvec concat(std::complex<double> a, const cvec &v);
template ITPP_EXPORT ivec concat(int a, const ivec &v);
template ITPP_EXPORT svec concat(short a, const svec &v);
template ITPP_EXPORT bvec concat(bin a, const bvec &v);

template ITPP_EXPORT vec concat(const vec &v1, const vec &v2);
template ITPP_EXPORT cvec concat(const cvec &v1, const cvec &v2);
template ITPP_EXPORT ivec concat(const ivec &v1, const ivec &v2);
template ITPP_EXPORT svec concat(const svec &v1, const svec &v2);
template ITPP_EXPORT bvec concat(const bvec &v1, const bvec &v2);

template ITPP_EXPORT vec concat(const vec &v1, const vec &v2, const vec &v3);
template ITPP_EXPORT cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
template ITPP_EXPORT ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
template ITPP_EXPORT svec concat(const svec &v1, const svec &v2, const svec &v3);
template ITPP_EXPORT bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

template ITPP_EXPORT vec concat(const vec &v1, const vec &v2,
                    const vec &v3, const vec &v4);
template ITPP_EXPORT cvec concat(const cvec &v1, const cvec &v2,
                     const cvec &v3, const cvec &v4);
template ITPP_EXPORT ivec concat(const ivec &v1, const ivec &v2,
                     const ivec &v3, const ivec &v4);
template ITPP_EXPORT svec concat(const svec &v1, const svec &v2,
                     const svec &v3, const svec &v4);
template ITPP_EXPORT bvec concat(const bvec &v1, const bvec &v2,
                     const bvec &v3, const bvec &v4);

template ITPP_EXPORT vec concat(const vec &v1, const vec &v2, const vec &v3,
                    const vec &v4, const vec &v5);
template ITPP_EXPORT cvec concat(const cvec &v1, const cvec &v2, const cvec &v3,
                     const cvec &v4, const cvec &v5);
template ITPP_EXPORT ivec concat(const ivec &v1, const ivec &v2, const ivec &v3,
                     const ivec &v4, const ivec &v5);
template ITPP_EXPORT svec concat(const svec &v1, const svec &v2, const svec &v3,
                     const svec &v4, const svec &v5);
template ITPP_EXPORT bvec concat(const bvec &v1, const bvec &v2, const bvec &v3,
                     const bvec &v4, const bvec &v5);

// I/O streams

template ITPP_EXPORT std::ostream &operator<<(std::ostream& os, const vec &vect);
template ITPP_EXPORT std::ostream &operator<<(std::ostream& os, const cvec &vect);
template ITPP_EXPORT std::ostream &operator<<(std::ostream& os, const svec &vect);
template ITPP_EXPORT std::ostream &operator<<(std::ostream& os, const ivec &vect);
template ITPP_EXPORT std::ostream &operator<<(std::ostream& os, const bvec &vect);
template ITPP_EXPORT std::istream &operator>>(std::istream& is, vec &vect);
template ITPP_EXPORT std::istream &operator>>(std::istream& is, cvec &vect);
template ITPP_EXPORT std::istream &operator>>(std::istream& is, svec &vect);
template ITPP_EXPORT std::istream &operator>>(std::istream& is, ivec &vect);
template ITPP_EXPORT std::istream &operator>>(std::istream& is, bvec &vect);

} // namespace itpp

//! \endcond
