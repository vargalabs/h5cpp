/*!
 * \file
 * \brief Implementation of converters between different vector and matrix types
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger and Adam Piatyszek
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

#include <itpp/base/converters.h>
#include <itpp/base/itcompat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/math/log_exp.h>

//! \cond

namespace itpp
{

// ----------------------------------------------------------------------
// Vector converters
// ----------------------------------------------------------------------

ivec to_ivec(int s) { ivec out(1); out(0) = s; return out; }

vec to_vec(double s) { vec out(1); out(0) = s; return out; }

cvec to_cvec(double real, double imag)
{
  cvec out(1);
  out(0) = std::complex<double>(real, imag);
  return out;
}

// ----------------------------------------------------------------------
// Miscellaneous converters
// ----------------------------------------------------------------------

bvec dec2bin(int length, int index)
{
  int i, bintemp = index;
  bvec temp(length);

  for (i = length - 1; i >= 0; i--) {
    temp(i) = bin(bintemp & 1);
    bintemp = (bintemp >> 1);
  }
  return temp;
}

bvec dec2bin(int index, bool msb_first)
{
  int length = int2bits(index);
  int i, bintemp = index;
  bvec temp(length);

  for (i = length - 1; i >= 0; i--) {
    temp(i) = bin(bintemp & 1);
    bintemp = (bintemp >> 1);
  }
  if (msb_first) {
    return temp;
  }
  else {
    return reverse(temp);
  }
}

void dec2bin(int index, bvec &v)
{
  int i, bintemp = index;
  v.set_size(int2bits(index), false);

  for (i = v.size() - 1; i >= 0; i--) {
    v(i) = bin(bintemp & 1);
    bintemp = (bintemp >> 1);
  }
}

int bin2dec(const bvec &inbvec, bool msb_first)
{
  int i, temp = 0;
  int sizebvec = inbvec.length();
  if (msb_first) {
    for (i = 0; i < sizebvec; i++) {
      temp += pow2i(sizebvec - i - 1) * int(inbvec(i));
    }
  }
  else {
    for (i = 0; i < sizebvec; i++) {
      temp += pow2i(i) * int(inbvec(i));
    }
  }
  return temp;
}

bvec oct2bin(const ivec &octalindex, short keepzeros)
{
  int length = octalindex.length(), i;
  bvec out(3*length);
  for (i = 0; i < length; i++) {
    out.replace_mid(3*i, dec2bin(3, octalindex(i)));
  }
  //remove zeros if keepzeros = 0
  if (keepzeros == 0) {
    for (i = 0; i < out.length(); i++) {
      if ((short)out(i) != 0) {
        return out.right(out.length() - i);
        break;
      }
    }
    return bvec("0");
  }
  else {
    return out;
  }
}

ivec bin2oct(const bvec &inbits)
{
  int start, Itterations = ceil_i(inbits.length() / 3.0);
  ivec out(Itterations);
  for (int i = Itterations - 1; i > 0; i--) {
    start = 3 * i - (3 * Itterations - inbits.length());
    out(i) = bin2dec(inbits.mid(start, 3));
  }
  out(0) = bin2dec(inbits.left(inbits.length() - ((Itterations - 1) * 3)));
  return out;
}

ivec bin2pol(const bvec &inbvec)
{
  return 1 -2*to_ivec(inbvec);
}

bvec pol2bin(const ivec &inpol)
{
  return to_bvec((1 -inpol) / 2);
}


// Round to nearest integer, return result in double
double round(double x) { return ::rint(x); }
// Round to nearest integer
vec round(const vec &x) { return apply_function<double>(::rint, x); }
// Round to nearest integer
mat round(const mat &x) { return apply_function<double>(::rint, x); }
// Round to nearest integer
int round_i(double x) { return static_cast<int>(::rint(x)); }

// Round to nearest integer and return ivec
ivec round_i(const vec &x) { return to_ivec(round(x)); }
// Round to nearest integer and return imat
imat round_i(const mat &x) { return to_imat(round(x)); }

// Round to nearest upper integer
ivec ceil_i(const vec &x) { return to_ivec(ceil(x)); }
// Round to nearest upper integer
imat ceil_i(const mat &x) { return to_imat(ceil(x)); }

// Round to nearest lower integer
ivec floor_i(const vec &x) { return to_ivec(floor(x)); }
// Round to nearest lower integer
imat floor_i(const mat &x) { return to_imat(floor(x)); }


cvec round_to_zero(const cvec &x, double threshold)
{
  cvec temp(x.length());

  for (int i = 0; i < x.length(); i++)
    temp(i) = round_to_zero(x(i), threshold);

  return temp;
}

cmat round_to_zero(const cmat &x, double threshold)
{
  cmat temp(x.rows(), x.cols());

  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      temp(i, j) = round_to_zero(x(i, j), threshold);
    }
  }

  return temp;
}

cvec round_to_infty(const cvec &in, const double threshold)
{
  cvec temp(in.length());

  for (int i = 0; i < in.length(); i++)
    temp(i) = round_to_infty(in(i), threshold);

  return temp;
}

cmat round_to_infty(const cmat &in, const double threshold)
{
  cmat temp(in.rows(), in.cols());

  for (int i = 0; i < in.rows(); i++) {
    for (int j = 0; j < in.cols(); j++) {
      temp(i, j) = round_to_infty(in(i, j), threshold);
    }
  }

  return temp;
}

std::string to_str(const double &i, const int precision)
{
  std::ostringstream ss;
  ss.precision(precision);
  ss.setf(std::ostringstream::scientific, std::ostringstream::floatfield);
  ss << i;
  return ss.str();
}

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template ITPP_EXPORT bvec to_bvec(const svec &v);
template ITPP_EXPORT bvec to_bvec(const ivec &v);

template ITPP_EXPORT svec to_svec(const bvec &v);
template ITPP_EXPORT svec to_svec(const ivec &v);
template ITPP_EXPORT svec to_svec(const vec &v);

template ITPP_EXPORT ivec to_ivec(const bvec &v);
template ITPP_EXPORT ivec to_ivec(const svec &v);
template ITPP_EXPORT ivec to_ivec(const vec &v);

template ITPP_EXPORT vec to_vec(const bvec &v);
template ITPP_EXPORT vec to_vec(const svec &v);
template ITPP_EXPORT vec to_vec(const ivec &v);

template ITPP_EXPORT cvec to_cvec(const bvec &v);
template ITPP_EXPORT cvec to_cvec(const svec &v);
template ITPP_EXPORT cvec to_cvec(const ivec &v);
template ITPP_EXPORT cvec to_cvec(const vec &v);

template ITPP_EXPORT cvec to_cvec(const bvec &real, const bvec &imag);
template ITPP_EXPORT cvec to_cvec(const svec &real, const svec &imag);
template ITPP_EXPORT cvec to_cvec(const ivec &real, const ivec &imag);
template ITPP_EXPORT cvec to_cvec(const vec &real, const vec &imag);

template ITPP_EXPORT bmat to_bmat(const smat &m);
template ITPP_EXPORT bmat to_bmat(const imat &m);

template ITPP_EXPORT smat to_smat(const bmat &m);
template ITPP_EXPORT smat to_smat(const imat &m);
template ITPP_EXPORT smat to_smat(const mat &m);

template ITPP_EXPORT imat to_imat(const bmat &m);
template ITPP_EXPORT imat to_imat(const smat &m);
template ITPP_EXPORT imat to_imat(const mat &m);

template ITPP_EXPORT mat to_mat(const bmat &m);
template ITPP_EXPORT mat to_mat(const smat &m);
template ITPP_EXPORT mat to_mat(const imat &m);

template ITPP_EXPORT cmat to_cmat(const bmat &m);
template ITPP_EXPORT cmat to_cmat(const smat &m);
template ITPP_EXPORT cmat to_cmat(const imat &m);
template ITPP_EXPORT cmat to_cmat(const mat &m);

template ITPP_EXPORT cmat to_cmat(const bmat &real, const bmat &imag);
template ITPP_EXPORT cmat to_cmat(const smat &real, const smat &imag);
template ITPP_EXPORT cmat to_cmat(const imat &real, const imat &imag);
template ITPP_EXPORT cmat to_cmat(const mat &real, const mat &imag);

} // namespace itpp

//! \endcond
