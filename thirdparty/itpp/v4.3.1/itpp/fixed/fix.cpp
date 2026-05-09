/*!
 * \file
 * \brief Implementation of a fixed-point data type Fix
 * \author Johan Bergman
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

#include <itpp/fixed/fix.h>
#include <itpp/base/itassert.h>
#include <cstdio>
#include <iostream>


namespace itpp
{

Fix& Fix::operator=(const Fix &x)
{
  shift = x.shift;
  re = apply_o_mode(x.re);
  return *this;
}

Fix& Fix::operator=(const int x)
{
  shift = 0;
  re = apply_o_mode(x);
  return *this;
}

Fix& Fix::operator+=(const Fix &x)
{
  shift = assert_shifts(*this, x);
  re = apply_o_mode(re + x.re);
  return *this;
}

Fix& Fix::operator+=(const int x)
{
  assert_shifts(*this, x);
  re = apply_o_mode(re + x);
  return *this;
}

Fix& Fix::operator-=(const Fix &x)
{
  shift = assert_shifts(*this, x);
  re = apply_o_mode(re - x.re);
  return *this;
}

Fix& Fix::operator-=(const int x)
{
  assert_shifts(*this, x);
  re = apply_o_mode(re - x);
  return *this;
}

Fix& Fix::operator*=(const Fix &x)
{
  shift += x.shift;
  re = apply_o_mode(re * x.re);
  return *this;
}

Fix& Fix::operator*=(const int x)
{
  re = apply_o_mode(re * x);
  return *this;
}

Fix& Fix::operator/=(const Fix &x)
{
  shift -= x.shift;
  re = apply_o_mode(re / x.re);
  return *this;
}

Fix& Fix::operator/=(const int x)
{
  re = apply_o_mode(re / x);
  return *this;
}

Fix Fix::operator-() const
{
  return Fix(-re, shift, 0, 0);
}

Fix& Fix::operator<<=(const int n)
{
  it_assert_debug(n >= 0, "Fix::operator<<=: n cannot be negative!");
  shift += n;
  re = apply_o_mode(re << n);
  return *this;
}

Fix& Fix::operator>>=(const int n)
{
  shift -= n;
  re = rshift_and_apply_q_mode(re, n);
  return *this;
}

void Fix::set(double x, int n)
{
  shift = n;
  re = scale_and_apply_modes(x);
}

void Fix::set(double x, int n, q_mode q)
{
  shift = n;
  re = scale_and_apply_modes(x, q);
}

void Fix::lshift(int n)
{
  it_assert_debug(n >= 0, "Fix::lshift: n cannot be negative!");
  shift += n;
  re = apply_o_mode(re << n);
}

void Fix::rshift(int n)
{
  shift -= n;
  re = rshift_and_apply_q_mode(re, n);
}

void Fix::rshift(int n, q_mode q)
{
  shift -= n;
  re = rshift_and_apply_q_mode(re, n, q);
}

double Fix::unfix() const
{
  it_assert_debug(shift >= -63 && shift <= 64, "Fix::unfix: Illegal shift!");
  return double(re)*DOUBLE_POW2[64 - shift];
}

void Fix::print() const
{
  Fix_Base::print();
  std::cout << "re = " << re << std::endl;
}

int assert_shifts(const Fix &x, const Fix &y)
{
  int ret = 0;

  if (x.shift == y.shift)
    ret = x.shift;
  else if (x.re == 0)
    ret = y.shift;
  else if (y.re == 0)
    ret = x.shift;
  else
    it_error("assert_shifts: Different shifts not allowed!");

  return ret;
}

int assert_shifts(const Fix &x, int y)
{
  it_error_if((x.shift != 0) && (x.re != 0) && (y != 0),
              "assert_shifts: Different shifts not allowed!");
  return x.shift;
}

std::istream &operator>>(std::istream &is, Fix &x)
{
  double value;
  is >> value;
  if (!is.eof() && (is.peek() == '<')) {
    int shift;
    is.get();  // Swallow '<' sign
    if (is.peek() == '<') {
      is.get();  // Swallow '<' sign
      is >> shift;
      x.set(value, shift);
    }
    else {
      is >> shift;
      is.get();  // Swallow '>' sign
      x.set_re(fixrep(value));
      x.set_shift(shift);
    }
  }
  else {
    // Change data representation but keep shift
    x.set_re(fixrep(value));
  }
  return is;
}

std::ostream &operator<<(std::ostream &os, const Fix &x)
{
  switch (x.get_output_mode()) {
  case OUTPUT_FIX:
    os << x.get_re();
    break;
  case OUTPUT_FIX_SHIFT:
    os << x.get_re() << '<' << x.get_shift() << '>';
    break;
  case OUTPUT_FLOAT:
    os << double(x);
    break;
  case OUTPUT_FLOAT_SHIFT:
    os << double(x) << "<<" << x.get_shift();
    break;
  default:
    it_error("operator<<: Illegal output mode!");
  }
  return os;
}

// Specialization of template definition in vec.cpp
template<>
void fixvec::set(const char *values)
{
  std::istringstream buffer(values);
  int b = 0, c = 0;
  int default_shift = 0, pos = 0, maxpos = 10;
  if (datasize > 0) {
    // Assume that all elements have the same shift
    default_shift = data[0].get_shift();
  }
  alloc(maxpos);
  while (buffer.peek() != EOF) {
    switch (buffer.peek()) {
    case ':': // reads format a:b:c or a:b
      buffer.get();
      if (!buffer.eof()) {
        buffer >> b;
      }
      if (!buffer.eof() && buffer.peek() == ':') {
        buffer.get();
        if (!buffer.eof()) {
          buffer >> c;
          while (int(double(data[pos-1])) + b - c <= 0) {
            pos++;
            if (pos > maxpos) {
              maxpos = maxpos * 2;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2];
            data[pos-1] += b;
          }
        }
      }
      else {
        while (int(double(data[pos-1])) < b) {
          pos++;
          if (pos > maxpos) {
            maxpos = maxpos * 2;
            set_size(maxpos, true);
          }
          data[pos-1] = data[pos-2];
          data[pos-1] += 1;
        }
      }
      break;
    case ',':
      buffer.get();
      break;
    default:
      pos++;
      if (pos > maxpos) {
        maxpos *= 2;
        set_size(maxpos, true);
      }
      data[pos-1].set_shift(default_shift);
      buffer >> data[pos-1];  // May override default_shift
      while (buffer.peek() == ' ') { buffer.get(); }
      break;
    }
  }
  set_size(pos, true);
}

// Specialization of template definition in mat.cpp
template<>
void fixmat::set(const char *values)
{
  std::istringstream buffer(values);
  int default_shift = 0, rows = 0, maxrows = 10, cols = 0, nocols = 0, maxcols = 10;
  if (datasize > 0) {
    // Assume that all elements have the same shift
    default_shift = data[0].get_shift();
  }
  alloc(maxrows, maxcols);
  while (buffer.peek() != EOF) {
    rows++;
    if (rows > maxrows) {
      maxrows = maxrows * 2;
      set_size(maxrows, maxcols, true);
    }
    cols = 0;
    while ((buffer.peek() != ';') && (buffer.peek() != EOF)) {
      if (buffer.peek() == ',') {
        buffer.get();
      }
      else {
        cols++;
        if (cols > nocols) {
          nocols = cols;
          if (cols > maxcols) {
            maxcols = maxcols * 2;
            set_size(maxrows, maxcols, true);
          }
        }
        this->operator()(rows-1, cols - 1).set_shift(default_shift);
        buffer >> this->operator()(rows-1, cols - 1);  // May override default_shift
        while (buffer.peek() == ' ') { buffer.get(); }
      }
    }
    if (!buffer.eof())
      buffer.get();
  }
  set_size(rows, nocols, true);
}

} //namespace itpp
