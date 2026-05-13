/*!
 * \file
 * \brief Implementation of Galois Field algebra classes and functions
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

#include <itpp/comm/galois.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/itcompat.h>
#include <string>
#include <iostream>


namespace itpp
{

Array<Array<int> > GF::alphapow;
Array<Array<int> > GF::logalpha;
ivec GF::q = "1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536";

// set q=2^mvalue
void GF::set_size(int qvalue)
{
  m = static_cast<char>(round_i(::log2(static_cast<double>(qvalue))));
  it_assert((1 << m) == qvalue, "GF::setsize : q is not a power of 2");
  it_assert((m > 0) && (m <= 16), "GF::setsize : q must be positive and "
                                  "less than or equal to 2^16");

  /* Construct GF(q), q=2^m. From Wicker, "Error Control Systems
     for digital communication and storage" pp. 463-465 */

  int reduce, temp, n;
  const int reducetable[] = {3, 3, 3, 5, 3, 9, 29, 17, 9, 5, 83, 27, 43, 3, 4107}; // starts at m=2,..,16

  if (alphapow.size() < (m + 1)) {
    alphapow.set_size(m + 1, true);
    logalpha.set_size(m + 1, true);
  }

  if (alphapow(m).size() == 0) {
    alphapow(m).set_size(qvalue);
    logalpha(m).set_size(qvalue);
    alphapow(m) = 0;
    logalpha(m) = 0;
    if (m == 1) { // GF(2), special case
      alphapow(1)(0) = 1;
      logalpha(1)(0) = -1;
      logalpha(1)(1) = 0;
    }
    else {
      reduce = reducetable[m-2];
      alphapow(m)(0) = 1; // alpha^0 = 1
      for (n = 1; n < (1 << m) - 1; n++) {
        temp = alphapow(m)(n - 1);
        temp = (temp << 1); // multiply by alpha
        if (temp & (1 << m)) // contains alpha**m term
          alphapow(m)(n) = (temp & ~(1 << m)) ^ reduce;
        else
          alphapow(m)(n) = temp; // if no alpha**m term, store as is

        // create table to go in opposite direction
        logalpha(m)(0) = -1; // special case, actually log(0)=-inf
      }

      for (n = 0;n < (1 << m) - 1;n++)
        logalpha(m)(alphapow(m)(n)) = n;
    }
  }
}
//! Input stream operator for GF
std::istream &operator>>(std::istream &is, GF &ingf)
{
  int val; char c;
  static const std::string prefix("alpha^");
  c = is.get();
  if(c == 'a') {
  //read alpha^pow form from stream
    std::string::const_iterator pr_it = prefix.begin(); pr_it++;
    for(; pr_it < prefix.end(); ++pr_it) {
        c = is.get();
        if(*pr_it != c) {
          is.setstate(std::ios_base::failbit);
          return is;
        }
    }
    is >> val;
    if(is) ingf.set(ingf.get_size(),val);
  }
  else {
  //try to read 0 from stream
    is >> val;
    if(is && (val==0)) {
      ingf.set(ingf.get_size(),0);
    }
    else {
      is.setstate(std::ios_base::failbit);
    }
  }
  return is;
}

//! Output stream operator for GF
std::ostream &operator<<(std::ostream &os, const GF &ingf)
{
  if (ingf.value == -1)
    os << "0";
  else
    os << "alpha^" << ingf.value;
  return os;
}

//! Output stream operator for GFX
std::ostream &operator<<(std::ostream &os, const GFX &ingfx)
{
  int terms = 0;
  for (int i = 0; i < ingfx.degree + 1; i++) {
    if (ingfx.coeffs(i) != GF(ingfx.q, -1)) {
      if (terms != 0) os << " + ";
      terms++;
      if (ingfx.coeffs(i) == GF(ingfx.q, 0)) {// is the coefficient an one (=alpha^0=1)
        os  << "x^" << i;
      }
      else {
        os  << ingfx.coeffs(i) << "*x^" << i;
      }
    }
  }
  if (terms == 0) os << "0";
  return os;
}

//----------------- Help Functions -----------------

//! Division of two GFX (local help function)
GFX divgfx(const GFX &c, const GFX &g)
{
  int q = c.get_size();
  GFX temp = c;
  int tempdegree = temp.get_true_degree();
  int gdegree = g.get_true_degree();
  int degreedif = tempdegree - gdegree;
  if (degreedif < 0) return GFX(q, 0); // denominator larger than nominator. Return zero polynomial.
  GFX m(q, degreedif), divisor(q);

  for (int i = 0; i < c.get_degree(); i++) {
    m[degreedif] = temp[tempdegree] / g[gdegree];
    divisor.set_degree(degreedif);
    divisor.clear();
    divisor[degreedif] = m[degreedif];
    temp -= divisor * g;
    tempdegree = temp.get_true_degree();
    degreedif = tempdegree - gdegree;
    if ((degreedif < 0) || (temp.get_true_degree() == 0 && temp[0] == GF(q, -1))) {
      break;
    }
  }
  return m;
}

//! Modulo function of two GFX (local help function)
GFX modgfx(const GFX &a, const GFX &b)
{
  int q = a.get_size();
  GFX temp = a;
  int tempdegree = temp.get_true_degree();
  int bdegree = b.get_true_degree();
  int degreedif = a.get_true_degree() - b.get_true_degree();
  if (degreedif < 0) return temp; // Denominator larger than nominator. Return nominator.
  GFX m(q, degreedif), divisor(q);

  for (int i = 0; i < a.get_degree(); i++) {
    m[degreedif] = temp[tempdegree] / b[bdegree];
    divisor.set_degree(degreedif);
    divisor.clear();
    divisor[degreedif] =  m[degreedif];
    temp -= divisor * b; // Bug-fixed. Used to be: temp -= divisor*a;
    tempdegree = temp.get_true_degree();
    degreedif = temp.get_true_degree() - bdegree;
    if ((degreedif < 0) || (temp.get_true_degree() == 0 && temp[0] == GF(q, -1))) {
      break;
    }
  }
  return temp;
}

} // namespace itpp
