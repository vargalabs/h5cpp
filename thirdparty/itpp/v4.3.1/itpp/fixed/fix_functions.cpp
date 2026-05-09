/*!
 * \file
 * \brief Implementation of a set of functions for Fix, Fixed, CFix and
 * CFixed classes
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

#include <itpp/fixed/fix_functions.h>


namespace itpp
{

vec to_vec(const fixvec &v)
{
  vec temp(v.length());
  for (int i = 0; i < v.length(); i++) {
    temp(i) = v(i).unfix();
  }
  return temp;
}

cvec to_cvec(const cfixvec &v)
{
  cvec temp(v.length());
  for (int i = 0; i < v.length(); i++) {
    temp(i) = v(i).unfix();
  }
  return temp;
}

mat to_mat(const fixmat &m)
{
  mat temp(m.rows(), m.cols());
  for (int i = 0; i < m.rows(); i++) {
    for (int j = 0; j < m.cols(); j++) {
      temp(i, j) = m(i, j).unfix();
    }
  }
  return temp;
}

cmat to_cmat(const cfixmat &m)
{
  cmat temp(m.rows(), m.cols());
  for (int i = 0; i < m.rows(); i++) {
    for (int j = 0; j < m.cols(); j++) {
      temp(i, j) = m(i, j).unfix();
    }
  }
  return temp;
}

Fix abs(const Fix &x)
{
  fixrep tmp = x.get_re();
  return Fix((tmp >= 0 ? tmp : -tmp),  // Risk for overflow!
             x.get_shift(),
             0, 0);
}

Fix real(const CFix &x)
{
  return Fix(x.get_re(),
             x.get_shift(),
             0, 0);
}

Fix imag(const CFix &x)
{
  return Fix(x.get_im(),
             x.get_shift(),
             0, 0);
}

CFix conj(const CFix &x)
{
  return CFix(x.get_re(),
              -x.get_im(),
              x.get_shift(),
              0, 0);
}

} // namespace itpp
