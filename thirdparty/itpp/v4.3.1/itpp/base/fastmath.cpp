/*!
 * \file
 * \brief Implementation of special operations on vectors and matricies
 * optimized for speed
 * \author Tony Ottosson and Tobias Ringstrom
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

#include <itpp/base/fastmath.h>


namespace itpp
{

// m=m-v*v'*m
void sub_v_vT_m(mat &m, const vec &v)
{
  vec v2(m.cols());
  double tmp, *v2p;
  const double *vp;
  int i, j;

  it_assert(v.size() == m.rows(), "sub_v_vT_m()");

  v2p = v2._data();
  for (j = 0; j < m.cols(); j++) {
    tmp = 0.0;
    vp = v._data();
    for (i = 0; i < m.rows(); i++)
      tmp += *(vp++) * m._elem(i, j);
    *(v2p++) = tmp;
  }

  vp = v._data();
  for (i = 0; i < m.rows(); i++) {
    v2p = v2._data();
    for (j = 0; j < m.cols(); j++)
      m._elem(i, j) -= *vp * *(v2p++);
    vp++;
  }
}

// m=m-m*v*v'
void sub_m_v_vT(mat &m, const vec &v)
{
  vec v2(m.rows());
  double tmp, *v2p;
  const double *vp;
  int i, j;

  it_assert(v.size() == m.cols(), "sub_m_v_vT()");

  v2p = v2._data();
  for (i = 0; i < m.rows(); i++) {
    tmp = 0.0;
    vp = v._data();
    for (j = 0; j < m.cols(); j++)
      tmp += *(vp++) * m._elem(i, j);
    *(v2p++) = tmp;
  }

  v2p = v2._data();
  for (i = 0; i < m.rows(); i++) {
    vp = v._data();
    for (j = 0; j < m.cols(); j++)
      m._elem(i, j) -= *v2p * *(vp++);
    v2p++;
  }
}

} // namespace itpp
