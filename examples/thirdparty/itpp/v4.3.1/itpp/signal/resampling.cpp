/*!
 * \file
 * \brief Resampling functions - source file
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/signal/resampling.h>


namespace itpp
{

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template ITPP_EXPORT vec repeat(const vec &v, int norepeats);
template ITPP_EXPORT cvec repeat(const cvec &v, int norepeats);
template ITPP_EXPORT svec repeat(const svec &v, int norepeats);
template ITPP_EXPORT ivec repeat(const ivec &v, int norepeats);
template ITPP_EXPORT bvec repeat(const bvec &v, int norepeats);

template ITPP_EXPORT mat repeat(const mat &m, int norepeats);
template ITPP_EXPORT cmat repeat(const cmat &m, int norepeats);
template ITPP_EXPORT smat repeat(const smat &m, int norepeats);
template ITPP_EXPORT imat repeat(const imat &m, int norepeats);
template ITPP_EXPORT bmat repeat(const bmat &m, int norepeats);

template ITPP_EXPORT vec upsample(const vec &v, int usf);
template ITPP_EXPORT cvec upsample(const cvec &v, int usf);
template ITPP_EXPORT svec upsample(const svec &v, int usf);
template ITPP_EXPORT ivec upsample(const ivec &v, int usf);
template ITPP_EXPORT bvec upsample(const bvec &v, int usf);

template ITPP_EXPORT mat upsample(const mat &v, int usf);
template ITPP_EXPORT cmat upsample(const cmat &v, int usf);
template ITPP_EXPORT smat upsample(const smat &v, int usf);
template ITPP_EXPORT imat upsample(const imat &v, int usf);
template ITPP_EXPORT bmat upsample(const bmat &v, int usf);

template ITPP_EXPORT void upsample(const vec &v, int usf,  vec & u);
template ITPP_EXPORT void upsample(const cvec &v, int usf,  cvec & u);
template ITPP_EXPORT void upsample(const svec &v, int usf,  svec & u);
template ITPP_EXPORT void upsample(const ivec &v, int usf,  ivec & u);
template ITPP_EXPORT void upsample(const bvec &v, int usf,  bvec & u);

template ITPP_EXPORT void upsample(const mat &v, int usf,  mat & u);
template ITPP_EXPORT void upsample(const cmat &v, int usf,  cmat & u);
template ITPP_EXPORT void upsample(const smat &v, int usf,  smat & u);
template ITPP_EXPORT void upsample(const imat &v, int usf,  imat & u);
template ITPP_EXPORT void upsample(const bmat &v, int usf,  bmat & u);

template ITPP_EXPORT vec lininterp(const vec &v, int usf);
template ITPP_EXPORT cvec lininterp(const cvec &v, int usf);

template ITPP_EXPORT mat lininterp(const mat &v, int usf);
template ITPP_EXPORT cmat lininterp(const cmat &v, int usf);

template ITPP_EXPORT void lininterp(const vec &v, int usf,  vec & u);
template ITPP_EXPORT void lininterp(const cvec &v, int usf,  cvec & u);

template ITPP_EXPORT void lininterp(const mat &v, int usf,  mat & u);
template ITPP_EXPORT void lininterp(const cmat &v, int usf,  cmat & u);

template ITPP_EXPORT mat lininterp(const mat &m, double f_base, double f_ups, int nrof_samples, double t_start);
template ITPP_EXPORT cmat lininterp(const cmat &m, double f_base, double f_ups, int nrof_samples, double t_start);

template ITPP_EXPORT vec lininterp(const vec &v, double f_base, double f_ups, int nrof_samples, double t_start);
template ITPP_EXPORT cvec lininterp(const cvec &v, double f_base, double f_ups, int nrof_samples, double t_start);

} // namespace itpp
