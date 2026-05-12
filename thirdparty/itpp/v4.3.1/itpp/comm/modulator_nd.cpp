/*!
 * \file
 * \brief Implementation of vector (MIMO) modulator classes
 * \author Mirsad Cirkic, Erik G. Larsson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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

#include <itpp/comm/modulator_nd.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/algebra/cholesky.h>
#include <itpp/base/algebra/inv.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/converters.h>
#include <itpp/base/itcompat.h>
#include <itpp/base/sort.h>
#include <itpp/stat/misc_stat.h>

#include <cmath>
#include <iostream>
#include <iomanip>

namespace itpp
{

// ----------------------------------------------------------------------
// Modulator_ND
// ----------------------------------------------------------------------


QLLRvec Modulator_ND::probabilities(QLLR l)
{
  QLLRvec result(2);

  if(l < 0) {         // this can be done more efficiently
    result(1) = -llrcalc.jaclog(0, -l);
    result(0) = result(1) - l;
  }
  else {
    result(0) = -llrcalc.jaclog(0, l);
    result(1) = result(0) + l;
  }
  return result;
}

	void Modulator_ND::update_LLR(const Array<QLLRvec> &logP_apriori, int s,
                              QLLR scaled_norm, int j, QLLRvec &p1,
                              QLLRvec &p0)
{
  QLLR log_apriori_prob_const_point = 0;
  int b = 0;
  for (int i = 0; i < k(j); i++) {
    log_apriori_prob_const_point +=
      ((bitmap(j)(s, i) == 0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
    b++;
  }

  b = 0;
  for (int i = 0; i < k(j); i++) {
    if (bitmap(j)(s, i) == 0) {
      p1(b) = llrcalc.jaclog(p1(b), scaled_norm
                             + log_apriori_prob_const_point);
    }
    else {
      p0(b) = llrcalc.jaclog(p0(b), scaled_norm
                             + log_apriori_prob_const_point);
    }
    b++;
  }
}

void Modulator_ND::update_LLR(const Array<QLLRvec> &logP_apriori,
                              const ivec &s, QLLR scaled_norm,
                              QLLRvec &p1, QLLRvec &p0)
{
  QLLR log_apriori_prob_const_point = 0;
  int b = 0;
  for (int j = 0; j < nt; j++) {
    for (int i = 0; i < k(j); i++) {
      log_apriori_prob_const_point +=
        ((bitmap(j)(s[j], i) == 0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
      b++;
    }
  }

  b = 0;
  for (int j = 0; j < nt; j++) {
    for (int i = 0; i < k(j); i++) {
      if (bitmap(j)(s[j], i) == 0) {
        p1(b) = llrcalc.jaclog(p1(b), scaled_norm
                               + log_apriori_prob_const_point);
      }
      else {
        p0(b) = llrcalc.jaclog(p0(b), scaled_norm
                               + log_apriori_prob_const_point);
      }
      b++;
    }
  }
}

	void Modulator_ND::marginalize_bits(itpp::QLLRvec& llr, Soft_Demod_Method method) const
{
  if(method == FULL_ENUM_LOGMAP) {
    // -- Demodulate the last 3 bits. The demodulation is hardcoded
    // -- to avoid initialization of many but tiny inner-loops
    demodllrbit0(llr[0]);
    if(nb > 1) demodllrbit1(llr[1]);
    if(nb > 2) demodllrbit2(llr[2]);
    // -- Demodulate the remaining bits except the first one
    QLLR logsum0, logsum1;
    const QLLR *const addrfirst = Qnorms._data();
    const QLLR *const addrsemilast = addrfirst + (1 << (nb - 1)), *const addrlast = addrfirst + (1 << nb);
    const QLLR *Qptr;
    for(int bi = 3; bi < nb - 1 ; bi++) { // Run the loops for bits 3,...,nb-1.
      logsum0 = -QLLR_MAX;
      logsum1 = -QLLR_MAX;
      const int forhalfdiff = 1 << bi, fordiff = 2 * forhalfdiff, fordbldiff = 2 * fordiff;
      Qptr = addrfirst;
      const QLLR *const addr1 = addrfirst + forhalfdiff, *const addr2 = addr1 + fordiff, *const addr3 = addrlast - fordiff;
      while(Qptr < addr1) logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
      while(Qptr < addr2) logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
      const QLLR *addrdyn0, *addrdyn1;
      while(Qptr < addr3) {
        addrdyn0 = Qptr + fordiff;
        addrdyn1 = Qptr + fordbldiff;
        while(Qptr < addrdyn0) logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
        while(Qptr < addrdyn1) logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
      }
      while(Qptr < addrlast) logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
      llr[bi] = logsum0 - logsum1;
    }
    // -- Demodulate the first bit
    logsum0 = -QLLR_MAX;
    logsum1 = -QLLR_MAX;
    Qptr = addrfirst;
    while(Qptr < addrsemilast) logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    while(Qptr < addrlast) logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    llr[nb - 1] = logsum0 - logsum1;
  }
  else if(method == FULL_ENUM_MAXLOG) {
    // -- Demodulate the last 3 bits. The demodulation is hardcoded to
    // -- avoid initialization of many but tiny inner-loops.
    demodmaxbit0(llr[0]);
    if(nb > 1) demodmaxbit1(llr[1]);
    if(nb > 2) demodmaxbit2(llr[2]);
    // -- Demodulate the remaining bits except the first one
    QLLR logmax0, logmax1;
    const QLLR *const addrfirst = Qnorms._data();
    const QLLR *const addrsemilast = addrfirst + (1 << (nb - 1)), *const addrlast = addrfirst + (1 << nb);
    const QLLR *Qptr;
    for(int bi = 3; bi < nb - 1; bi++) { // Run the loops for bits nb-3,nb-4,...,2.
      logmax0 = -QLLR_MAX;
      logmax1 = -QLLR_MAX;
      const int forhalfdiff = 1 << bi, fordiff = 2 * forhalfdiff, fordbldiff = 2 * fordiff;
      Qptr = addrfirst;
      const QLLR *const addr1 = addrfirst + forhalfdiff, *const addr2 = addr1 + fordiff, *const addr3 = addrlast - fordiff;
      for(; Qptr < addr1; Qptr++) logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
      for(; Qptr < addr2; Qptr++) logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
      const QLLR *addrdyn0, *addrdyn1;
      while(Qptr < addr3) {
        addrdyn0 = Qptr + fordiff;
        addrdyn1 = Qptr + fordbldiff;
        for(; Qptr < addrdyn0; Qptr++) logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
        for(; Qptr < addrdyn1; Qptr++) logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
      }
      for(; Qptr < addrlast; Qptr++) logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
      llr[bi] = logmax0 - logmax1;
    }
    // -- Demodulate the first bit
    logmax0 = -QLLR_MAX;
    logmax1 = -QLLR_MAX;
    Qptr = addrfirst;
    for(; Qptr < addrsemilast; Qptr++) logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    for(; Qptr < addrlast; Qptr++) logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    llr[nb - 1] = logmax0 - logmax1;
  }
  else it_error("Improper soft demodulation method\n.");	
}

	void Modulator_ND::demodllrbit0(itpp::QLLR& llr) const
{
  using namespace itpp;
  QLLR logsum0 = -QLLR_MAX, logsum1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 1;
  const QLLR *Qptr = addrfirst;
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  while(Qptr < addr3) {
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  }
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  llr = logsum0 - logsum1;
}

void Modulator_ND::demodllrbit1(itpp::QLLR& llr) const
{
  using namespace itpp;
  QLLR logsum0 = -QLLR_MAX, logsum1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 2;
  const QLLR *Qptr = addrfirst;

  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  while(Qptr < addr3) {
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  }
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  llr = logsum0 - logsum1;
}

void Modulator_ND::demodllrbit2(itpp::QLLR& llr) const
{
  using namespace itpp;
  QLLR logsum0 = -QLLR_MAX, logsum1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 4;
  const QLLR *Qptr = addrfirst;

  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  while(Qptr < addr3) {
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
    logsum1 = llrcalc.jaclog(*(Qptr++), logsum1);
  }
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  logsum0 = llrcalc.jaclog(*(Qptr++), logsum0);
  llr = logsum0 - logsum1;
}

void Modulator_ND::demodmaxbit0(itpp::QLLR& maxllr) const
{
  using namespace itpp;
  QLLR logmax0 = -QLLR_MAX, logmax1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 1;
  const QLLR *Qptr = addrfirst;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  while(Qptr < addr3) {
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
  }
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  maxllr = logmax0 - logmax1;
}

void Modulator_ND::demodmaxbit1(itpp::QLLR& maxllr) const
{
  using namespace itpp;
  QLLR logmax0 = -QLLR_MAX, logmax1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 2;
  const QLLR *Qptr = addrfirst;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  while(Qptr < addr3) {
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
  }
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  maxllr = logmax0 - logmax1;
}

void Modulator_ND::demodmaxbit2(itpp::QLLR& maxllr) const
{
  using namespace itpp;
  QLLR logmax0 = -QLLR_MAX, logmax1 = -QLLR_MAX;
  const QLLR *const addrfirst = Qnorms._data(), *const addr3 = addrfirst + (1 << nb) - 4;
  const QLLR *Qptr = addrfirst;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
  Qptr++;
  while(Qptr < addr3) {
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
    logmax1 = *Qptr > logmax1 ? *Qptr : logmax1;
    Qptr++;
  }
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  Qptr++;
  logmax0 = *Qptr > logmax0 ? *Qptr : logmax0;
  maxllr = logmax0 - logmax1;
}


Array<QLLRvec> Modulator_ND::probabilities(const QLLRvec &l)
{
  Array<QLLRvec> result(length(l));
  for(int i = 0; i < length(l); i++) {
    result(i) = probabilities(l(i));
  }
  return result;
}


// ----------------------------------------------------------------------
// Modulator_NRD
// ----------------------------------------------------------------------


Array<vec> Modulator_NRD::get_symbols() const
{
  Array<vec> retvec(nt);
  for(int i = 0; i < nt; ++i) {
    it_assert(M.length() == symbols.length(), "Modulator_NRD::get_symbols(): "
              "The length of M vector is different than length of the symbols vector.");
    retvec(i) = symbols(i).left(M(i));
  }
  return retvec;
}

void Modulator_NRD::modulate_bits(const bvec &bits, vec &out_symbols) const
{
  it_assert(length(bits) == sum(k), "Modulator_NRD::modulate_bits(): "
            "The number of input bits does not match.");

  out_symbols.set_size(nt);

  int b = 0;
  for(int i = 0; i < nt; ++i) {
    int symb = bin2dec(bits.mid(b, k(i)));
    out_symbols(i) = symbols(i)(bits2symbols(i)(symb));
    b += k(i);
  }
}

vec Modulator_NRD::modulate_bits(const bvec &bits) const
{
  vec result(nt);
  modulate_bits(bits, result);
  return result;
}


void Modulator_NRD::init_soft_demodulator(const itpp::mat& H_in, const double& sigma2)
{
  using namespace itpp;
  it_assert(H_in.cols() == nt, "Number of Tx antennas is wrong.\n");
  it_assert(sum(k) < 32, "Number of total bits per transmission can not be larger than 32.\n");
  it_assert(pow2i(sum(k)) == prod(M), "Modulator must use exhaustive constellations, i.e., #bits=log2(#symbs).\n");
  H = H_in;
  bitcumsum = reverse(cumsum(reverse(k)) - reverse(k)); // Shifted cummulative sum
  nb = sum(k);
  hnorms.set_size(1 << nb);
  Qnorms.set_size(1 << nb);
  hspacings.set_size(nt);
  yspacings.set_size(nt);
  bpos2cpos.set_size(nb);
  gray2dec.set_size(nt);
  gaussnorm = 2 * sigma2;
  vec startsymbvec(nt);
  for(int ci = 0; ci < nt; ci++) startsymbvec[ci] = symbols(ci)[0];
  itpp::vec Hx = H * startsymbvec;
  for(int ci = 0, bcs = 0; ci < nt; bcs += k[ci++]) {
    for(int bi = 0; bi < k[ci]; bi++) bpos2cpos[bcs + bi] = ci;
    gray2dec(ci).set_size(M[ci]);
    for(int si = 0; si < M[ci]; si++) gray2dec(ci)[si ^(si >> 1)] = si;
    yspacings(ci).set_size(M[ci] - 1);
    hspacings(ci).set_size(M[ci] - 1);
    for(int si = 0; si < M[ci] - 1; si++) {
      double xspacing = symbols(ci)[bits2symbols(ci)[(si + 1) ^((si + 1) >> 1)]];
      xspacing -= symbols(ci)[bits2symbols(ci)[si ^(si >> 1)]];
      hspacings(ci)(si) = H.get_col(ci) * xspacing;
    }
  }
  bpos2cpos = reverse(bpos2cpos);
  unsigned bitstring = 0, ind = 0;
  hxnormupdate(Hx, bitstring, ind, nb - 1);
  demod_initialized = true;
}



void Modulator_NRD::demodulate_soft_bits(const itpp::vec& y,
    const itpp::QLLRvec& llr_apr,
    itpp::QLLRvec& llr,
    Soft_Demod_Method method)
{
  using namespace itpp;

  it_assert_debug(demod_initialized, "You have to first run init_soft_demodulator().\n");
  it_assert_debug(H.rows() == y.length(), "The dimensions are not correct.\n");
  it_assert_debug(llr_apr.length() == nb, "The LLR_apr length is not correct.\n");

  // -- Prepare all the norms with the newly received vectory y
  llr.set_size(nb);
  llrapr = reverse(llr_apr); /* The bits are reversed due to the
            norm-updating functions having the rightmost bit
            as the least significant*/
  
  vec ytil = H.T() * y;
  vec startsymbvec(nt);
  for(int ci = 0; ci < nt; ci++) startsymbvec[ci] = symbols(ci)[0];
  double yx = 2*(ytil * startsymbvec);
  QLLR lapr = 0;
  for(int bi = 0; bi < nb; lapr -= llrcalc.jaclog(0, -llrapr[bi++]));

  for(int ci = 0; ci < nt; ci++)  for(int si = 0; si < M[ci] - 1; si++) {
      double xspacing = symbols(ci)[bits2symbols(ci)[(si + 1) ^((si + 1) >> 1)]];
      xspacing -= symbols(ci)[bits2symbols(ci)[si ^(si >> 1)]];
      yspacings(ci)[si] = 2*(ytil(ci) * xspacing);
    }
  unsigned bitstring = 0, ind = 0;
  yxnormupdate(yx, lapr, bitstring, ind, nb - 1); // Recursive update of all the norms
  marginalize_bits(llr,method); // Perform the appropriate bit marginalization
  llr = reverse(llr);
}

void Modulator_NRD::demodulate_soft_bits(const vec &y, const mat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori,
    Soft_Demod_Method method)
{
  switch(method) {
  case ZF_LOGMAP: {
    it_assert(H.rows() >= H.cols(), "Modulator_NRD::demodulate_soft_bits():"
              " ZF demodulation impossible for undetermined systems");
    // Set up the ZF detector
    mat Ht = H.T();
    mat inv_HtH = inv(Ht * H);
    vec shat = inv_HtH * Ht * y;
    vec h = ones(shat.size());
    for(int i = 0; i < shat.size(); ++i) {
      // noise covariance of shat
      double sigma_zf = std::sqrt(inv_HtH(i, i) * sigma2);
      shat(i) /= sigma_zf;
      h(i) /= sigma_zf;
    }
    demodulate_soft_bits(shat, h, 1.0, zeros_i(sum(k)), LLR_aposteriori);
  }
  break;
  default: {
    init_soft_demodulator(H, sigma2);
    demodulate_soft_bits(y, LLR_apriori, LLR_aposteriori, method);
  }
  }
}

QLLRvec Modulator_NRD::demodulate_soft_bits(const vec &y, const mat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    Soft_Demod_Method method)
{
  QLLRvec result;
  demodulate_soft_bits(y, H, sigma2, LLR_apriori, result, method);
  return result;
}

void Modulator_NRD::demodulate_soft_bits(const vec &y, const vec &h,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  it_assert(length(LLR_apriori) == sum(k),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");
  it_assert((length(h) == length(y)) && (length(h) == nt),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over two sigma^2"
  double moo2s2 = -1.0 / (2.0 * sigma2);

  int b = 0;
  for(int i = 0; i < nt; ++i) {
    QLLRvec bnum = -QLLR_MAX * ones_i(k(i));
    QLLRvec bdenom = bnum;
    Array<QLLRvec> logP_apriori = probabilities(LLR_apriori(b, b + k(i) - 1));
    for(int j = 0; j < M(i); ++j) {
      double norm2 = moo2s2 * sqr(y(i) - h(i) * symbols(i)(j));
      QLLR scaled_norm = llrcalc.to_qllr(norm2);
      update_LLR(logP_apriori, j, scaled_norm, i, bnum, bdenom);
    }
    LLR_aposteriori.set_subvector(b, bnum - bdenom);
    b += k(i);
  }
}

void Modulator_NRD::hxnormupdate(itpp::vec& Hx, unsigned& bitstring, unsigned& ind, unsigned bit)
{
	using namespace itpp;
	const unsigned col = bpos2cpos[bit];
	if(bit < 1) {
		hnorms[ind++] = Hx * Hx;
		unsigned oldi = gray2dec(col)[bitstring & (M[col] - 1)];
		bitstring ^= 1;
		unsigned newi = gray2dec(col)[bitstring & (M[col] - 1)];
		Hx += oldi > newi ? -hspacings(col)(newi) : hspacings(col)(oldi);
		hnorms[ind++] = Hx * Hx;
		return;
	}
	hxnormupdate(Hx, bitstring, ind, bit - 1);
	unsigned oldi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
	bitstring ^= 1 << bit;
	unsigned newi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
	Hx += oldi > newi ? -hspacings(col)(newi) : hspacings(col)(oldi);
	hxnormupdate(Hx, bitstring, ind, bit - 1);
}

void Modulator_NRD::yxnormupdate(double& yx, itpp::QLLR& lapr, unsigned& bitstring, unsigned& ind, unsigned bit)
{
  using namespace itpp;
  const unsigned col = bpos2cpos[bit];
  if(bit < 1) {
	 Qnorms[ind] = llrcalc.to_qllr((yx - hnorms[ind]) / gaussnorm) + lapr;
    ind++;
    unsigned oldi = gray2dec(col)[bitstring & (M[col] - 1)];
    bitstring ^= 1;
    unsigned newi = gray2dec(col)[bitstring & (M[col] - 1)];
    yx += oldi > newi ? -yspacings(col)[newi] : yspacings(col)[oldi];
    lapr += (bitstring & 1) ? -llrapr[bit] : llrapr[bit];
    Qnorms[ind] = llrcalc.to_qllr((yx - hnorms[ind]) / gaussnorm) + lapr;
    ind++;
    return;
  }
  yxnormupdate(yx, lapr, bitstring, ind, bit - 1);
  unsigned oldi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  bitstring ^= 1 << bit;
  unsigned newi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  yx += oldi > newi ? -yspacings(col)[newi] : yspacings(col)[oldi];
  lapr += ((bitstring >> bit) & 1) ? -llrapr[bit] : llrapr[bit];
  yxnormupdate(yx, lapr, bitstring, ind, bit - 1);
}


std::ostream &operator<<(std::ostream &os, const Modulator_NRD &mod)
{
  os << "--- REAL MIMO (NRD) CHANNEL ---------" << std::endl;
  os << "Dimension (nt):           " << mod.nt << std::endl;
  os << "Bits per dimension (k):   " << mod.k << std::endl;
  os << "Symbols per dimension (M):" << mod.M << std::endl;
  for(int i = 0; i < mod.nt; i++) {
    os << "Bitmap for dimension " << i << ": " << mod.bitmap(i) << std::endl;
    // skip printing the trailing zero
    os << "Symbol coordinates for dimension " << i << ": " << mod.symbols(i).left(mod.M(i)) << std::endl;
  }
  os << mod.get_llrcalc() << std::endl;
  return os;
}

// ----------------------------------------------------------------------
// Modulator_NCD
// ----------------------------------------------------------------------

Array<cvec> Modulator_NCD::get_symbols() const
{
  Array<cvec> retvec(nt);
  for(int i = 0; i < nt; ++i) {
    it_assert(M.length() == symbols.length(), "Modulator_NRD::get_symbols(): "
              "The length of M vector is different than length of the symbols vector.");
    retvec(i) = symbols(i).left(M(i));
  }
  return retvec;
}

void Modulator_NCD::modulate_bits(const bvec &bits, cvec &out_symbols) const
{
  it_assert(length(bits) == sum(k), "Modulator_NCD::modulate_bits(): "
            "The number of input bits does not match.");

  out_symbols.set_size(nt);

  int b = 0;
  for(int i = 0; i < nt; ++i) {
    int symb = bin2dec(bits.mid(b, k(i)));
    out_symbols(i) = symbols(i)(bits2symbols(i)(symb));
    b += k(i);
  }
}

cvec Modulator_NCD::modulate_bits(const bvec &bits) const
{
  cvec result(nt);
  modulate_bits(bits, result);
  return result;
}

void Modulator_NCD::init_soft_demodulator(const itpp::cmat& H_in, const double& sigma2)
{
  using namespace itpp;
  it_assert_debug(H_in.cols() == nt, "The number of Tx antennas is wrong.\n");
  it_assert_debug(sum(k) < 32, "Number of total bits per transmission can not be larger than 32.\n");
  it_assert_debug(pow2i(sum(k)) == prod(M), "The modulater must use exhaustive constellations, i.e., #bits=log2(#symbs).\n");
  H = H_in;
  bitcumsum = reverse(cumsum(reverse(k)) - reverse(k)); // Shifted cummulative sum
  nb = sum(k);
  hnorms.set_size(1 << nb);
  Qnorms.set_size(1 << nb);
  hspacings.set_size(nt);
  yspacings.set_size(nt);
  bpos2cpos.set_size(nb);
  gray2dec.set_size(nt);
  gaussnorm = sigma2;
  cvec startsymbvec(nt);
  for(int ci = 0; ci < nt; ci++) startsymbvec[ci] = symbols(ci)[0];
  cvec Hx = H * startsymbvec;
  for(int ci = 0, bcs = 0; ci < nt; bcs += k[ci++]) {
    for(int bi = 0; bi < k[ci]; bi++) bpos2cpos[bcs + bi] = ci;
    gray2dec(ci).set_size(M[ci]);
    for(int si = 0; si < M[ci]; si++) gray2dec(ci)[si ^(si >> 1)] = si;
    yspacings(ci).set_size(M[ci] - 1);
    hspacings(ci).set_size(M[ci] - 1);
    for(int si = 0; si < M[ci] - 1; si++) {
		 std::complex<double> xspacing = symbols(ci)[bits2symbols(ci)[(si + 1) ^((si + 1) >> 1)]];
		 xspacing -= symbols(ci)[bits2symbols(ci)[si ^(si >> 1)]];
		 hspacings(ci)(si) = H.get_col(ci) * xspacing;
    }
  }
  bpos2cpos = reverse(bpos2cpos);
  unsigned bitstring = 0, ind = 0;
  hxnormupdate(Hx, bitstring, ind, nb - 1);
  demod_initialized = true;
}

void Modulator_NCD::demodulate_soft_bits(const itpp::cvec& y,
    const itpp::QLLRvec& llr_apr,
    itpp::QLLRvec& llr,
    Soft_Demod_Method method)
{
  using namespace itpp;

  it_assert_debug(demod_initialized, "You have to first run init_soft_demodulator().\n");
  it_assert_debug(H.rows() == y.length(), "The dimensions are not correct.\n");
  it_assert_debug(llr_apr.length() == nb, "The LLR_apr length is not correct.\n");

  // -- Prepare all the norms with the newly received vectory y
  llr.set_size(nb);
  llrapr = reverse(llr_apr); /* The bits are reversed due to the
            norm-updating functions having the rightmost bit
            as the least significant*/
  cvec ytil = conj(H.H() * y);
  cvec startsymbvec(nt);
  for(int ci = 0; ci < nt; ci++) startsymbvec[ci] = symbols(ci)[0];
  double yx = 2*(ytil * startsymbvec).real();
  QLLR lapr = 0;
  for(int bi = 0; bi < nb; lapr -= llrcalc.jaclog(0, -llrapr[bi++]));
  for(int ci = 0; ci < nt; ci++)  for(int si = 0; si < M[ci] - 1; si++) {
		  std::complex<double> xspacing = symbols(ci)[bits2symbols(ci)[(si + 1) ^((si + 1) >> 1)]];
		  xspacing -= symbols(ci)[bits2symbols(ci)[si ^(si >> 1)]];
		  yspacings(ci)[si] = 2*(ytil[ci] * xspacing).real();
    }
  unsigned bitstring = 0, ind = 0;
  yxnormupdate(yx, lapr, bitstring, ind, nb - 1); // Recursive update of all the norms
  marginalize_bits(llr,method);
  llr=reverse(llr);
}

void Modulator_NCD::demodulate_soft_bits(const cvec &y, const cmat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori,
    Soft_Demod_Method method)
{
  switch(method) {
  case ZF_LOGMAP: {
	      it_assert(H.rows() >= H.cols(), "Modulator_NCD::demodulate_soft_bits():"
              " ZF demodulation impossible for undetermined systems");
    // Set up the ZF detector
    cmat Hht = H.H();
    cmat inv_HhtH = inv(Hht * H);
    cvec shat = inv_HhtH * Hht * y;
    cvec h = ones_c(shat.size());
    for(int i = 0; i < shat.size(); ++i) {
      double sigma_zf = std::sqrt(real(inv_HhtH(i, i)) * sigma2);
      shat(i) /= sigma_zf;
      h(i) /= sigma_zf;
    }
    demodulate_soft_bits(shat, h, 1.0, zeros_i(sum(k)), LLR_aposteriori);
  }
  break;
  default: {
    init_soft_demodulator(H, sigma2);
    demodulate_soft_bits(y, LLR_apriori, LLR_aposteriori, method);
  }
  }
}

QLLRvec Modulator_NCD::demodulate_soft_bits(const cvec &y, const cmat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    Soft_Demod_Method method)
{
  QLLRvec result;
  demodulate_soft_bits(y, H, sigma2, LLR_apriori, result, method);
  return result;
}


void Modulator_NCD::demodulate_soft_bits(const cvec &y, const cvec &h,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  it_assert(length(LLR_apriori) == sum(k),
            "Modulator_NCD::demodulate_soft_bits(): Wrong sizes");
  it_assert((length(h) == length(y)) && (length(h) == nt),
            "Modulator_NCD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over sigma^2"
  double moos2 = -1.0 / sigma2;

  int b = 0;
  for(int i = 0; i < nt; ++i) {
    QLLRvec bnum = -QLLR_MAX * ones_i(k(i));
    QLLRvec bdenom = -QLLR_MAX * ones_i(k(i));
    Array<QLLRvec> logP_apriori = probabilities(LLR_apriori(b, b + k(i) - 1));
    for(int j = 0; j < M(i); ++j) {
      double norm2 = moos2 * sqr(y(i) - h(i) * symbols(i)(j));
      QLLR scaled_norm = llrcalc.to_qllr(norm2);
      update_LLR(logP_apriori, j, scaled_norm, i, bnum, bdenom);
    }
    LLR_aposteriori.set_subvector(b, bnum - bdenom);
    b += k(i);
  }
}


void Modulator_NCD::hxnormupdate(itpp::cvec& Hx, unsigned& bitstring, unsigned& ind, unsigned bit)
{
  using namespace itpp;
  const unsigned col = bpos2cpos[bit];
  if(bit < 1) {
	  hnorms[ind++] = sqr(norm(Hx));
	  unsigned oldi = gray2dec(col)[bitstring & (M[col] - 1)];
    bitstring ^= 1;
    unsigned newi = gray2dec(col)[bitstring & (M[col] - 1)];
    Hx += oldi > newi ? -hspacings(col)(newi) : hspacings(col)(oldi);
    hnorms[ind++] = sqr(norm(Hx));
    return;
  }
  hxnormupdate(Hx, bitstring, ind, bit - 1);
  unsigned oldi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  bitstring ^= 1 << bit;
  unsigned newi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  Hx += oldi > newi ? -hspacings(col)(newi) : hspacings(col)(oldi);
  hxnormupdate(Hx, bitstring, ind, bit - 1);
}

void Modulator_NCD::yxnormupdate(double& yx, itpp::QLLR& lapr, unsigned& bitstring, unsigned& ind, unsigned bit)
{
  using namespace itpp;
  const unsigned col = bpos2cpos[bit];
  if(bit < 1) {
    Qnorms[ind] =  llrcalc.to_qllr((yx - hnorms[ind]) / gaussnorm) + lapr;
	 //std::cerr << dec2bin(sum(k),(int)bitstring) << " " << Qnorms[ind] << " "
		 //	 			  << llrcalc.to_qllr((2*(rec.H()*H*modulate_bits(dec2bin(sum(k),(int)bitstring)))[0].real() - hnorms[ind]) / gaussnorm) + lapr << std::endl;
    ind++;
    unsigned oldi = gray2dec(col)[bitstring & (M[col] - 1)];
    bitstring ^= 1;
    unsigned newi = gray2dec(col)[bitstring & (M[col] - 1)];
    yx += oldi > newi ? -yspacings(col)[newi] : yspacings(col)[oldi];
    lapr += (bitstring & 1) ? -llrapr[bit] : llrapr[bit];
    Qnorms[ind] = llrcalc.to_qllr((yx - hnorms[ind]) / gaussnorm) + lapr;
    ind++;
    return;
  }
  yxnormupdate(yx, lapr, bitstring, ind, bit - 1);
  unsigned oldi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  bitstring ^= 1 << bit;
  unsigned newi = gray2dec(col)[(bitstring >> bitcumsum[col]) & (M[col] - 1)];
  yx += oldi > newi ? -yspacings(col)[newi] : yspacings(col)[oldi];
  lapr += ((bitstring >> bit) & 1) ? -llrapr[bit] : llrapr[bit];
  yxnormupdate(yx, lapr, bitstring, ind, bit - 1);
}


std::ostream &operator<<(std::ostream &os, const Modulator_NCD &mod)
{
  os << "--- COMPLEX MIMO (NCD) CHANNEL --------" << std::endl;
  os << "Dimension (nt):           " << mod.nt << std::endl;
  os << "Bits per dimension (k):   " << mod.k << std::endl;
  os << "Symbols per dimension (M):" << mod.M << std::endl;
  for(int i = 0; i < mod.nt; i++) {
    os << "Bitmap for dimension " << i << ": "
       << mod.bitmap(i) << std::endl;
    os << "Symbol coordinates for dimension " << i << ": "
       << mod.symbols(i).left(mod.M(i)) << std::endl;
  }
  os << mod.get_llrcalc() << std::endl;
  return os;
}

// ----------------------------------------------------------------------
// ND_UPAM
// ----------------------------------------------------------------------

ND_UPAM::ND_UPAM(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UPAM::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UPAM::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UPAM::set_M(): Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);
  spacing.set_size(nt);

  for(int i = 0; i < nt; i++) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UPAM::set_M(): M is not a power of 2.");

    symbols(i).set_size(M(i) + 1);
    bits2symbols(i).set_size(M(i));
    bitmap(i) = graycode(k(i));
    double average_energy = (M(i) * M(i) - 1) / 3.0;
    double scaling_factor = std::sqrt(average_energy);

    for(int j = 0; j < M(i); ++j) {
      symbols(i)(j) = ((M(i) - 1) - j * 2) / scaling_factor;
      bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
    }

    // the "symbols" vector must end with a zero; only for a trick
    // exploited in update_norm()
    symbols(i)(M(i)) = 0.0;

    spacing(i) = 2.0 / scaling_factor;
  }
}

int ND_UPAM::sphere_search_SE(const vec &y_in, const mat &H,
                              const imat &zrange, double r, ivec &zhat)
{
  // The implementation of this function basically follows the
  // Schnorr-Eucner algorithm described in Agrell et al. (IEEE
  // Trans.  IT, 2002), but taking into account constellation
  // boundaries, see the "accelerated sphere decoder" in Boutros et
  // al. (IEEE Globecom, 2003).  No lattice reduction is performed.
  // Potentially the function can be speeded up by performing
  // lattice reduction, but it seems difficult to keep track of
  // constellation boundaries.

  mat R = chol(H.transpose() * H);
  mat Ri = inv(R);
  mat Q = H * Ri;
  vec y = Q.transpose() * y_in;
  mat Vi = Ri.transpose();

  int n = H.cols();
  vec dist(n);
  dist(n - 1) = 0;
  double bestdist = r * r;
  int status = -1; // search failed

  mat E = zeros(n, n);
  for(int i = 0; i < n; i++) {    // E(k,:) = y*Vi;
    for(int j = 0; j < n; j++) {
      E(i * n + n - 1) += y(j) * Vi(j + n * i);
    }
  }

  ivec z(n);
  zhat.set_size(n);
  z(n - 1) = floor_i(0.5 + E(n * n - 1));
  z(n - 1) = std::max(z(n - 1), zrange(n - 1, 0));
  z(n - 1) = std::min(z(n - 1), zrange(n - 1, 1));
  double p = (E(n * n - 1) - z(n - 1)) / Vi(n * n - 1);
  ivec step(n);
  step(n - 1) = sign_nozero_i(p);

  // Run search loop
  int k = n - 1;  // k uses natural indexing, goes from 0 to n-1

  while(true) {
    double newdist = dist(k) + p * p;

    if((newdist < bestdist) && (k != 0)) {
      for(int i = 0; i < k; i++) {
        E(k - 1 + i * n) = E(k + i * n) - p * Vi(k + i * n);
      }

      k--;
      dist(k) = newdist;
      z(k) = floor_i(0.5 + E(k + k * n));
      z(k) = std::max(z(k), zrange(k, 0));
      z(k) = std::min(z(k), zrange(k, 1));
      p = (E(k + k * n) - z(k)) / Vi(k + k * n);

      step(k) = sign_nozero_i(p);
    }
    else {
      while(true) {
        if(newdist < bestdist) {
          zhat = z;
          bestdist = newdist;
          status = 0;
        }
        else if(k == n - 1) {
          goto exit_point;
        }
        else {
          k++;
        }

        z(k) += step(k);

        if((z(k) < zrange(k, 0)) || (z(k) > zrange(k, 1))) {
          step(k) = (-step(k) - sign_nozero_i(step(k)));
          z(k) += step(k);
        }

        if((z(k) >= zrange(k, 0)) && (z(k) <= zrange(k, 1))) {
          break;
        }
      }

      p = (E(k + k * n) - z(k)) / Vi(k + k * n);
      step(k) = (-step(k) - sign_nozero_i(step(k)));
    }
  }

exit_point:
  return status;
}


int ND_UPAM::sphere_decoding(const vec &y, const mat &H, double rstart,
                             double rmax, double stepup,
                             QLLRvec &detected_bits)
{
  it_assert(H.rows() == length(y),
            "ND_UPAM::sphere_decoding(): dimension mismatch");
  it_assert(H.cols() == nt,
            "ND_UPAM::sphere_decoding(): dimension mismatch");
  it_assert(rstart > 0, "ND_UPAM::sphere_decoding(): radius error");
  it_assert(rmax > rstart, "ND_UPAM::sphere_decoding(): radius error");

  // This function can be improved, e.g., by using an ordered search.

  vec ytemp = y;
  mat Htemp(H.rows(), H.cols());
  for(int i = 0; i < H.cols(); i++) {
    Htemp.set_col(i, H.get_col(i)*spacing(i));
    ytemp += Htemp.get_col(i) * 0.5 * (M(i) - 1.0);
  }

  imat crange(nt, 2);
  for(int i = 0; i < nt; i++) {
    crange(i, 0) = 0;
    crange(i, 1) = M(i) - 1;
  }

  int status = 0;
  double r = rstart;
  ivec s(sum(M));
  while(r <= rmax) {
    status = sphere_search_SE(ytemp, Htemp, crange, r, s);

    if(status == 0) {  // search successful
      detected_bits.set_size(sum(k));
      int b = 0;
      for(int j = 0; j < nt; j++) {
        for(int i = 0; i < k(j); i++) {
          if(bitmap(j)((M(j) - 1 - s[j]), i) == 0) {
            detected_bits(b) = 1000;
          }
          else {
            detected_bits(b) = -1000;
          }
          b++;
        }
      }

      return status;
    }
    r = r * stepup;
  }

  return status;
}

// ----------------------------------------------------------------------
// ND_UQAM
// ----------------------------------------------------------------------

// The ND_UQAM (MIMO with uniform QAM) class could alternatively
// have been implemented by using a ND_UPAM class of twice the
// dimension, but this does not fit as elegantly into the class
// structure

ND_UQAM::ND_UQAM(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UQAM::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UQAM::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UQAM::set_M(): Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  L.set_size(nt);
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);

  for(int i = 0; i < nt; ++i) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UQAM::set_M(): M is not a power of 2");

    L(i) = round_i(std::sqrt(static_cast<double>(M(i))));
    it_assert(L(i)*L(i) == M(i), "ND_UQAM: constellation M must be square");

    symbols(i).set_size(M(i) + 1);
    bitmap(i).set_size(M(i), k(i));
    bits2symbols(i).set_size(M(i));
    double average_energy = (M(i) - 1) * 2.0 / 3.0;
    double scaling_factor = std::sqrt(average_energy);
    bmat gray_code = graycode(levels2bits(L(i)));

    for(int j1 = 0; j1 < L(i); ++j1) {
      for(int j2 = 0; j2 < L(i); ++j2) {
        symbols(i)(j1 * L(i) + j2) =
          std::complex<double>(((L(i) - 1) - j2 * 2.0) / scaling_factor,
                               ((L(i) - 1) - j1 * 2.0) / scaling_factor);
        bitmap(i).set_row(j1 * L(i) + j2, concat(gray_code.get_row(j1),
                          gray_code.get_row(j2)));
        bits2symbols(i)(bin2dec(bitmap(i).get_row(j1 * L(i) + j2)))
        = j1 * L(i) + j2;
      }
    }

    // must end with a zero; only for a trick exploited in
    // update_norm()
    symbols(i)(M(i)) = 0.0;
  }
}

void ND_UQAM::set_constellation_points(const int nth, const cvec& inConstellation, const ivec& in_bit2symbols)
{
  it_assert(nt > nth, "ND_UQAM::set_constellation_points(): Number of input to change is out of the size");
  it_assert(inConstellation.size() == in_bit2symbols.size(),
            "ND_UQAM::set_constellation_points(): Number of constellation and bits2symbols does not match");
  it_assert(is_even(inConstellation.size()) && (inConstellation.size() > 0),
            "ND_UQAM::set_constellation_points(): Number of symbols needs to be even and non-zero");

  symbols(nth).replace_mid(0, inConstellation);

  bits2symbols(nth) = in_bit2symbols;

  for(int m = 0; m < M(nth); ++m) {
    bitmap(nth).set_row(bits2symbols(nth)(m), dec2bin(k(nth), m));
  }

  // must end with a zero; only for a trick exploited in
  // update_norm()
  symbols(nth)(M(nth)) = 0.0;
};

// ----------------------------------------------------------------------
// ND_UPSK
// ----------------------------------------------------------------------

ND_UPSK::ND_UPSK(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UPSK::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UPSK::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UPSK::set_M() Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);

  for(int i = 0; i < nt; ++i) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UPSK::set_M(): M is not a power of 2");

    symbols(i).set_size(M(i) + 1);
    bits2symbols(i).set_size(M(i));
    bitmap(i) = graycode(k(i));

    double delta = 2.0 * pi / M(i);
    double epsilon = delta / 10000.0;

    for(int j = 0; j < M(i); ++j) {
      std::complex<double> symb
      = std::complex<double>(std::polar(1.0, delta * j));

      if(std::abs(std::real(symb)) < epsilon) {
        symbols(i)(j) = std::complex<double>(0.0, std::imag(symb));
      }
      else if(std::abs(std::imag(symb)) < epsilon) {
        symbols(i)(j) = std::complex<double>(std::real(symb), 0.0);
      }
      else {
        symbols(i)(j) = symb;
      }

      bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
    }

    // must end with a zero; only for a trick exploited in
    // update_norm()
    symbols(i)(M(i)) = 0.0;
  }
}

} // namespace itpp
