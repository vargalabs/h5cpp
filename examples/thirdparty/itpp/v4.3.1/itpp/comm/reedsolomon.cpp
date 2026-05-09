/*!
 * \file
 * \brief Implementation of a Reed-Solomon codec class
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
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

#include <itpp/comm/reedsolomon.h>
#include <itpp/base/specmat.h>
#include <itpp/base/math/log_exp.h>

namespace itpp
{

//-------------------- Help Function ----------------------------

//! Local help function
GFX formal_derivate(const GFX &f)
{
  int degree = f.get_true_degree();
  int q = f.get_size();
  int i;
  GFX fprim(q, degree);
  fprim.clear();
  for (i = 0; i <= degree - 1; i += 2) {
    fprim[i] = f[i + 1];
  }
  return fprim;
}

//-------------------- Reed-Solomon ----------------------------
//A Reed-Solomon code is a q^m-ary BCH code of length n = pow(q,m)-1.
//k = pow(q,m)-1-t. This class works for q==2.
Reed_Solomon::Reed_Solomon(int in_m, int in_t, bool sys, int in_b):
  m(in_m), t(in_t), b(in_b), systematic(sys)
{
  n = pow2i(m) - 1;
  k = pow2i(m) - 1 - 2 * t;
  q = pow2i(m);
  it_assert( (b >= 0) && (b < n), "Reed_Solomon::Reed_Solomon: narrow-sense parameter restricted to 0 <= b <= n.");
  GFX x(q, (char *)"-1 0");
  ivec alphapow(1);
  g.set(q, (char *)"0");
  for (int i = 1; i <= 2 * t; i++) {
    alphapow(0) = b + i - 1;
    g *= (x - GFX(q, alphapow));
  }
}

void Reed_Solomon::encode(const bvec &uncoded_bits, bvec &coded_bits)
{
  int i, j, iterations = floor_i(static_cast<double>(uncoded_bits.length())
                                 / (k * m));
  GFX mx(q, k), cx(q, n);
  GFX r(n + 1, n - k);
  GFX uncoded_shifted(n + 1, n);
  GF mpow;
  bvec mbit(k * m), cbit(m);

  coded_bits.set_size(iterations * n * m, false);

  if (systematic)
    for (i = 0; i < n - k; i++)
      uncoded_shifted[i] = GF(n + 1, -1);

  for (i = 0; i < iterations; i++) {
    //Fix the message polynom m(x).
    for (j = 0; j < k; j++) {
      mpow.set(q, uncoded_bits.mid((i * m * k) + (j * m), m));
      mx[j] = mpow;
      if (systematic) {
        cx[j] = mx[j];
        uncoded_shifted[j + n - k] = mx[j];
      }
    }
    //Fix the outputbits cbit.
    if (systematic) {
      r = modgfx(uncoded_shifted, g);
      for (j = k; j < n; j++) {
        cx[j] = r[j - k];
      }
    }
    else {
      cx = g * mx;
    }
    for (j = 0; j < n; j++) {
      cbit = cx[j].get_vectorspace();
      coded_bits.replace_mid((i * n * m) + (j * m), cbit);
    }
  }
}

bvec Reed_Solomon::encode(const bvec &uncoded_bits)
{
  bvec coded_bits;
  encode(uncoded_bits, coded_bits);
  return coded_bits;
}

bool Reed_Solomon::decode(const bvec &coded_bits, const ivec &erasure_positions, bvec &decoded_message, bvec &cw_isvalid)
{
  bool decoderfailure, no_dec_failure;
  int j, i, kk, l, L, foundzeros, iterations = floor_i(static_cast<double>(coded_bits.length()) / (n * m));
  bvec mbit(m * k);
  decoded_message.set_size(iterations * k * m, false);
  cw_isvalid.set_length(iterations);

  GFX rx(q, n - 1), cx(q, n - 1), mx(q, k - 1), ex(q, n - 1), S(q, 2 * t), Xi(q, 2 * t), Gamma(q), Lambda(q),
      Psiprime(q), OldLambda(q), T(q), Omega(q);
  GFX dummy(q), One(q, (char*)"0"), Omegatemp(q);
  GF delta(q), tempsum(q), rtemp(q), temp(q), Xk(q), Xkinv(q);
  ivec errorpos;

  if ( erasure_positions.length() ) {
    it_assert(max(erasure_positions) < iterations*n, "Reed_Solomon::decode: erasure position is invalid.");
  }
  
  no_dec_failure = true;
  for (i = 0; i < iterations; i++) {
    decoderfailure = false;
    //Fix the received polynomial r(x)
    for (j = 0; j < n; j++) {
      rtemp.set(q, coded_bits.mid(i * n * m + j * m, m));
      rx[j] = rtemp;
    }
    // Fix the Erasure polynomial Gamma(x)
    // and replace erased coordinates with zeros
    rtemp.set(q, -1);
    ivec alphapow = - ones_i(2);
    Gamma = One;
    for (j = 0; j < erasure_positions.length(); j++) {
      rx[erasure_positions(j)] = rtemp;
      alphapow(1) = erasure_positions(j);
      Gamma *= (One - GFX(q, alphapow));
    }
    //Fix the syndrome polynomial S(x).
    S.clear();
    for (j = 1; j <= 2 * t; j++) {
      S[j] = rx(GF(q, b + j - 1));
    }
    // calculate the modified syndrome polynomial Xi(x) = Gamma * (1+S) - 1
    Xi = Gamma * (One + S) - One;
    // Apply Berlekam-Massey algorithm
    if (Xi.get_true_degree() >= 1) { //Errors in the received word
      // Iterate to find Lambda(x), which hold all error locations
      kk = 0;
      Lambda = One;
      L = 0;
      T = GFX(q, (char*)"-1 0");
      while (kk < 2 * t) {
        kk = kk + 1;
        tempsum = GF(q, -1);
        for (l = 1; l <= L; l++) {
          tempsum += Lambda[l] * Xi[kk - l];
        }
        delta = Xi[kk] - tempsum;
        if (delta != GF(q, -1)) {
          OldLambda = Lambda;
          Lambda -= delta * T;
          if (2 * L < kk) {
            L = kk - L;
            T = OldLambda / delta;
          }
        }
        T = GFX(q, (char*)"-1 0") * T;
      }
      // Find the zeros to Lambda(x)
      errorpos.set_size(Lambda.get_true_degree());
      foundzeros = 0;
      for (j = q - 2; j >= 0; j--) {
        if (Lambda(GF(q, j)) == GF(q, -1)) {
          errorpos(foundzeros) = (n - j) % n;
          foundzeros += 1;
          if (foundzeros >= Lambda.get_true_degree()) {
            break;
          }
        }
      }
      if (foundzeros != Lambda.get_true_degree()) {
        decoderfailure = true;
      }
      else { // Forney algorithm...
        //Compute Omega(x) using the key equation for RS-decoding
        Omega.set_degree(2 * t);
        Omegatemp = Lambda * (One + Xi);
        for (j = 0; j <= 2 * t; j++) {
          Omega[j] = Omegatemp[j];
        }
        //Find the error/erasure magnitude polynomial by treating them the same
        Psiprime = formal_derivate(Lambda*Gamma);
        errorpos = concat(errorpos, erasure_positions);
        ex.clear();
        for (j = 0; j < errorpos.length(); j++) {
          Xk = GF(q, errorpos(j));
          Xkinv = GF(q, 0) / Xk;
          // we calculate ex = - error polynomial, in order to avoid the 
          // subtraction when recunstructing the corrected codeword
          ex[errorpos(j)] = (Xk * Omega(Xkinv)) / Psiprime(Xkinv);
          if (b != 1) { // non-narrow-sense code needs corrected error magnitudes
            int correction_exp = ( errorpos(j)*(1-b) ) % n;
            ex[errorpos(j)] *= GF(q, correction_exp + ( (correction_exp < 0) ? n : 0 ));
          }
        }
        //Reconstruct the corrected codeword.
        // instead of subtracting the error/erasures, we calculated 
        // the negative error with 'ex' above
        cx = rx + ex;
        //Code word validation
        S.clear();
        for (j = 1; j <= 2 * t; j++) {
          S[j] = cx(GF(q, b + j - 1));
        }
        if (S.get_true_degree() >= 1) {
          decoderfailure = true;
        }
      }
    }
    else {
      cx = rx;
      decoderfailure = false;
    }
    //Find the message polynomial
    mbit.clear();
    if (decoderfailure == false) {
      if (cx.get_true_degree() >= 1) { // A nonzero codeword was transmitted
        if (systematic) {
          for (j = 0; j < k; j++) {
            mx[j] = cx[j];
          }
        }
        else {
          mx = divgfx(cx, g);
        }
        for (j = 0; j <= mx.get_true_degree(); j++) {
          mbit.replace_mid(j * m, mx[j].get_vectorspace());
        }
      }
    }
    else { //Decoder failure.
      // for a systematic code it is better to extract the undecoded message
      // from the received code word, i.e. obtaining a bit error
      // prob. p_b << 1/2, than setting all-zero (p_b = 1/2)
      if (systematic) {
        mbit = coded_bits.mid(i * n * m, k * m);
      }
      else {
        mbit = zeros_b(k);
      }
      no_dec_failure = false;
    }
    decoded_message.replace_mid(i * m * k, mbit);
    cw_isvalid(i) = (!decoderfailure);
  }
  return no_dec_failure;
}

bool Reed_Solomon::decode(const bvec &coded_bits, bvec &decoded_message, bvec &cw_isvalid)
{
  ivec   erasures(0);
  return decode(coded_bits, erasures, decoded_message, cw_isvalid);
}

void Reed_Solomon::decode(const bvec &coded_bits, bvec &decoded_bits)
{
  bvec   cw_isvalid;
  ivec   erasures(0);
  if (!decode(coded_bits, erasures, decoded_bits, cw_isvalid)) {
    for (int i = 0; i < cw_isvalid.length(); i++) {
      if (!cw_isvalid(i)) {
        decoded_bits.replace_mid(i * k * m, zeros_b(k * m));
      }
    }
  }
}

bvec Reed_Solomon::decode(const bvec &coded_bits)
{
  bvec decoded_bits;
  decode(coded_bits, decoded_bits);
  return decoded_bits;
}

// --- Soft-decision decoding is not implemented ---

void Reed_Solomon::decode(const vec &, bvec &)
{
  it_error("Reed_Solomon::decode(): Soft-decision decoding not implemented");
}

bvec Reed_Solomon::decode(const vec &)
{
  it_error("Reed_Solomon::decode(): Soft-decision decoding not implemented");
  return bvec();
}

} // namespace itpp
