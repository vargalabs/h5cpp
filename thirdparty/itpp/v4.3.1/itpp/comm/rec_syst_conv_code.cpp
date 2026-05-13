/*!
 * \file
 * \brief Implementation of a Recursive Systematic Convolutional codec class
 * \author Pal Frenger.  QLLR support by Erik G. Larsson.
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

#include <itpp/comm/rec_syst_conv_code.h>


namespace itpp
{

//! Pointer to logarithmic branch metric function
double(*com_log)(double, double) = NULL;

//! \cond
// This wrapper is because "com_log = std::max;" below caused an error
inline double max(double x, double y) { return std::max(x, y); }
//! \endcond

// ----------------- Protected functions -----------------------------

int Rec_Syst_Conv_Code::calc_state_transition(const int instate, const int input, ivec &parity)
{
  int i, j, in = 0, temp = (gen_pol_rev(0) & (instate << 1)), parity_temp, parity_bit;

  for (i = 0; i < K; i++) {
    in = (temp & 1) ^ in;
    temp = temp >> 1;
  }
  in = in ^ input;

  parity.set_size(n - 1, false);
  for (j = 0; j < (n - 1); j++) {
    parity_temp = ((instate << 1) + in) & gen_pol_rev(j + 1);
    parity_bit = 0;
    for (i = 0; i < K; i++) {
      parity_bit = (parity_temp & 1) ^ parity_bit;
      parity_temp = parity_temp >> 1;
    }
    parity(j) = parity_bit;
  }
  return in + ((instate << 1) & ((1 << m) - 1));
}

// --------------- Public functions -------------------------
void Rec_Syst_Conv_Code::set_generator_polynomials(const ivec &gen, int constraint_length)
{
  int j;
  gen_pol = gen;
  n = gen.size();
  K = constraint_length;
  m = K - 1;
  rate = 1.0 / n;

  gen_pol_rev.set_size(n, false);
  for (int i = 0; i < n; i++) {
    gen_pol_rev(i) = reverse_int(K, gen_pol(i));
  }

  Nstates = (1 << m);
  state_trans.set_size(Nstates, 2, false);
  rev_state_trans.set_size(Nstates, 2, false);
  output_parity.set_size(Nstates, 2*(n - 1), false);
  rev_output_parity.set_size(Nstates, 2*(n - 1), false);
  int s0, s1, s_prim;
  ivec p0, p1;
  for (s_prim = 0; s_prim < Nstates; s_prim++) {
    s0 = calc_state_transition(s_prim, 0, p0);
    state_trans(s_prim, 0) = s0;
    rev_state_trans(s0, 0) = s_prim;
    for (j = 0; j < (n - 1); j++) {
      output_parity(s_prim, 2*j + 0) = p0(j);
      rev_output_parity(s0, 2*j + 0) = p0(j);
    }

    s1 = calc_state_transition(s_prim, 1, p1);
    state_trans(s_prim, 1) = s1;
    rev_state_trans(s1, 1) = s_prim;
    for (j = 0; j < (n - 1); j++) {
      output_parity(s_prim, 2*j + 1) = p1(j);
      rev_output_parity(s1, 2*j + 1) = p1(j);
    }
  }

  ln2 = std::log(2.0);

  //The default value of Lc is 1:
  Lc = 1.0;
}

void Rec_Syst_Conv_Code::set_awgn_channel_parameters(double Ec, double N0)
{
  Lc = 4.0 * std::sqrt(Ec) / N0;
}

void Rec_Syst_Conv_Code::set_scaling_factor(double in_Lc)
{
  Lc = in_Lc;
}

void Rec_Syst_Conv_Code::encode_tail(const bvec &input, bvec &tail, bmat &parity_bits)
{
  int i, j, length = input.size(), target_state;
  parity_bits.set_size(length + m, n - 1, false);
  tail.set_size(m, false);

  encoder_state = 0;
  for (i = 0; i < length; i++) {
    for (j = 0; j < (n - 1); j++) {
      parity_bits(i, j) = output_parity(encoder_state, 2 * j + int(input(i)));
    }
    encoder_state = state_trans(encoder_state, int(input(i)));
  }

  // add tail of m=K-1 zeros
  for (i = 0; i < m; i++) {
    target_state = (encoder_state << 1) & ((1 << m) - 1);
    if (state_trans(encoder_state, 0) == target_state) { tail(i) = bin(0); }
    else { tail(i) = bin(1); }
    for (j = 0; j < (n - 1); j++) {
      parity_bits(length + i, j) = output_parity(encoder_state, 2 * j + int(tail(i)));
    }
    encoder_state = target_state;
  }
  terminated = true;
}

void Rec_Syst_Conv_Code::encode(const bvec &input, bmat &parity_bits)
{
  int i, j, length = input.size();
  parity_bits.set_size(length, n - 1, false);

  encoder_state = 0;
  for (i = 0; i < length; i++) {
    for (j = 0; j < (n - 1); j++) {
      parity_bits(i, j) = output_parity(encoder_state, 2 * j + int(input(i)));
    }
    encoder_state = state_trans(encoder_state, int(input(i)));
  }
  terminated = false;
}

void Rec_Syst_Conv_Code::map_decode(const vec &rec_systematic, const mat &rec_parity, const vec &extrinsic_input,
                                    vec &extrinsic_output, bool in_terminated)
{
  double gamma_k_e, nom, den, temp0, temp1, exp_temp0, exp_temp1;
  int j, s0, s1, k, kk, s, s_prim, s_prim0, s_prim1, block_length = rec_systematic.length();
  ivec p0, p1;

  mat alpha(Nstates, block_length + 1);
  mat beta(Nstates, block_length + 1);
  mat gamma(2*Nstates, block_length + 1);
  vec denom(block_length + 1);
  denom.clear();

  extrinsic_output.set_size(block_length, false);

  if (in_terminated) { terminated = true; }

  //Calculate gamma
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      exp_temp0 = 0.0;
      exp_temp1 = 0.0;
      for (j = 0; j < (n - 1); j++) {
        exp_temp0 += 0.5 * Lc * rec_parity(kk, j) * double(1 - 2 * output_parity(s_prim, 2 * j + 0)); /* Is this OK? */
        exp_temp1 += 0.5 * Lc * rec_parity(kk, j) * double(1 - 2 * output_parity(s_prim, 2 * j + 1));
      }
      // gamma(2*s_prim+0,k) = std::exp( 0.5*(extrinsic_input(kk) + Lc*rec_systematic(kk))) * std::exp( exp_temp0 );
      // gamma(2*s_prim+1,k) = std::exp(-0.5*(extrinsic_input(kk) + Lc*rec_systematic(kk))) * std::exp( exp_temp1 );
      /* == Changed to trunc_exp() to address bug 1088420 which
         described a numerical instability problem in map_decode()
         at high SNR. This should be regarded as a temporary fix and
         it is not necessarily a waterproof one: multiplication of
         probabilities still can result in values out of
         range. (Range checking for the multiplication operator was
         not implemented as it was felt that this would sacrifice
         too much runtime efficiency.  Some margin was added to the
         numerical hardlimits below to reflect this. The hardlimit
         values below were taken as the minimum range that a
         "double" should support reduced by a few orders of
         magnitude to make sure multiplication of several values
         does not exceed the limits.)

         It is suggested to use the QLLR based log-domain() decoders
         instead of map_decode() as they are much faster and more
         numerically stable.

         EGL 8/06. == */
      gamma(2*s_prim + 0, k) = trunc_exp(0.5 * (extrinsic_input(kk) + Lc * rec_systematic(kk)) + exp_temp0);
      gamma(2*s_prim + 1, k) = trunc_exp(-0.5 * (extrinsic_input(kk) + Lc * rec_systematic(kk)) + exp_temp1);
    }
  }

  //Initiate alpha
  alpha.set_col(0, zeros(Nstates));
  alpha(0, 0) = 1.0;

  //Calculate alpha and denom going forward through the trellis
  for (k = 1; k <= block_length; k++) {
    for (s = 0; s < Nstates; s++) {
      s_prim0 = rev_state_trans(s, 0);
      s_prim1 = rev_state_trans(s, 1);
      temp0 = alpha(s_prim0, k - 1) * gamma(2 * s_prim0 + 0, k);
      temp1 = alpha(s_prim1, k - 1) * gamma(2 * s_prim1 + 1, k);
      alpha(s, k) = temp0 + temp1;
      denom(k)  += temp0 + temp1;
    }
    alpha.set_col(k, alpha.get_col(k) / denom(k));
  }

  //Initiate beta
  if (terminated) {
    beta.set_col(block_length, zeros(Nstates));
    beta(0, block_length) = 1.0;
  }
  else {
    beta.set_col(block_length, alpha.get_col(block_length));
  }

  //Calculate beta going backward in the trellis
  for (k = block_length; k >= 2; k--) {
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      beta(s_prim, k - 1) = (beta(s0, k) * gamma(2 * s_prim + 0, k)) + (beta(s1, k) * gamma(2 * s_prim + 1, k));
    }
    beta.set_col(k - 1, beta.get_col(k - 1) / denom(k));
  }

  //Calculate extrinsic output for each bit
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    nom = 0;
    den = 0;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      exp_temp0 = 0.0;
      exp_temp1 = 0.0;
      for (j = 0; j < (n - 1); j++) {
        exp_temp0 += 0.5 * Lc * rec_parity(kk, j) * double(1 - 2 * output_parity(s_prim, 2 * j + 0));
        exp_temp1 += 0.5 * Lc * rec_parity(kk, j) * double(1 - 2 * output_parity(s_prim, 2 * j + 1));
      }
      // gamma_k_e = std::exp( exp_temp0 );
      gamma_k_e = trunc_exp(exp_temp0);
      nom += alpha(s_prim, k - 1) * gamma_k_e * beta(s0, k);

      // gamma_k_e = std::exp( exp_temp1 );
      gamma_k_e = trunc_exp(exp_temp1);
      den += alpha(s_prim, k - 1) * gamma_k_e * beta(s1, k);
    }
    //      extrinsic_output(kk) = std::log(nom/den);
    extrinsic_output(kk) = trunc_log(nom / den);
  }

}

void Rec_Syst_Conv_Code::log_decode(const vec &rec_systematic, const mat &rec_parity,
                                    const vec &extrinsic_input, vec &extrinsic_output, bool in_terminated, std::string metric)
{
  if (metric == "TABLE") {
    /* Use the QLLR decoder.  This can probably be done more
    efficiently since it converts floating point vectors to QLLR.
    However we have to live with this for the time being. */
    QLLRvec rec_systematic_q  = llrcalc.to_qllr(rec_systematic);
    QLLRmat rec_parity_q = llrcalc.to_qllr(rec_parity);
    QLLRvec extrinsic_input_q = llrcalc.to_qllr(extrinsic_input);
    QLLRvec extrinsic_output_q(length(extrinsic_output));
    log_decode(rec_systematic_q, rec_parity_q, extrinsic_input_q,
               extrinsic_output_q, in_terminated);
    extrinsic_output = llrcalc.to_double(extrinsic_output_q);
    return;
  }

  double nom, den, exp_temp0, exp_temp1, rp, temp0, temp1;
  int i, j, s0, s1, k, kk, l, s, s_prim, s_prim0, s_prim1, block_length = rec_systematic.length();
  ivec p0, p1;

  //Set the internal metric:
  if (metric == "LOGMAX") { com_log = max; }
  else if (metric == "LOGMAP") { com_log = log_add; }
  else {
    it_error("Rec_Syst_Conv_Code::log_decode: Illegal metric parameter");
  }

  mat alpha(Nstates, block_length + 1);
  mat beta(Nstates, block_length + 1);
  mat gamma(2*Nstates, block_length + 1);
  extrinsic_output.set_size(block_length, false);
  vec denom(block_length + 1);
  for (k = 0; k <= block_length; k++) { denom(k) = -infinity; }

  if (in_terminated) { terminated = true; }

  //Check that Lc = 1.0
  it_assert(Lc == 1.0,
            "Rec_Syst_Conv_Code::log_decode: This function assumes that Lc = 1.0. Please use proper scaling of the input data");

  //Calculate gamma
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      exp_temp0 = 0.0;
      exp_temp1 = 0.0;
      for (j = 0; j < (n - 1); j++) {
        rp = rec_parity(kk, j);
        if (output_parity(s_prim , 2*j + 0) == 0) { exp_temp0 += rp; }
        else { exp_temp0 -= rp; }
        if (output_parity(s_prim , 2*j + 1) == 0) { exp_temp1 += rp; }
        else { exp_temp1 -= rp; }
      }
      gamma(2*s_prim + 0, k) =  0.5 * ((extrinsic_input(kk) + rec_systematic(kk)) + exp_temp0);
      gamma(2*s_prim + 1, k) = -0.5 * ((extrinsic_input(kk) + rec_systematic(kk)) - exp_temp1);
    }
  }

  //Initiate alpha
  for (j = 1; j < Nstates; j++) { alpha(j, 0) = -infinity; }
  alpha(0, 0) = 0.0;

  //Calculate alpha, going forward through the trellis
  for (k = 1; k <= block_length; k++) {
    for (s = 0; s < Nstates; s++) {
      s_prim0 = rev_state_trans(s, 0);
      s_prim1 = rev_state_trans(s, 1);
      temp0 = alpha(s_prim0, k - 1) + gamma(2 * s_prim0 + 0, k);
      temp1 = alpha(s_prim1, k - 1) + gamma(2 * s_prim1 + 1, k);
      alpha(s, k) = com_log(temp0, temp1);
      denom(k)   = com_log(alpha(s, k), denom(k));
    }
    //Normalization of alpha
    for (l = 0; l < Nstates; l++) { alpha(l, k) -= denom(k); }
  }

  //Initiate beta
  if (terminated) {
    for (i = 1; i < Nstates; i++) { beta(i, block_length) = -infinity; }
    beta(0, block_length) = 0.0;
  }
  else {
    beta.set_col(block_length, alpha.get_col(block_length));
  }

  //Calculate beta going backward in the trellis
  for (k = block_length; k >= 1; k--) {
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      beta(s_prim, k - 1) = com_log(beta(s0, k) + gamma(2 * s_prim + 0, k) , beta(s1, k) + gamma(2 * s_prim + 1, k));
    }
    //Normalization of beta
    for (l = 0; l < Nstates; l++) { beta(l, k - 1) -= denom(k); }
  }

  //Calculate extrinsic output for each bit
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    nom = -infinity;
    den = -infinity;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      exp_temp0 = 0.0;
      exp_temp1 = 0.0;
      for (j = 0; j < (n - 1); j++) {
        rp = rec_parity(kk, j);
        if (output_parity(s_prim , 2*j + 0) == 0) { exp_temp0 += rp; }
        else { exp_temp0 -= rp; }
        if (output_parity(s_prim , 2*j + 1) == 0) { exp_temp1 += rp; }
        else { exp_temp1 -= rp; }
      }
      nom = com_log(nom, alpha(s_prim, kk) + 0.5 * exp_temp0 + beta(s0, k));
      den = com_log(den, alpha(s_prim, kk) + 0.5 * exp_temp1 + beta(s1, k));
    }
    extrinsic_output(kk) = nom - den;
  }
}

void Rec_Syst_Conv_Code::log_decode_n2(const vec &rec_systematic, const vec &rec_parity,
                                       const vec &extrinsic_input, vec &extrinsic_output, bool in_terminated, std::string metric)
{
  if (metric == "TABLE") {  // use the QLLR decoder; also see comment under log_decode()
    QLLRvec rec_systematic_q  = llrcalc.to_qllr(rec_systematic);
    QLLRvec rec_parity_q = llrcalc.to_qllr(rec_parity);
    QLLRvec extrinsic_input_q = llrcalc.to_qllr(extrinsic_input);
    QLLRvec extrinsic_output_q(length(extrinsic_output));
    log_decode_n2(rec_systematic_q, rec_parity_q, extrinsic_input_q,
                  extrinsic_output_q, in_terminated);
    extrinsic_output = llrcalc.to_double(extrinsic_output_q);
    return;
  }

  //    const double INF = 10e300;  // replaced by DEFINE to be file-wide in scope
  double nom, den, exp_temp0, exp_temp1, rp;
  int k, kk, l, s, s_prim, s_prim0, s_prim1, block_length = rec_systematic.length();
  int ext_info_length = extrinsic_input.length();
  ivec p0, p1;
  double ex, norm;

  //Set the internal metric:
  if (metric == "LOGMAX") { com_log = max; }
  else if (metric == "LOGMAP") { com_log = log_add; }
  else {
    it_error("Rec_Syst_Conv_Code::log_decode_n2: Illegal metric parameter");
  }

  mat alpha(Nstates, block_length + 1);
  mat beta(Nstates, block_length + 1);
  mat gamma(2*Nstates, block_length + 1);
  extrinsic_output.set_size(ext_info_length, false);
  //denom.set_size(block_length+1,false); for (k=0; k<=block_length; k++) { denom(k) = -infinity; }

  if (in_terminated) { terminated = true; }

  //Check that Lc = 1.0
  it_assert(Lc == 1.0,
            "Rec_Syst_Conv_Code::log_decode_n2: This function assumes that Lc = 1.0. Please use proper scaling of the input data");

  //Initiate alpha
  for (s = 1; s < Nstates; s++) { alpha(s, 0) = -infinity; }
  alpha(0, 0) = 0.0;

  //Calculate alpha and gamma going forward through the trellis
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    if (kk < ext_info_length) {
      ex = 0.5 * (extrinsic_input(kk) + rec_systematic(kk));
    }
    else {
      ex = 0.5 * rec_systematic(kk);
    }
    rp = 0.5 * rec_parity(kk);
    for (s = 0; s < Nstates; s++) {
      s_prim0 = rev_state_trans(s, 0);
      s_prim1 = rev_state_trans(s, 1);
      if (output_parity(s_prim0 , 0)) { exp_temp0 = -rp; }
      else { exp_temp0 = rp; }
      if (output_parity(s_prim1 , 1)) { exp_temp1 = -rp; }
      else { exp_temp1 = rp; }
      gamma(2*s_prim0  , k) =   ex + exp_temp0;
      gamma(2*s_prim1 + 1, k) =  -ex + exp_temp1;
      alpha(s, k) = com_log(alpha(s_prim0, kk) + gamma(2 * s_prim0  , k),
                            alpha(s_prim1, kk) + gamma(2 * s_prim1 + 1, k));
      //denom(k)   = com_log( alpha(s,k), denom(k) );
    }
    norm = alpha(0, k); //norm = denom(k);
    for (l = 0; l < Nstates; l++) { alpha(l, k) -= norm; }
  }

  //Initiate beta
  if (terminated) {
    for (s = 1; s < Nstates; s++) { beta(s, block_length) = -infinity; }
    beta(0, block_length) = 0.0;
  }
  else {
    beta.set_col(block_length, alpha.get_col(block_length));
  }

  //Calculate beta going backward in the trellis
  for (k = block_length; k >= 1; k--) {
    kk = k - 1;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      beta(s_prim, kk) = com_log(beta(state_trans(s_prim, 0), k) + gamma(2 * s_prim, k),
                                 beta(state_trans(s_prim, 1), k) + gamma(2 * s_prim + 1, k));
    }
    norm = beta(0, k); //norm = denom(k);
    for (l = 0; l < Nstates; l++) { beta(l, k) -= norm; }
  }

  //Calculate extrinsic output for each bit
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    if (kk < ext_info_length) {
      nom = -infinity;
      den = -infinity;
      rp = 0.5 * rec_parity(kk);
      for (s_prim = 0; s_prim < Nstates; s_prim++) {
        if (output_parity(s_prim , 0)) { exp_temp0 = -rp; }
        else { exp_temp0 = rp; }
        if (output_parity(s_prim , 1)) { exp_temp1 = -rp; }
        else { exp_temp1 = rp; }
        nom = com_log(nom, alpha(s_prim, kk) + exp_temp0 + beta(state_trans(s_prim, 0), k));
        den = com_log(den, alpha(s_prim, kk) + exp_temp1 + beta(state_trans(s_prim, 1), k));
      }
      extrinsic_output(kk) = nom - den;
    }
  }
}

// === Below new decoder functions by EGL, using QLLR arithmetics ===========

void Rec_Syst_Conv_Code::log_decode(const QLLRvec &rec_systematic, const QLLRmat &rec_parity,
                                    const QLLRvec &extrinsic_input,
                                    QLLRvec &extrinsic_output, bool in_terminated)
{

  int nom, den, exp_temp0, exp_temp1, rp, temp0, temp1;
  int i, j, s0, s1, k, kk, l, s, s_prim, s_prim0, s_prim1, block_length = rec_systematic.length();
  //    ivec p0, p1;

  QLLRmat alpha_q(Nstates, block_length + 1);
  QLLRmat beta_q(Nstates, block_length + 1);
  QLLRmat gamma_q(2*Nstates, block_length + 1);
  extrinsic_output.set_size(block_length, false);
  QLLRvec denom_q(block_length + 1);
  for (k = 0; k <= block_length; k++) { denom_q(k) = -QLLR_MAX; }

  if (in_terminated) { terminated = true; }

  //Check that Lc = 1.0
  it_assert(Lc == 1.0,
            "Rec_Syst_Conv_Code::log_decode: This function assumes that Lc = 1.0. Please use proper scaling of the input data");

  //Calculate gamma_q
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      exp_temp0 = 0;
      exp_temp1 = 0;
      for (j = 0; j < (n - 1); j++) {
        rp = rec_parity(kk, j);
        if (output_parity(s_prim , 2*j + 0) == 0) { exp_temp0 += rp; }
        else { exp_temp0 -= rp; }
        if (output_parity(s_prim , 2*j + 1) == 0) { exp_temp1 += rp; }
        else { exp_temp1 -= rp; }
      }
      // right shift cannot be used due to implementation dependancy of how sign is handled?
      gamma_q(2*s_prim + 0, k) = ((extrinsic_input(kk) + rec_systematic(kk)) + exp_temp0) / 2;
      gamma_q(2*s_prim + 1, k) = - ((extrinsic_input(kk) + rec_systematic(kk)) - exp_temp1) / 2;
    }
  }

  //Initiate alpha_q
  for (j = 1; j < Nstates; j++) { alpha_q(j, 0) = -QLLR_MAX; }
  alpha_q(0, 0) = 0;

  //Calculate alpha_q, going forward through the trellis
  for (k = 1; k <= block_length; k++) {
    for (s = 0; s < Nstates; s++) {
      s_prim0 = rev_state_trans(s, 0);
      s_prim1 = rev_state_trans(s, 1);
      temp0 = alpha_q(s_prim0, k - 1) + gamma_q(2 * s_prim0 + 0, k);
      temp1 = alpha_q(s_prim1, k - 1) + gamma_q(2 * s_prim1 + 1, k);
      // alpha_q(s,k) = com_log( temp0, temp1 );
      // denom_q(k)   = com_log( alpha_q(s,k), denom_q(k) );
      alpha_q(s, k) = llrcalc.jaclog(temp0, temp1);
      denom_q(k)   = llrcalc.jaclog(alpha_q(s, k), denom_q(k));
    }
    //Normalization of alpha_q
    for (l = 0; l < Nstates; l++) { alpha_q(l, k) -= denom_q(k); }
  }

  //Initiate beta_q
  if (terminated) {
    for (i = 1; i < Nstates; i++) { beta_q(i, block_length) = -QLLR_MAX; }
    beta_q(0, block_length) = 0;
  }
  else {
    beta_q.set_col(block_length, alpha_q.get_col(block_length));
  }

  //Calculate beta_q going backward in the trellis
  for (k = block_length; k >= 1; k--) {
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      // beta_q(s_prim,k-1) = com_log( beta_q(s0,k) + gamma_q(2*s_prim+0,k) , beta_q(s1,k) + gamma_q(2*s_prim+1,k) );
      beta_q(s_prim, k - 1) = llrcalc.jaclog(beta_q(s0, k) + gamma_q(2 * s_prim + 0, k) , beta_q(s1, k) + gamma_q(2 * s_prim + 1, k));
    }
    //Normalization of beta_q
    for (l = 0; l < Nstates; l++) { beta_q(l, k - 1) -= denom_q(k); }
  }

  //Calculate extrinsic output for each bit
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    nom = -QLLR_MAX;
    den = -QLLR_MAX;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      s0 = state_trans(s_prim, 0);
      s1 = state_trans(s_prim, 1);
      exp_temp0 = 0;
      exp_temp1 = 0;
      for (j = 0; j < (n - 1); j++) {
        rp = rec_parity(kk, j);
        if (output_parity(s_prim , 2*j + 0) == 0) { exp_temp0 += rp; }
        else { exp_temp0 -= rp; }
        if (output_parity(s_prim , 2*j + 1) == 0) { exp_temp1 += rp; }
        else { exp_temp1 -= rp; }
      }
      // nom = com_log(nom, alpha_q(s_prim,kk) + 0.5*exp_temp0 + beta_q(s0,k) );
      // den = com_log(den, alpha_q(s_prim,kk) + 0.5*exp_temp1 + beta_q(s1,k) );
      nom = llrcalc.jaclog(nom, alpha_q(s_prim, kk) + exp_temp0 / 2 + beta_q(s0, k));
      den = llrcalc.jaclog(den, alpha_q(s_prim, kk) + exp_temp1 / 2 + beta_q(s1, k));
    }
    extrinsic_output(kk) = nom - den;
  }

}



void Rec_Syst_Conv_Code::log_decode_n2(const QLLRvec &rec_systematic,
                                       const QLLRvec &rec_parity,
                                       const QLLRvec &extrinsic_input,
                                       QLLRvec &extrinsic_output,
                                       bool in_terminated)
{
  int nom, den, exp_temp0, exp_temp1, rp;
  int k, kk, l, s, s_prim, s_prim0, s_prim1, block_length = rec_systematic.length();
  int ext_info_length = extrinsic_input.length();
  ivec p0, p1;
  int ex, norm;


  QLLRmat alpha_q(Nstates, block_length + 1);
  QLLRmat beta_q(Nstates, block_length + 1);
  QLLRmat gamma_q(2*Nstates, block_length + 1);
  extrinsic_output.set_size(ext_info_length, false);
  //denom.set_size(block_length+1,false); for (k=0; k<=block_length; k++) { denom(k) = -infinity; }

  if (in_terminated) { terminated = true; }

  //Check that Lc = 1.0
  it_assert(Lc == 1.0,
            "Rec_Syst_Conv_Code::log_decode_n2: This function assumes that Lc = 1.0. Please use proper scaling of the input data");

  //Initiate alpha
  for (s = 1; s < Nstates; s++) { alpha_q(s, 0) = -QLLR_MAX; }
  alpha_q(0, 0) = 0;

  //Calculate alpha and gamma going forward through the trellis
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    if (kk < ext_info_length) {
      ex = (extrinsic_input(kk) + rec_systematic(kk)) / 2;
    }
    else {
      ex =  rec_systematic(kk) / 2;
    }
    rp =  rec_parity(kk) / 2;
    for (s = 0; s < Nstates; s++) {
      s_prim0 = rev_state_trans(s, 0);
      s_prim1 = rev_state_trans(s, 1);
      if (output_parity(s_prim0 , 0)) { exp_temp0 = -rp; }
      else { exp_temp0 = rp; }
      if (output_parity(s_prim1 , 1)) { exp_temp1 = -rp; }
      else { exp_temp1 = rp; }
      gamma_q(2*s_prim0  , k) =   ex + exp_temp0;
      gamma_q(2*s_prim1 + 1, k) =  -ex + exp_temp1;
      alpha_q(s, k) = llrcalc.jaclog(alpha_q(s_prim0, kk) + gamma_q(2 * s_prim0  , k),
                                     alpha_q(s_prim1, kk) + gamma_q(2 * s_prim1 + 1, k));
      //denom(k)   = com_log( alpha(s,k), denom(k) );
    }
    norm = alpha_q(0, k); //norm = denom(k);
    for (l = 0; l < Nstates; l++) { alpha_q(l, k) -= norm; }
  }

  //Initiate beta
  if (terminated) {
    for (s = 1; s < Nstates; s++) { beta_q(s, block_length) = -QLLR_MAX; }
    beta_q(0, block_length) = 0;
  }
  else {
    beta_q.set_col(block_length, alpha_q.get_col(block_length));
  }

  //Calculate beta going backward in the trellis
  for (k = block_length; k >= 1; k--) {
    kk = k - 1;
    for (s_prim = 0; s_prim < Nstates; s_prim++) {
      beta_q(s_prim, kk) = llrcalc.jaclog(beta_q(state_trans(s_prim, 0), k) + gamma_q(2 * s_prim, k),
                                          beta_q(state_trans(s_prim, 1), k) + gamma_q(2 * s_prim + 1, k));
    }
    norm = beta_q(0, k); //norm = denom(k);
    for (l = 0; l < Nstates; l++) { beta_q(l, k) -= norm; }
  }

  //Calculate extrinsic output for each bit
  for (k = 1; k <= block_length; k++) {
    kk = k - 1;
    if (kk < ext_info_length) {
      nom = -QLLR_MAX;
      den = -QLLR_MAX;
      rp =  rec_parity(kk) / 2;
      for (s_prim = 0; s_prim < Nstates; s_prim++) {
        if (output_parity(s_prim , 0)) { exp_temp0 = -rp; }
        else { exp_temp0 = rp; }
        if (output_parity(s_prim , 1)) { exp_temp1 = -rp; }
        else { exp_temp1 = rp; }
        nom = llrcalc.jaclog(nom, alpha_q(s_prim, kk) + exp_temp0 + beta_q(state_trans(s_prim, 0), k));
        den = llrcalc.jaclog(den, alpha_q(s_prim, kk) + exp_temp1 + beta_q(state_trans(s_prim, 1), k));
      }
      extrinsic_output(kk) = nom - den;
    }
  }
}

void Rec_Syst_Conv_Code::set_llrcalc(LLR_calc_unit in_llrcalc)
{
  llrcalc = in_llrcalc;
}


} // namespace itpp
