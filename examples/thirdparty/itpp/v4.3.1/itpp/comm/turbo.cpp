/*!
 * \file
 * \brief Implementation of a turbo encoder/decoder class
 * \author Pal Frenger. QLLR support by Erik G. Larsson.
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

#include <itpp/comm/turbo.h>


namespace itpp
{
// -------------------------------------------------------------------------------------
// Turbo Codec
// -------------------------------------------------------------------------------------
std::string Turbo_Codec::string_from_metric(const Turbo_Codec::Metric& in_metric)
{
  if(in_metric == Metric::LOGMAX) {
    return std::string("LOGMAX");
  }
  else if(in_metric == Metric::LOGMAP) {
    return std::string("LOGMAP");
  }
  else if(in_metric == Metric::MAP) {
    return std::string("MAP");
  }
  else if(in_metric == Metric::TABLE) {
    return std::string("TABLE");
  }
  else {
    return std::string("UNKNOWN");
  }
}
void Turbo_Codec::set_parameters(ivec gen1, ivec gen2, int constraint_length, const ivec &interleaver_sequence,
                                 int in_iterations, const std::string &in_metric, double in_logmax_scale_factor,
                                 bool in_adaptive_stop,  LLR_calc_unit in_llrcalc)
{
  //Set the input parameters:
  iterations          = in_iterations;
  interleaver_size    = interleaver_sequence.size();
  Nuncoded            = interleaver_size;
  logmax_scale_factor = in_logmax_scale_factor;
  adaptive_stop       = in_adaptive_stop;

  //Check the decoding metric
  if(in_metric == "LOGMAX") {
    metric = Metric::LOGMAX;
  }
  else if(in_metric == "LOGMAP") {
    metric = Metric::LOGMAP;
  }
  else if(in_metric == "MAP") {
    metric = Metric::MAP;
  }
  else if(in_metric == "TABLE") {
    metric = Metric::TABLE;
  }
  else {
    it_error("Turbo_Codec::set_parameters: The decoder metric must be either MAP, LOGMAP, LOGMAX or TABLE");
  }

  if(logmax_scale_factor != 1.0) {
    it_assert(metric == Metric::LOGMAX, "Turbo_Codec::set_parameters: logmax_scale_factor can only be used together with LOGMAX decoding");
  }

  //The RSC Encoders:
  rscc1.set_generator_polynomials(gen1, constraint_length);
  rscc2.set_generator_polynomials(gen2, constraint_length);
  n1 = gen1.length() - 1; //Number of parity bits from rscc1
  n2 = gen2.length() - 1; //Number of parity bits from rscc2
  n_tot = 1 + n1 + n2;  //Total number of parity bits and systematic bits

  //Set the number of tail bits:
  m_tail = constraint_length - 1;

  //Calculate the number of coded bits per code-block:
  Ncoded = Nuncoded * n_tot + m_tail * (1 + n1) + m_tail * (1 + n2);

  //Set the interleaver sequence
  bit_interleaver.set_interleaver_depth(interleaver_size);
  float_interleaver.set_interleaver_depth(interleaver_size);
  bit_interleaver.set_interleaver_sequence(interleaver_sequence);
  float_interleaver.set_interleaver_sequence(interleaver_sequence);

  //Default value of the channel reliability scaling factor is 1
  Lc = 1.0;

  // LLR algebra table
  rscc1.set_llrcalc(in_llrcalc);
  rscc2.set_llrcalc(in_llrcalc);

}

void Turbo_Codec::set_interleaver(const ivec &interleaver_sequence)
{
  interleaver_size = interleaver_sequence.size();
  Nuncoded = interleaver_size;

  //Calculate the number of coded bits per code-block:
  Ncoded = Nuncoded * n_tot + m_tail * (1 + n1) + m_tail * (1 + n2);

  //Set the interleaver sequence
  bit_interleaver.set_interleaver_depth(interleaver_size);
  float_interleaver.set_interleaver_depth(interleaver_size);
  bit_interleaver.set_interleaver_sequence(interleaver_sequence);
  float_interleaver.set_interleaver_sequence(interleaver_sequence);
}

void Turbo_Codec::set_metric(std::string in_metric, double in_logmax_scale_factor, LLR_calc_unit in_llrcalc)
{
  logmax_scale_factor = in_logmax_scale_factor;

  //Check the decoding metric
  if(in_metric == "LOGMAX") {
    metric = Metric::LOGMAX;
  }
  else if(in_metric == "LOGMAP") {
    metric = Metric::LOGMAP;
  }
  else if(in_metric == "MAP") {
    metric = Metric::MAP;
  }
  else if(in_metric == "TABLE") {
    metric = Metric::TABLE;
  }
  else {
    it_error("Turbo_Codec::set_metric: The decoder metric must be either MAP, LOGMAP, LOGMAX or TABLE");
  }

  rscc1.set_llrcalc(in_llrcalc);
  rscc2.set_llrcalc(in_llrcalc);
}

void Turbo_Codec::set_iterations(int in_iterations)
{
  iterations = in_iterations;
}

void Turbo_Codec::set_adaptive_stop(bool in_adaptive_stop)
{
  adaptive_stop = in_adaptive_stop;
}

void Turbo_Codec::set_awgn_channel_parameters(double in_Ec, double in_N0)
{
  Ec = in_Ec;
  N0 = in_N0;
  Lc = 4.0 * std::sqrt(Ec) / N0;
}

void Turbo_Codec::set_scaling_factor(double in_Lc)
{
  Lc = in_Lc;
}


void Turbo_Codec::encode(const bvec &input, bvec &output)
{
  //Local variables:
  int i, k, j, no_blocks;
  int count;
  bvec input_bits, in1, in2, tail1, tail2, out;
  bmat parity1, parity2;

  //Initializations:
  no_blocks = input.length() / Nuncoded;
  output.set_size(no_blocks * Ncoded, false);

  //Set the bit counter to zero:
  count = 0;

  //Encode all code blocks:
  for(i = 0; i < no_blocks; i++) {

    //Encode one block
    input_bits = input.mid(i * Nuncoded, Nuncoded);
    encode_block(input_bits, in1, in2, parity1, parity2);

    //The data part:
    for(k = 0; k < Nuncoded; k++) {
      output(count) = in1(k);
      count++;                                //Systematic bits
      for(j = 0; j < n1; j++) { output(count) = parity1(k, j); count++; }  //Parity-1 bits
      for(j = 0; j < n2; j++) { output(count) = parity2(k, j); count++; }  //Parity-2 bits
    }

    //The first tail:
    for(k = 0; k < m_tail; k++) {
      output(count) = in1(Nuncoded + k);
      count++;                                //First systematic tail bit
      for(j = 0; j < n1; j++) { output(count) = parity1(Nuncoded + k, j); count++; }  //Parity-1 tail bits
    }

    //The second tail:
    for(k = 0; k < m_tail; k++) {
      output(count) = in2(Nuncoded + k);
      count++;                                //Second systematic tail bit
      for(j = 0; j < n2; j++) { output(count) = parity2(Nuncoded + k, j); count++; }  //Parity-2 tail bits
    }

  }

}

void Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits)
{
  ivec nrof_used_iterations;
  decode(received_signal, decoded_bits, nrof_used_iterations, true_bits);
}

void Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                         const bvec &true_bits)
{

  if((n1 == 1) && (n2 == 1) && (metric != Metric::MAP)) {
    //This is a speed optimized decoder for R=1/3 (log domain metrics only)
    decode_n3(received_signal, decoded_bits, nrof_used_iterations, true_bits);
  }
  else {

    //Local variables:
    vec rec, rec_syst1, rec_syst2;
    mat rec_parity1, rec_parity2;
    bmat decoded_bits_i;
    int no_blocks, i, j, k, nrof_used_iterations_i;
    int count;
    bool CHECK_TRUE_BITS;

    //Initilaizations:
    no_blocks = received_signal.length() / Ncoded;
    decoded_bits.set_size(no_blocks * Nuncoded, false);
    decoded_bits_i.set_size(iterations, no_blocks * Nuncoded, false);
    rec_syst1.set_size(Nuncoded + m_tail, false);
    rec_syst2.set_size(Nuncoded + m_tail, false);
    rec_syst2.clear();
    rec_parity1.set_size(Nuncoded + m_tail, n1, false);
    rec_parity2.set_size(Nuncoded + m_tail, n2, false);
    nrof_used_iterations.set_size(no_blocks, false);

    //Check the vector true_bits:
    if(true_bits.size() > 1) {
      it_assert(true_bits.size() == (Nuncoded * no_blocks), "Turbo_Codec::decode: Wrong size of input vectors");
      CHECK_TRUE_BITS = true;
    }
    else {
      CHECK_TRUE_BITS = false;
    }

    //Set the bit counter to zero:
    count = 0;

    //Itterate over all received code blocks:
    for(i = 0; i < no_blocks; i++) {

      //The data part:
      for(k = 0; k < Nuncoded; k++) {
        rec_syst1(k) = received_signal(count);
        count++;                               //Systematic bit
        for(j = 0; j < n1; j++) { rec_parity1(k, j) = received_signal(count); count++; }   //Parity-1 bits
        for(j = 0; j < n2; j++) { rec_parity2(k, j) = received_signal(count); count++; }   //Parity-2 bits
      }

      //The first tail:
      for(k = 0; k < m_tail; k++) {
        rec_syst1(Nuncoded + k) = received_signal(count);
        count++;                               //Tail 1 systematic bit
        for(j = 0; j < n1; j++) { rec_parity1(Nuncoded + k, j) = received_signal(count); count++; }   //Tail 1 parity-1 bits
      }

      //The second tail:
      for(k = 0; k < m_tail; k++) {
        rec_syst2(Nuncoded + k) = received_signal(count);
        count++;                              //Tail2 systematic bit
        for(j = 0; j < n2; j++) { rec_parity2(Nuncoded + k, j) = received_signal(count); count++; }  //Tali2 parity-2 bits
      }

      //Scale the input data if necessary:
      if(Lc != 1.0) {
        rec_syst1 *= Lc;
        rec_syst2 *= Lc;
        rec_parity1 *= Lc;
        rec_parity2 *= Lc;
      }

      //Decode the block:
      if(CHECK_TRUE_BITS) {
        decode_block(rec_syst1, rec_syst2, rec_parity1, rec_parity2, decoded_bits_i,
                     nrof_used_iterations_i, true_bits.mid(i * Nuncoded, Nuncoded));
        nrof_used_iterations(i) = nrof_used_iterations_i;
      }
      else {
        decode_block(rec_syst1, rec_syst2, rec_parity1, rec_parity2, decoded_bits_i, nrof_used_iterations_i);
        nrof_used_iterations(i) = nrof_used_iterations_i;
      }

      //Put the decoded bits in the output vector:
      decoded_bits.replace_mid(i * Nuncoded, decoded_bits_i.get_row(iterations - 1));

    }

  }

}

void Turbo_Codec::encode_block(const bvec &input, bvec &in1, bvec &in2, bmat &parity1, bmat &parity2)
{
  //Local variables:
  bvec tail1, tail2, interleaved_input;

  //Error check:
  it_assert(input.length() == Nuncoded, "Turbo_Codec::encode_block: Parameter error in Nuncoded.");

  //Initializations:
  tail1.set_size(m_tail, false);
  tail1.clear();
  tail2.set_size(m_tail, false);
  tail2.clear();
  parity1.set_size(Nuncoded + m_tail, n1, false);
  parity1.clear();
  parity2.set_size(Nuncoded + m_tail, n2, false);
  parity2.clear();
  interleaved_input.set_size(Nuncoded, false);
  interleaved_input.clear();

  //The first encoder:
  rscc1.encode_tail(input, tail1, parity1);

  //The interleaver:
  bit_interleaver.interleave(input, interleaved_input);

  //The second encoder:
  rscc2.encode_tail(interleaved_input, tail2, parity2);

  //The input vectors used to the two constituent encoders:
  in1 = concat(input, tail1);
  in2 = concat(interleaved_input, tail2);

}

void Turbo_Codec::decode_block(const vec &rec_syst1, const vec &rec_syst2, const mat &rec_parity1,
                               const mat &rec_parity2, bmat &decoded_bits_i, int &nrof_used_iterations_i,
                               const bvec &true_bits)
{
  //Local variables:
  int i;
  int count, l, k;
  vec extrinsic_input, extrinsic_output, int_rec_syst1, int_rec_syst, tmp;
  vec deint_rec_syst2, rec_syst, sub_rec_syst, Le12, Le21, Le12_int, Le21_int, L, tail1, tail2;
  bool CHECK_TRUE_BITS, CONTINUE;

  //Size initializations:
  decoded_bits_i.set_size(iterations, Nuncoded, false);
  Le12.set_size(Nuncoded + m_tail, false);
  Le21.set_size(Nuncoded + m_tail, false);
  Le21.zeros();

  //Calculate the interleaved and the deinterleaved sequences:
  float_interleaver.interleave(rec_syst1.left(interleaver_size), int_rec_syst1);
  float_interleaver.deinterleave(rec_syst2.left(interleaver_size), deint_rec_syst2);

  //Combine the results from rec_syst1 and rec_syst2 (in case some bits are transmitted several times)
  rec_syst = rec_syst1.left(interleaver_size) + deint_rec_syst2;
  int_rec_syst = rec_syst2.left(interleaver_size) + int_rec_syst1;

  //Get the two tails
  tail1 = rec_syst1.right(m_tail);
  tail2 = rec_syst2.right(m_tail);

  //Form the input vectors (including tails) to the two decoders:
  rec_syst = concat(rec_syst, tail1);
  int_rec_syst = concat(int_rec_syst, tail2);

  // Check the vector true_bits
  if(true_bits.size() > 1) {
    it_assert(true_bits.size() == Nuncoded, "Turbo_Codec::decode_block: Illegal size of input vector true_bits");
    CHECK_TRUE_BITS = true;
  }
  else {
    CHECK_TRUE_BITS = false;
  }

  if(CHECK_TRUE_BITS) {
    it_assert(adaptive_stop == false,
              "Turbo_Codec::decode_block: You can not stop iterations both adaptively and on true bits");
  }

  // Do the iterative decoding:
  nrof_used_iterations_i = iterations;
  for(i = 0; i < iterations; i++) {

    // Decode Code 1
    if(metric == Metric::MAP) {
      rscc1.map_decode(rec_syst, rec_parity1, Le21, Le12, true);
    }
    else if((metric == Metric::LOGMAX) || (metric == Metric::LOGMAP) || (metric == Metric::TABLE)) {
      rscc1.log_decode(rec_syst, rec_parity1, Le21, Le12, true, string_from_metric(metric));
      if(logmax_scale_factor != 1.0) {
        Le12 *= logmax_scale_factor;
      }
    }
    else {
      it_error("Turbo_Codec::decode_block: Illegal metric value");
    }

    // Interleave the extrinsic information:
    float_interleaver.interleave(Le12.left(interleaver_size), tmp);
    Le12_int = concat(tmp, zeros(Le12.size() - interleaver_size));

    // Decode Code 2
    if(metric == Metric::MAP) {
      rscc2.map_decode(int_rec_syst, rec_parity2, Le12_int, Le21_int, true);
    }
    else if((metric == Metric::LOGMAX) || (metric == Metric::LOGMAP) || (metric == Metric::TABLE)) {
      rscc2.log_decode(int_rec_syst, rec_parity2, Le12_int, Le21_int, true, string_from_metric(metric));
      if(logmax_scale_factor != 1.0) {
        Le21_int *= logmax_scale_factor;
      }
    }
    else {
      it_error("Turbo_Codec::decode_block: Illegal metric value");
    }

    // De-interleave the extrinsic information:
    float_interleaver.deinterleave(Le21_int.left(interleaver_size), tmp);
    Le21 = concat(tmp, zeros(Le21_int.size() - interleaver_size));

    // Take bit decisions
    L = rec_syst + Le21 + Le12;
    count = 0;
    for(l = 0; l < Nuncoded; l++) {
      (L(l) > 0.0) ? (decoded_bits_i(i, count) = bin(0)) : (decoded_bits_i(i, count) = bin(1));
      count++;
    }

    //Check if it is possible to stop iterating early:
    CONTINUE = true;
    if(i < (iterations - 1)) {

      if(CHECK_TRUE_BITS) {
        CONTINUE = false;
        for(k = 0; k < Nuncoded; k++) { if(true_bits(k) != decoded_bits_i(i, k)) { CONTINUE = true; break; } }
      }

      if((adaptive_stop) && (i > 0)) {
        CONTINUE = false;
        for(k = 0; k < Nuncoded; k++) { if(decoded_bits_i(i - 1, k) != decoded_bits_i(i, k)) { CONTINUE = true; break; } }
      }

    }

    //Check if iterations shall continue:
    if(CONTINUE == false) {
      //Copy the results from current iteration to all following iterations:
      for(k = (i + 1); k < iterations; k++) {
        decoded_bits_i.set_row(k, decoded_bits_i.get_row(i));
        nrof_used_iterations_i = i + 1;
      }
      break;
    }

  }

}

void Turbo_Codec::decode_n3(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                            const bvec &true_bits)
{
  //Local variables:
  vec rec, rec_syst1, int_rec_syst1, rec_syst2;
  vec rec_parity1, rec_parity2;
  vec extrinsic_input, extrinsic_output, Le12, Le21, Le12_int, Le21_int, L;
  bvec temp_decoded_bits;
  int no_blocks, i, j, k, l, nrof_used_iterations_i;
  int count, count_out;
  bool CHECK_TRUE_BITS, CONTINUE;

  //Initializations:
  no_blocks = received_signal.length() / Ncoded;
  decoded_bits.set_size(no_blocks * Nuncoded, false);
  rec_syst1.set_size(Nuncoded + m_tail, false);
  rec_syst2.set_size(Nuncoded + m_tail, false);
  rec_syst2.clear();
  rec_parity1.set_size(Nuncoded + m_tail, false);
  rec_parity2.set_size(Nuncoded + m_tail, false);
  temp_decoded_bits.set_size(Nuncoded, false);
  decoded_bits_previous_iteration.set_size(Nuncoded, false);
  nrof_used_iterations.set_size(no_blocks, false);

  //Size initializations:
  Le12.set_size(Nuncoded, false);
  Le21.set_size(Nuncoded, false);

  //Set the bit counter to zero:
  count = 0;
  count_out = 0;

  // Check the vector true_bits
  if(true_bits.size() > 1) {
    it_assert(true_bits.size() == Nuncoded * no_blocks, "Turbo_Codec::decode_n3: Illegal size of input vector true_bits");
    CHECK_TRUE_BITS = true;
  }
  else {
    CHECK_TRUE_BITS = false;
  }

  if(CHECK_TRUE_BITS) {
    it_assert(adaptive_stop == false,
              "Turbo_Codec::decode_block: You can not stop iterations both adaptively and on true bits");
  }

  //Iterate over all received code blocks:
  for(i = 0; i < no_blocks; i++) {

    //Reset extrinsic data:
    Le21.zeros();

    //The data part:
    for(k = 0; k < Nuncoded; k++) {
      rec_syst1(k)   = received_signal(count);
      count++; //Systematic bit
      rec_parity1(k) = received_signal(count);
      count++; //Parity-1 bits
      rec_parity2(k) = received_signal(count);
      count++; //Parity-2 bits
    }

    //The first tail:
    for(k = 0; k < m_tail; k++) {
      rec_syst1(Nuncoded + k)   = received_signal(count);
      count++; //Tail 1 systematic bit
      rec_parity1(Nuncoded + k) = received_signal(count);
      count++; //Tail 1 parity-1 bits
    }

    //The second tail:
    for(k = 0; k < m_tail; k++) {
      rec_syst2(Nuncoded + k)   = received_signal(count);
      count++; //Tail2 systematic bit
      rec_parity2(Nuncoded + k) = received_signal(count);
      count++; //Tali2 parity-2 bits
    }

    float_interleaver.interleave(rec_syst1.left(Nuncoded), int_rec_syst1);
    rec_syst2.replace_mid(0, int_rec_syst1);

    //Scale the input data if necessary:
    if(Lc != 1.0) {
      rec_syst1   *= Lc;
      rec_syst2   *= Lc;
      rec_parity1 *= Lc;
      rec_parity2 *= Lc;
    }

    //Decode the block:
    CONTINUE = true;
    nrof_used_iterations_i = iterations;
    for(j = 0; j < iterations; j++) {

      rscc1.log_decode_n2(rec_syst1, rec_parity1, Le21, Le12, true, string_from_metric(metric));
      if(logmax_scale_factor != 1.0) { Le12 *= logmax_scale_factor; }
      float_interleaver.interleave(Le12, Le12_int);

      rscc2.log_decode_n2(rec_syst2, rec_parity2, Le12_int, Le21_int, true, string_from_metric(metric));
      if(logmax_scale_factor != 1.0) { Le21_int *= logmax_scale_factor; }
      float_interleaver.deinterleave(Le21_int, Le21);

      if(adaptive_stop) {
        L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
        for(l = 0; l < Nuncoded; l++) {(L(l) > 0.0) ? (temp_decoded_bits(l) = bin(0)) : (temp_decoded_bits(l) = bin(1)); }
        if(j == 0) { decoded_bits_previous_iteration = temp_decoded_bits; }
        else {
          if(temp_decoded_bits == decoded_bits_previous_iteration) {
            CONTINUE = false;
          }
          else if(j < (iterations - 1)) {
            decoded_bits_previous_iteration = temp_decoded_bits;
          }
        }
      }

      if(CHECK_TRUE_BITS) {
        L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
        for(l = 0; l < Nuncoded; l++) {(L(l) > 0.0) ? (temp_decoded_bits(l) = bin(0)) : (temp_decoded_bits(l) = bin(1)); }
        if(temp_decoded_bits == true_bits.mid(i * Nuncoded, Nuncoded)) {
          CONTINUE = false;
        }
      }

      if(CONTINUE == false) { nrof_used_iterations_i = j + 1; break; }

    }

    //Take final bit decisions
    L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
    for(l = 0; l < Nuncoded; l++) {
      (L(l) > 0.0) ? (decoded_bits(count_out) = bin(0)) : (decoded_bits(count_out) = bin(1));
      count_out++;
    }

    nrof_used_iterations(i) = nrof_used_iterations_i;

  }

}

// -------------------------------------------------------------------------------------
// Punctured Turbo Codec
// -------------------------------------------------------------------------------------

void Punctured_Turbo_Codec::set_parameters(ivec gen1, ivec gen2, int constraint_length, const ivec &interleaver_sequence, bmat &pmatrix, int in_iterations, std::string in_metric, double in_logmax_scale_factor, bool in_adaptive_stop, itpp::LLR_calc_unit lcalc)
{
  Turbo_Codec::set_parameters(gen1, gen2, constraint_length, interleaver_sequence, in_iterations, in_metric, in_logmax_scale_factor, in_adaptive_stop, lcalc);
  set_puncture_matrix(pmatrix);
}

void Punctured_Turbo_Codec::set_puncture_matrix(const bmat &pmatrix)
{
  int p, j;

  punct_total = 0;
  punct_total2 = 0;

  it_error_if(pmatrix.rows() != n_tot || pmatrix.cols() == 0, "Wrong size of puncture matrix");
  puncture_matrix = pmatrix;
  Period = puncture_matrix.cols();

  // all rows
  for(j = 0; j < n_tot; j++) {
    for(p = 0; p < Period; p++)
      punct_total += static_cast<int>(puncture_matrix(j, p));
  }
  // systematic bits
  for(p = 0; p < Period; p++)
    punct_total2 += static_cast<int>(puncture_matrix(0, p));
  punct_total1 = punct_total2;
  // 1st code parity bits
  for(j = 1; j < n1 + 1; j++) {
    for(p = 0; p < Period; p++)
      punct_total1 += static_cast<int>(puncture_matrix(j, p));
  }
  // 2nd code parity bits
  for(j = 1 + n1; j < n_tot; j++)  {
    for(p = 0; p < Period; p++)
      punct_total2 += static_cast<int>(puncture_matrix(j, p));
  }

  // nominal rate
  rate = Period / static_cast<double>(punct_total);
  calculate_punctured_size();
}


double Punctured_Turbo_Codec::get_rate(bool nominal)
{
  if(nominal) return rate;
  else {
    if(Period == 0)
      return static_cast<double>(Nuncoded) / Ncoded;
    else
      return static_cast<double>(Nuncoded) / pNcoded;
  }
}


bvec Punctured_Turbo_Codec::encode(const bvec &input)
{
  bvec coded_bits;

  encode(input, coded_bits);
  return coded_bits;
}

bvec Punctured_Turbo_Codec::decode(const vec &received_signal)
{
  bvec decoded_bits;

  decode(received_signal, decoded_bits);
  return decoded_bits;
}

void Punctured_Turbo_Codec::encode(const bvec &input, bvec &output)
{
  it_assert(Period != 0, "Punctured_Turbo_Codec: puncture matrix is not set");

  Turbo_Codec::encode(input, output);

  int i, k, p, j, p1;
  int no_blocks = output.size() / Ncoded;
  int count = 0, count_p = 0;

  for(k = 0; k < no_blocks; k++) {
    p = 0;
    // data
    for(i = 0; i < Nuncoded; i++) {
      for(j = 0; j < n_tot; j++) {
        if(puncture_matrix(j, p) == bin(1)) {
          output(count_p) = output(count);
          count_p++;
        }
        count++;
      }
      p = (p + 1) % Period;
    }
    p1 = p;

    //The first tail:
    for(i = 0; i < m_tail; i++) {
      for(j = 0; j < n1 + 1; j++) {
        if(puncture_matrix(j, p) == bin(1)) {
          output(count_p) = output(count);
          count_p++;
        }
        count++;
      }
      p = (p + 1) % Period;
    }
    //The second tail:
    for(i = 0; i < m_tail; i++) {
      // systematic bit
      if(puncture_matrix(0, p1) == bin(1)) {
        output(count_p) = output(count);
        count_p++;
      }
      count++;
      // parity
      for(j = n1 + 1; j < n_tot; j++) {
        if(puncture_matrix(j, p1) == bin(1)) {
          output(count_p) = output(count);
          count_p++;
        }
        count++;
      }
      p1 = (p1 + 1) % Period;
    }
  }  //for

  output.set_size(count_p, true);
}


void Punctured_Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations, const bvec &true_bits)
{
  int i, k, p, j, p1;
  int index = 0, index_p = 0;
  int no_blocks = received_signal.size() / pNcoded;
  vec temp(no_blocks * Ncoded);

  it_assert(Period != 0, "Punctured_Turbo_Codec: puncture matrix is not set");
  it_assert(no_blocks * pNcoded == received_signal.size(), "Punctured_Turbo_Codec: received vector is not an integer multiple of encoded block");
  for(i = 0; i < no_blocks; i++)  {
    p = 0;
    // data
    for(k = 0; k < Nuncoded; k++)  {
      for(j = 0; j < n_tot; j++) {
        if(puncture_matrix(j, p) == bin(1)) {
          temp(index) = received_signal(index_p);
          index_p++;
        }
        else { // insert dummy symbols with same contribution for 0 and 1
          temp(index) = 0;
        }
        index++;
      }
      p = (p + 1) % Period;
    } // for
    p1 = p;

    // 1st code tail
    for(k = 0; k < m_tail; k++) {
      for(j = 0; j < n1 + 1; j++) {
        if(puncture_matrix(j, p) == bin(1)) {
          temp(index) = received_signal(index_p);
          index_p++;
        }
        else { // insert dummy symbols with same contribution for 0 and 1
          temp(index) = 0;
        }
        index++;
      }
      p = (p + 1) % Period;
    } // for
    // 2nd code tail
    for(k = 0; k < m_tail; k++) {
      // systematic bits
      if(puncture_matrix(0, p1) == bin(1)) {
        temp(index) = received_signal(index_p);
        index_p++;
      }
      else { // insert dummy symbols with same contribution for 0 and 1
        temp(index) = 0;
      }
      index++;
      // parity bits
      for(j = n1 + 1; j < n_tot; j++) {
        if(puncture_matrix(j, p1) == bin(1)) {
          temp(index) = received_signal(index_p);
          index_p++;
        }
        else { // insert dummy symbols with same contribution for 0 and 1
          temp(index) = 0;
        }
        index++;
      }
      p1 = (p1 + 1) % Period;
    } //2nd tail
  }  // for

  Turbo_Codec::decode(temp, decoded_bits, nrof_used_iterations, true_bits);
}

void Punctured_Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits)
{
  ivec nrof_used_iterations;
  decode(received_signal, decoded_bits, nrof_used_iterations, true_bits);
}

void Punctured_Turbo_Codec::calculate_punctured_size(void)
{
  int i, j, ii, p = 0, p1;

  if(Period == 0)
    pNcoded = Ncoded;
  else  {
    i = (Nuncoded / Period);
    ii = i * punct_total;
    i *= Period;
    for(; i < Nuncoded; i++) {
      for(j = 0; j < n_tot; j++)
        if(puncture_matrix(j, p) == bin(1))  ii++;
      p = (p + 1) % Period;
    }
    p1 = p;

    // first tail
    for(i = 0; i < m_tail; i++) {
      for(j = 0; j < n1 + 1; j++)
        if(puncture_matrix(j, p) == bin(1))  ii++;
      p = (p + 1) % Period;
    }
    // second tail
    for(i = 0; i < m_tail; i++) {
      for(j = 0; j < n_tot; j++)  {
        if(puncture_matrix(j, p1) == bin(1))  ii++;
        if(j == 0) j += n1;
      }
      p1 = (p1 + 1) % Period;
    }

    pNcoded = ii;
  }
}

int calculate_uncoded_size(Punctured_Turbo_Codec &tc, int punctured_size, int &fill_bits)
{
  // fill_bits - number of bits that must be added at the end of encoded block (in order to obtain punctured_size length vector)

  int Nuncoded;

  if(tc.Period == 0) {
    Nuncoded = (punctured_size - tc.m_tail * (tc.n_tot + 1)) / tc.n_tot;
    fill_bits = punctured_size - (Nuncoded * tc.n_tot + tc.m_tail * (tc.n_tot + 1));
  }
  else  {
    int i, j, ii, p, p1, no_pblocks;
    // uncoded - // no_pblocks might be too small
    j = static_cast<int>(std::ceil(static_cast<double>(tc.m_tail * (tc.punct_total1 + tc.punct_total2)) / tc.Period));
    no_pblocks = (punctured_size - j) / tc.punct_total;
    ii = punctured_size - no_pblocks * tc.punct_total - j;

    for(i = 0; i < 2 * tc.Period; i++) {
      for(j = 0; j < tc.n_tot; j++)
        if(tc.puncture_matrix(j, i % tc.Period) == bin(1))  ii--;
      if(ii < 0) break;
    }
    Nuncoded = no_pblocks * tc.Period + i;

    // punctured (from uncoded)
    no_pblocks = (Nuncoded / tc.Period);
    ii = no_pblocks * tc.punct_total;
    p = 0;
    for(i = no_pblocks * tc.Period; i < Nuncoded; i++) {
      for(j = 0; j < tc.n_tot; j++)
        if(tc.puncture_matrix(j, p) == bin(1))  ii++;
      p = (p + 1) % tc.Period;
    }
    p1 = p;
    // first tail
    for(i = 0; i < tc.m_tail; i++) {
      for(j = 0; j < tc.n1 + 1; j++)
        if(tc.puncture_matrix(j, p1) == bin(1))  ii++;
      p1 = (p1 + 1) % tc.Period;
    }
    // second tail
    for(i = 0; i < tc.m_tail; i++) {
      for(j = 0; j < tc.n_tot; j++)  {
        if(tc.puncture_matrix(j, p1) == bin(1))  ii++;
        if(j == 0) j += tc.n1;
      }
      p1 = (p1 + 1) % tc.Period;
    }
    fill_bits =  punctured_size - ii;
  }

  return Nuncoded;
}


// -------------------------------------------------------------------------------------
// Special interleaver sequence generators
// -------------------------------------------------------------------------------------

ivec wcdma_turbo_interleaver_sequence(int interleaver_size)
{
  const int MAX_INTERLEAVER_SIZE = 5114;
  const int MIN_INTERLEAVER_SIZE = 40;
  int K;  //Interleaver size
  int R;  //Number of rows of rectangular matrix
  int C;  //Number of columns of rectangular matrix
  int p;  //Prime number
  int v;  //Primitive root
  ivec s; //Base sequence for intra-row permutation
  ivec q; //Minimum prime integers
  ivec r; //Permuted prime integers
  ivec T; //Inter-row permutation pattern
  imat U; //Intra-row permutation patter
  ivec I; //The interleaver sequence
  ivec primes, roots, Pat1, Pat2, Pat3, Pat4, Isort;
  int i, j, qj, temp, row, col, index, count;

  if(interleaver_size > MAX_INTERLEAVER_SIZE) {

    I = sort_index(randu(interleaver_size));
    return I;

  }
  else {

    p = 0;
    v = 0;

    //Check the range of the interleaver size:
    it_assert(interleaver_size <= MAX_INTERLEAVER_SIZE, "wcdma_turbo_interleaver_sequence: The interleaver size is to large");
    it_assert(interleaver_size >= MIN_INTERLEAVER_SIZE, "wcdma_turbo_interleaver_sequence: The interleaver size is to small");

    K = interleaver_size;

    //Definitions of primes and associated primitive roots:
    primes = "2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257";
    roots = "0 0 0 3 2 2 3 2 5 2 3 2 6 3 5 2 2 2 2 7 5 3 2 3 5 2 5 2 6 3 3 2 3 2 2 6 5 2 5 2 2 2 19 5 2 3 2 3 2 6 3 7 7 6 3";

    //Determine R
    if((K >= 40) && (K <= 159)) {
      R = 5;
    }
    else if(((K >= 160) && (K <= 200)) || ((K >= 481) && (K <= 530))) {
      R = 10;
    }
    else {
      R = 20;
    }

    //Determine C
    if((K >= 481) && (K <= 530)) {
      p = 53;
      v = 2;
      C = p;
    }
    else {
      //Find minimum prime p such that (p+1) - K/R >= 0 ...
      for(i = 0; i < primes.length(); i++) {
        if((double(primes(i) + 1) - double(K) / double(R)) >= 0.0) {
          p = primes(i);
          v = roots(i);
          break;
        }
      }
      //... and etermine C such that
      if((double(p) - double(K) / double(R)) >= 0.0) {
        if((double(p) - 1.0 - double(K) / double(R)) >= 0.0) {
          C = p - 1;
        }
        else {
          C = p;
        }
      }
      else {
        C = p + 1;
      }
    }

    //Construct the base sequencs s for intra-row permutaions
    s.set_size(p - 1, false);
    s.clear();
    s(0) = 1;
    for(i = 1; i <= (p - 2); i++) {
      s(i) = mod(v * s(i - 1), p);
    }

    //Let q(0) = 1 be the first prime integer in {q(j)}, and select the consecutive
    //minimum prime integers {q(j)}, j = 1, 2, ..., (R-1) such that gcd( q(j), p-1) == 1, q(j) > 6, and q(j) > q(j-1)
    q.set_size(R, false);
    q.clear();
    q(0) = 1;
    for(j = 1; j <= (R - 1); j++) {
      for(i = 0; i < primes.length(); i++) {
        qj = primes(i);
        if((qj > 6) && (qj > q(j - 1))) {
          if(gcd(qj, p - 1) == 1) {
            q(j) = qj;
            break;
          }
        }
      }
    }

    //Definitions of Pat1, Pat2, Pat3, and Pat4:
    Pat1 = "19 9 14 4 0 2 5 7 12 18 10 8 13 17 3 1 16 6 15 11";
    Pat2 = "19 9 14 4 0 2 5 7 12 18 16 13 17 15 3 1 6 11 8 10";
    Pat3 = "9 8 7 6 5 4 3 2 1 0";
    Pat4 = "4 3 2 1 0";

    //T(j) is the inter-row permutation patters defined as one of the following four
    //kinds of patterns: Pat1, Pat2, Pat3, and Pat4 depending on the number of input bits K
    if(K >= 3211) {
      T = Pat1;
    }
    else if(K >= 3161) {
      T = Pat2;
    }
    else if(K >= 2481) {
      T = Pat1;
    }
    else if(K >= 2281) {
      T = Pat2;
    }
    else if(K >= 531) {
      T = Pat1;
    }
    else if(K >= 481) {
      T = Pat3;
    }
    else if(K >= 201) {
      T = Pat1;
    }
    else if(K >= 160) {
      T = Pat3;
    }
    else {
      T = Pat4;
    }

    //Permute {q(j)} to make {r(j)} such that r(T(j)) = q(j), j = 0, 1, ..., (R-1),
    //where T(j) indicates the original row position of the j-th permuted row
    r.set_size(R, false);
    r.clear();
    for(j = 0; j <= (R - 1); j++) {
      r(T(j)) = q(j);
    }

    //U(j,i) is the input bit position of i-th output after the permutation of j-th row
    //Perform the j-th (j=0, 1, 2, ..., (R-1)) intra-row permutation as
    U.set_size(R, C, false);
    U.clear();
    if(C == p) {
      for(j = 0; j <= (R - 1); j++) {
        for(i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1));
        }
        U(j, p - 1) = 0;
      }
    }
    else if(C == (p + 1)) {
      for(j = 0; j <= (R - 1); j++) {
        for(i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1));
        }
        U(j, p - 1) = 0;
        U(j, p) = p;
      }
      if(K == (C * R)) {
        temp = U(R - 1, p);
        U(R - 1, p) = U(R - 1, 0);
        U(R - 1, 0) = temp;
      }
    }
    else if(C == (p - 1)) {
      for(j = 0; j <= (R - 1); j++) {
        for(i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1)) - 1;
        }
      }
    }

    //Calculate the interleaver sequence:
    I.set_size(K, false);
    I.clear();
    count = 0;
    for(i = 0; i < C; i++) {
      for(j = 0; j < R; j++) {
        row = T(j);
        col = U(row, i);
        index = row * C + col;
        if(index < K) {
          I(count) = index;
          count++;
        }
      }
    }

    return I;
  }
}

ivec lte_turbo_interleaver_sequence(int interleaver_size)
// for standard see pp. 14 http://www.3gpp.org/FTP/Specs/latest/Rel-10/36_series/36212-a50.zip
{

  //Definitions of block lengths and associated f1 and f2 factors:
  ivec block_lengths("40 48 56 64 72 80 88 96 104 112 120 128 136 144 152 160 168 176 184 192 200 208 216 224 232 240 248 256 264 272 280 288 296 304 312 320 328 336 344 352 360 368 376 384 392 400 408 416 424 432 440 448 456 464 472 480 488 496 504 512 528 544 560 576 592 608 624 640 656 672 688 704 720 736 752 768 784 800 816 832 848 864 880 896 912 928 944 960 976 992 1008 1024 1056 1088 1120 1152 1184 1216 1248 1280 1312 1344 1376 1408 1440 1472 1504 1536 1568 1600 1632 1664 1696 1728 1760 1792 1824 1856 1888 1920 1952 1984 2016 2048 2112 2176 2240 2304 2368 2432 2496 2560 2624 2688 2752 2816 2880 2944 3008 3072 3136 3200 3264 3328 3392 3456 3520 3584 3648 3712 3776 3840 3904 3968 4032 4096 4160 4224 4288 4352 4416 4480 4544 4608 4672 4736 4800 4864 4928 4992 5056 5120 5184 5248 5312 5376 5440 5504 5568 5632 5696 5760 5824 5888 5952 6016 6080 6144");
  ivec f1_factors(" 3  7 19  7  7 11  5 11   7  41 103  15   9  17   9  21 101  21  57  23  13  27  11  27  85  29  33  15  17  33 103  19  19  37  19  21  21 115 193  21 133  81  45  23 243 151 155  25  51  47  91  29  29 247  29  89  91 157  55  31  17  35 227  65  19  37  41  39 185  43  21 155  79 139  23 217  25  17 127  25 239  17 137 215  29  15 147  29  59  65   55   31   17  171   67   35   19   39   19  199   21  211   21   43  149   45   49   71   13   17   25  183   55  127   27   29   29   57   45   31   59  185  113   31   17  171  209  253  367  265  181   39   27  127  143   43   29   45  157   47   13  111  443   51   51  451  257   57  313  271  179  331  363  375  127   31   33   43   33  477   35  233  357  337   37   71   71   37   39  127   39   39   31  113   41  251   43   21   43   45   45  161   89  323   47   23   47  263");
  ivec f2_factors("10 12 42 16 18 20 22 24  26  84  90  32  34 108  38 120  84  44  46  48  50  52  36  56  58  60  62  32 198  68 210  36  74  76  78 120  82  84  86  44  90  46  94  48  98  40 102  52 106  72 110 168 114  58 118 180 122  62  84  64  66  68 420  96  74  76 234  80  82 252  86  44 120  92  94  48  98  80 102  52 106  48 110 112 114  58 118  60 122 124   84   64   66  204  140   72   74   76   78  240   82  252   86   88   60   92  846   48   28   80  102  104  954   96  110  112  114  116  354  120  610  124  420   64   66  136  420  216  444  456  468   80  164  504  172   88  300   92  188   96   28  240  204  104  212  192  220  336  228  232  236  120  244  248  168   64  130  264  134  408  138  280  142  480  146  444  120  152  462  234  158   80   96  902  166  336  170   86  174  176  178  120  182  184  186   94  190  480");
  const int MAX_INTERLEAVER_SIZE = 6144;
  const int MIN_INTERLEAVER_SIZE = 40;

  // Check the range of the interleaver size:
  it_assert(interleaver_size <= MAX_INTERLEAVER_SIZE, "lte_turbo_interleaver_sequence: The interleaver size is too large");
  it_assert(interleaver_size >= MIN_INTERLEAVER_SIZE, "lte_turbo_interleaver_sequence: The interleaver size is too small");

  // Check whether the given interleaver size is correct:
  int left, right, index, temp;
  bool search = true;

  // do a binary search for interleaver_size in block_lengths
  left = 0;
  right = block_lengths.size() - 1;
  temp = 0;
  while((search) && (left <= right)) {
    index = (left + right) / 2;
    temp = block_lengths(index);
    if(temp == interleaver_size) {
      search = false;
    }
    else {
      if(temp > interleaver_size) {
        right = index - 1;
      }
      else {
        left = index + 1;
      }
    }
  }
  it_assert(!search, "lte_turbo_interleaver_sequence: The interleaver size is incorrect!");

  // Definitions of key parameters:
  int K = interleaver_size; // Interleaver size
  int f1_factor = f1_factors(index);
  int f2_factor = f2_factors(index);
  ivec I(K); //The interleaver sequence

  // Calculate the interleaver sequence:
  for(int i = 0; i < K; i++) {
    I(i) = static_cast<int>((static_cast<int64_t>(i) * f1_factor + static_cast<int64_t>(i) * i * f2_factor) % K);
  }

  return I;
}


} // namespace itpp
