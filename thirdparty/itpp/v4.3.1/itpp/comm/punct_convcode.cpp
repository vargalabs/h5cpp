/*!
 * \file
 * \brief Implementation of a Binary Punctured Convolutional Encoder class
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

#include <itpp/comm/punct_convcode.h>


namespace itpp
{

// --------------------- Punctured_Conv_Code --------------------------------

//------- Protected functions -----------------------
int Punctured_Convolutional_Code::weight(int state, int input, int time)
{
  int i, j, temp, out, w = 0, shiftreg = state;

  shiftreg = shiftreg | (int(input) << m);
  for (j = 0; j < n; j++) {
    if (puncture_matrix(j, time) == bin(1)) {
      out = 0;
      temp = shiftreg & gen_pol(j);
      for (i = 0; i < K; i++) {
        out ^= (temp & 1);
        temp = temp >> 1;
      }
      w += out;
    }
  }
  return w;
}

int Punctured_Convolutional_Code::weight_reverse(int state, int input, int time)
{
  int i, j, temp, out, w = 0, shiftreg = state;

  shiftreg = shiftreg | (int(input) << m);
  for (j = 0; j < n; j++) {
    if (puncture_matrix(j, Period - 1 - time) == bin(1)) {
      out = 0;
      temp = shiftreg & gen_pol_rev(j);
      for (i = 0; i < K; i++) {
        out ^= (temp & 1);
        temp = temp >> 1;
      }
      w += out;
    }
  }
  return w;
}

void Punctured_Convolutional_Code::weight(int state, int &w0, int &w1, int time)
{
  int i, j, temp, out, shiftreg = state;
  w0 = 0;
  w1 = 0;

  shiftreg = shiftreg | (1 << m);
  for (j = 0; j < n; j++) {
    if (puncture_matrix(j, time) == bin(1)) {
      out = 0;
      temp = shiftreg & gen_pol(j);
      for (i = 0; i < m; i++) {
        out ^= (temp & 1);
        temp = temp >> 1;
      }
      w0 += out;
      w1 += out ^(temp & 1);
    }
  }
}

void Punctured_Convolutional_Code::weight_reverse(int state, int &w0, int &w1, int time)
{
  int i, j, temp, out, shiftreg = state;
  w0 = 0;
  w1 = 0;

  shiftreg = shiftreg | (1 << m);
  for (j = 0; j < n; j++) {
    if (puncture_matrix(j, Period - 1 - time) == bin(1)) {
      out = 0;
      temp = shiftreg & gen_pol_rev(j);
      for (i = 0; i < m; i++) {
        out ^= (temp & 1);
        temp = temp >> 1;
      }
      w0 += out;
      w1 += out ^(temp & 1);
    }
  }
}

//------- Public functions -----------------------

void Punctured_Convolutional_Code::set_puncture_matrix(const bmat &pmatrix)
{
  it_error_if((pmatrix.rows() != n) || (pmatrix.cols() == 0), "Wrong size of puncture matrix");
  puncture_matrix = pmatrix;
  Period = puncture_matrix.cols();

  int p, j;
  total = 0;

  for (j = 0; j < n; j++) {
    for (p = 0; p < Period; p++)
      total = total + (int)(puncture_matrix(j, p));
  }
  rate = (double)Period / total;
}

void Punctured_Convolutional_Code::encode(const bvec &input, bvec &output)
{
  switch (cc_method) {
  case Trunc:
    encode_trunc(input, output);
    break;
  case Tail:
    encode_tail(input, output);
    break;
  case Tailbite:
    encode_tailbite(input, output);
    break;
  default:
    encode_tail(input, output);
    break;
  };
}

void Punctured_Convolutional_Code::encode_trunc(const bvec &input, bvec &output)
{
  Convolutional_Code::encode_trunc(input, output);

  int nn = 0, i, p = 0, j;

  for (i = 0; i < int(output.size() / n); i++) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, p) == bin(1)) {
        output(nn) = output(i * n + j);
        nn++;
      }
    }
    p = (p + 1) % Period;
  }
  output.set_size(nn, true);
}

void Punctured_Convolutional_Code::encode_tail(const bvec &input, bvec &output)
{
  Convolutional_Code::encode_tail(input, output);

  int nn = 0, i, p = 0, j;

  for (i = 0; i < int(output.size() / n); i++) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, p) == bin(1)) {
        output(nn) = output(i * n + j);
        nn++;
      }
    }
    p = (p + 1) % Period;
  }
  output.set_size(nn, true);
}

void Punctured_Convolutional_Code::encode_tailbite(const bvec &input, bvec &output)
{
  Convolutional_Code::encode_tailbite(input, output);

  int nn = 0, i, p = 0, j;

  for (i = 0; i < int(output.size() / n); i++) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, p) == bin(1)) {
        output(nn) = output(i * n + j);
        nn++;
      }
    }
    p = (p + 1) % Period;
  }
  output.set_size(nn, true);
}


// --------------- Hard-decision decoding is not implemented --------------------------------
void Punctured_Convolutional_Code::decode(const bvec &, bvec &)
{
  it_error("Punctured_Convolutional_Code::decode(bvec, bvec); hard-decision decoding is not implemented");
}

bvec Punctured_Convolutional_Code::decode(const bvec &)
{
  it_error("Punctured_Convolutional_Code::decode(bvec, bvec); hard-decision decoding is not implemented");
  return bvec();
}

/*
  Decode a block of encoded data using specified method
*/
void Punctured_Convolutional_Code::decode(const vec &received_signal, bvec &output)
{
  switch (cc_method) {
  case Trunc:
    decode_trunc(received_signal, output);
    break;
  case Tail:
    decode_tail(received_signal, output);
    break;
  case Tailbite:
    decode_tailbite(received_signal, output);
    break;
  default:
    decode_tail(received_signal, output);
    break;
  };
}


// Viterbi decoder using TruncLength (=5*K if not specified)
void Punctured_Convolutional_Code::decode_trunc(const vec &received_signal, bvec &output)
{
  int nn = 0, i = 0, p = received_signal.size() / total, j;

  int temp_size = p * Period * n;
  // number of punctured bits in a fraction of the puncture matrix
  // correcponding to the end of the received_signal vector
  p = received_signal.size() - p * total;
  // precise calculation of the number of unpunctured bits
  // (in the above fraction of the puncture matrix)
  while (p > 0) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, nn) == bin(1))
        p--;
    }
    nn++;
  }
  temp_size += n * nn;
  if (p != 0) {
    it_warning("Punctured_Convolutional_Code::decode(): Improper length of "
               "the received punctured block, dummy bits have been added");
  }

  vec temp(temp_size);
  nn = 0;
  j = 0;
  p = 0;

  while (nn < temp.size()) {
    if ((puncture_matrix(j, p) == bin(1)) && (i < received_signal.size())) {
      temp(nn) = received_signal(i);
      i++;
    }
    else { // insert dummy symbols with the same contribution for 0 and 1
      temp(nn) = 0;
    }

    nn++;
    j++;

    if (j == n) {
      j = 0;
      p = (p + 1) % Period;
    }
  } // while

  Convolutional_Code::decode_trunc(temp, output);
}

// Viterbi decoder using TruncLength (=5*K if not specified)
void Punctured_Convolutional_Code::decode_tail(const vec &received_signal, bvec &output)
{
  int nn = 0, i = 0, p = received_signal.size() / total, j;

  int temp_size = p * Period * n;
  // number of punctured bits in a fraction of the puncture matrix
  // correcponding to the end of the received_signal vector
  p = received_signal.size() - p * total;
  // precise calculation of the number of unpunctured bits
  // (in the above fraction of the puncture matrix)
  while (p > 0) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, nn) == bin(1))
        p--;
    }
    nn++;
  }
  temp_size += n * nn;
  if (p != 0) {
    it_warning("Punctured_Convolutional_Code::decode_tail(): Improper length "
               "of the received punctured block, dummy bits have been added");
  }

  vec temp(temp_size);

  nn = 0;
  j = 0;
  p = 0;

  while (nn < temp.size()) {
    if ((puncture_matrix(j, p) == bin(1)) && (i < received_signal.size())) {
      temp(nn) = received_signal(i);
      i++;
    }
    else { // insert dummy symbols with same contribution for 0 and 1
      temp(nn) = 0;
    }

    nn++;
    j++;

    if (j == n) {
      j = 0;
      p = (p + 1) % Period;
    }
  } // while

  Convolutional_Code::decode_tail(temp, output);
}

// Decode a block of encoded data where encode_tailbite has been used. Tries all start states.
void Punctured_Convolutional_Code::decode_tailbite(const vec &received_signal, bvec &output)
{
  int nn = 0, i = 0, p = received_signal.size() / total, j;

  int temp_size = p * Period * n;
  // number of punctured bits in a fraction of the puncture matrix
  // correcponding to the end of the received_signal vector
  p = received_signal.size() - p * total;
  // precise calculation of the number of unpunctured bits
  // (in the above fraction of the puncture matrix)
  while (p > 0) {
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, nn) == bin(1))
        p--;
    }
    nn++;
  }
  temp_size += n * nn;
  if (p != 0) {
    it_warning("Punctured_Convolutional_Code::decode_tailbite(): Improper "
               "length of the received punctured block, dummy bits have been "
               "added");
  }

  vec temp(temp_size);

  nn = 0;
  j = 0;
  p = 0;

  while (nn < temp.size()) {
    if ((puncture_matrix(j, p) == bin(1)) && (i < received_signal.size())) {
      temp(nn) = received_signal(i);
      i++;
    }
    else { // insert dummy symbols with same contribution for 0 and 1
      temp(nn) = 0;
    }

    nn++;
    j++;

    if (j == n) {
      j = 0;
      p = (p + 1) % Period;
    }
  } // while

  Convolutional_Code::decode_tailbite(temp, output);
}


/*
  Calculate the inverse sequence

  Assumes that encode_tail is used in the encoding process. Returns false if there is an
  error in the coded sequence (not a valid codeword). Does not check that the tail forces
  the encoder into the zeroth state.
*/
bool Punctured_Convolutional_Code::inverse_tail(const bvec coded_sequence, bvec &input)
{
  int state = 0, zero_state, one_state, zero_temp, one_temp, i, j, p = 0, prev_pos = 0, no_symbols;
  int block_length = 0;
  bvec zero_output(n), one_output(n), temp_output(n);

  block_length = coded_sequence.size() * Period / total - m;

  it_error_if(block_length <= 0, "The input sequence is to short");
  input.set_length(block_length, false); // not include the tail in the output

  p = 0;

  for (i = 0; i < block_length; i++) {
    zero_state = state;
    one_state = state | (1 << m);
    no_symbols = 0;
    for (j = 0; j < n; j++) {
      if (puncture_matrix(j, p) == bin(1)) {
        zero_temp = zero_state & gen_pol(j);
        one_temp = one_state & gen_pol(j);
        zero_output(no_symbols) = xor_int_table(zero_temp);
        one_output(no_symbols) = xor_int_table(one_temp);
        no_symbols++;
      }
    }
    if (coded_sequence.mid(prev_pos, no_symbols) == zero_output.left(no_symbols)) {
      input(i) = bin(0);
      state = zero_state >> 1;
    }
    else if (coded_sequence.mid(prev_pos, no_symbols) == one_output.left(no_symbols)) {
      input(i) = bin(1);
      state = one_state >> 1;
    }
    else
      return false;

    prev_pos += no_symbols;
    p = (p + 1) % Period;
  }

  return true;
}




/*
   It is possible to do this more efficiently by registrating all (inputs,states,Periods)
   that has zero metric (just and with the generators). After that build paths from a input=1 state.
*/
bool Punctured_Convolutional_Code::catastrophic(void)
{
  int max_stack_size = 50000;
  ivec S_stack(max_stack_size), t_stack(max_stack_size);
  int start, S, W0, W1, S0, S1, pos, t = 0;
  int stack_pos = -1;

  for (pos = 0; pos < Period; pos++) { // test all start positions
    for (start = 0; start < (1 << m); start++) {
      stack_pos = -1;
      S = start;
      t = 0;

    node1:
      if (t > (1 << m) * Period) { return true; }
      S0 = next_state(S, 0);
      S1 = next_state(S, 1);
      weight(S, W0, W1, (pos + t) % Period);
      if (W1 > 0) { goto node0; }
      if ((W0 == 0) && (S0 == start) && (((pos + t + 1) % Period) == pos)) { return true; }
      if ((W0 == 0) && (S0 == 0) && (S != 0)) { return true; }
      if ((S0 != 0) && (W0 == 0)) {
        stack_pos++;
        if (stack_pos >= max_stack_size) {
          max_stack_size = 2 * max_stack_size;
          S_stack.set_size(max_stack_size, true);
          t_stack.set_size(max_stack_size, true);
        }
        S_stack(stack_pos) = S0;
        t_stack(stack_pos) = t + 1;
      }
      if ((W1 == 0) && (S1 == start) && (((pos + t + 1) % Period) == pos)) { return true; }
      S = S1;
      t++;
      goto node1;

    node0:
      if (W0 > 0) goto stack;
      if ((W0 == 0) && (S0 == start) && (((pos + t + 1) % Period) == pos)) { return true; }
      if ((W0 == 0) && (S0 == 0) && (S != 0)) { return true; }
      if (S0 != 0) {
        S = S0;
        t++;
        goto node1;
      }
      else {
        goto stack;
      }

    stack:
      if (stack_pos >= 0) {
        S = S_stack(stack_pos);
        t = t_stack(stack_pos);
        stack_pos--;
        goto node1;
      }
    }
  }
  return false;
}

void Punctured_Convolutional_Code::distance_profile(ivec &dist_prof, int start_time, int dmax, bool reverse)
{
  int max_stack_size = 50000;
  ivec S_stack(max_stack_size), W_stack(max_stack_size), t_stack(max_stack_size);

  int stack_pos = -1, t, S, W, W0, w0, w1;


  dist_prof.zeros();
  dist_prof += dmax; // just select a big number!
  if (reverse)
    W = weight_reverse(0, 1, start_time);
  else
    W = weight(0, 1, start_time);

  S = next_state(0, 1); // start from zero state and a one input
  t = 0;
  dist_prof(0) = W;

node1:
  if (reverse)
    weight_reverse(S, w0, w1, (start_time + t + 1) % Period);
  else
    weight(S, w0, w1, (start_time + t + 1) % Period);

  if (t < m) {
    W0 = W + w0;
    if (W0 < dist_prof(m)) { // store node0 (S, W0, and t+1)
      stack_pos++;
      if (stack_pos >= max_stack_size) {
        max_stack_size = 2 * max_stack_size;
        S_stack.set_size(max_stack_size, true);
        W_stack.set_size(max_stack_size, true);
        t_stack.set_size(max_stack_size, true);
      }
      S_stack(stack_pos) = next_state(S, 0);
      W_stack(stack_pos) = W0;
      t_stack(stack_pos) = t + 1;
    }
  }
  else goto stack;

  W += w1;
  if (W > dist_prof(m))
    goto stack;

  t++;
  S = next_state(S, 1);

  if (W < dist_prof(t))
    dist_prof(t) = W;

  if (t == m) goto stack;
  goto node1;


stack:
  if (stack_pos >= 0) {
    // get root node from stack
    S = S_stack(stack_pos);
    W = W_stack(stack_pos);
    t = t_stack(stack_pos);
    stack_pos--;

    if (W < dist_prof(t))
      dist_prof(t) = W;

    if (t == m) goto stack;

    goto node1;
  }
}

int Punctured_Convolutional_Code::fast(Array<ivec> &spectrum, int start_time, int dfree, int no_terms,
                                       int d_best_so_far, bool test_catastrophic)
{
  int cat_treshold = (1 << m) * Period;
  int i, j, t = 0;
  ivec dist_prof(K), dist_prof_rev(K), dist_prof_temp(K);

  //calculate distance profile
  distance_profile(dist_prof, start_time, dfree);
  distance_profile(dist_prof_rev, 0, dfree, true); // for the reverse code

  int dist_prof_rev0 = dist_prof_rev(0);

  bool reverse = true; // true = use reverse dist_prof

  // choose the lowest dist_prof over all periods
  for (i = 1; i < Period; i++) {
    distance_profile(dist_prof_temp, i, dfree, true);
    // switch if new is lower
    for (j = 0; j < K; j++) {
      if (dist_prof_temp(j) < dist_prof_rev(j)) {
        dist_prof_rev(j) = dist_prof_temp(j);
      }
    }
  }

  dist_prof_rev0 = dist_prof(0);
  dist_prof = dist_prof_rev;

  int d = dfree + no_terms - 1;
  int max_stack_size = 50000;
  ivec S_stack(max_stack_size), W_stack(max_stack_size), b_stack(max_stack_size);
  ivec c_stack(max_stack_size), t_stack(max_stack_size);
  int stack_pos = -1;

  // F1.
  int mf = 1, b = 1;
  spectrum.set_size(2);
  spectrum(0).set_size(dfree + no_terms, 0); // ad
  spectrum(1).set_size(dfree + no_terms, 0); // cd
  spectrum(0).zeros();
  spectrum(1).zeros();
  int S, S0, S1, W0, W1, w0, w1, c = 0;
  S = next_state(0, 1); // start in zero state and one input
  int W = d - dist_prof_rev0;
  t = 1;

F2:
  S0 = next_state(S, 0);
  S1 = next_state(S, 1);

  if (reverse) {
    weight(S, w0, w1, (start_time + t) % Period);
  }
  else {
    weight_reverse(S, w0, w1, (start_time + t) % Period);
  }
  W0 = W - w0;
  W1 = W - w1;

  if (mf < m) goto F6;

  //F3:
  if (W0 >= 0) {
    spectrum(0)(d - W0)++;
    spectrum(1)(d - W0) += b;
  }
  //Test on d_best_so_far
  if ((d - W0) < d_best_so_far) { return -1; }

F4:
  if ((W1 < dist_prof(m - 1)) || (W < dist_prof(m))) goto F5;
  // select node 1
  if (test_catastrophic && (W == W1)) {
    c++;
    if (c > cat_treshold)
      return 0;
  }
  else {
    c = 0;
  }

  W = W1;
  S = S1;
  t++;
  mf = 1;
  b++;
  goto F2;

F5:
  if (stack_pos == -1) goto end;
  // get node from stack
  S = S_stack(stack_pos);
  W = W_stack(stack_pos);
  b = b_stack(stack_pos);
  c = c_stack(stack_pos);
  t = t_stack(stack_pos);
  stack_pos--;
  mf = 1;
  goto F2;

F6:
  if (W0 < dist_prof(m - mf - 1)) goto F4;

  //F7:
  if ((W1 >= dist_prof(m - 1)) && (W >= dist_prof(m))) {
    // save node 1
    if (test_catastrophic && (stack_pos > 40000))
      return 0;

    stack_pos++;
    if (stack_pos >= max_stack_size) {
      max_stack_size = 2 * max_stack_size;
      S_stack.set_size(max_stack_size, true);
      W_stack.set_size(max_stack_size, true);
      b_stack.set_size(max_stack_size, true);
      c_stack.set_size(max_stack_size, true);
      t_stack.set_size(max_stack_size, true);
    }
    S_stack(stack_pos) = S1;
    W_stack(stack_pos) = W1;
    b_stack(stack_pos) = b + 1; // because gone to node 1
    c_stack(stack_pos) = c;
    t_stack(stack_pos) = t + 1;
  }
  // select node 0
  S = S0;
  t++;
  if (test_catastrophic && (W == W0)) {
    c++;
    if (c > cat_treshold)
      return false;
  }
  else {
    c = 0;
  }

  W = W0;
  mf++;
  goto F2;

end:
  return 1;
}

void Punctured_Convolutional_Code::calculate_spectrum(Array<ivec> &spectrum, int dmax, int no_terms)
{
  Array<ivec> temp_spectra(2);
  spectrum.set_size(2);
  spectrum(0).set_size(dmax + no_terms, false);
  spectrum(1).set_size(dmax + no_terms, false);
  spectrum(0).zeros();
  spectrum(1).zeros();

  for (int pos = 0; pos < Period; pos++) {
    calculate_spectrum(temp_spectra, pos, dmax, no_terms);
    spectrum(0) += temp_spectra(0);
    spectrum(1) += temp_spectra(1);
  }
}

void Punctured_Convolutional_Code::calculate_spectrum(Array<ivec> &spectrum, int time, int dmax, int no_terms, int block_length)
{
  imat Ad_states(1 << (K - 1), dmax + no_terms), Cd_states(1 << m, dmax + no_terms);
  imat Ad_temp(1 << (K - 1), dmax + no_terms), Cd_temp(1 << m, dmax + no_terms);
  ivec mindist(1 << (K - 1)), mindist_temp(1 << m);

  spectrum.set_size(2);
  spectrum(0).set_size(dmax + no_terms, false);
  spectrum(1).set_size(dmax + no_terms, false);
  spectrum(0).zeros();
  spectrum(1).zeros();
  Ad_states.zeros();
  Cd_states.zeros();
  mindist.zeros();
  int wmax = dmax + no_terms;
  ivec visited_states(1 << m), visited_states_temp(1 << m);
  bool proceede, expand_s1;
  int t, d, w0, w1, s, s0, s1 = 0, s_start;

  // initial phase. Evolv trellis up to time K.
  visited_states.zeros(); // 0 = false

  s1 = next_state(0, 1);   // Start in state 0 and 1 input
  visited_states(s1) = 1;  // 1 = true.
  w1 = weight(0, 1, time);
  Ad_states(s1, w1) = 1;
  Cd_states(s1, w1) = 1;
  mindist(s1) = w1;

  if (block_length > 0) {
    s0 = next_state(0, 0);
    visited_states(s0) = 1;  // 1 = true. Expand also the zero-path
    w0 = weight(0, 0, time);
    Ad_states(s0, w0) = 1;
    Cd_states(s0, w0) = 0;
    mindist(s0) = w0;
    s_start = 0;
  }
  else {
    s_start = 1;
  }

  t = 1;
  proceede = true;
  while (proceede) {
    proceede = false;
    Ad_temp.zeros();
    Cd_temp.zeros();
    mindist_temp.zeros();
    visited_states_temp.zeros(); //false

    for (s = s_start; s < (1 << m); s++) {

      if (visited_states(s) == 1) {
        if ((mindist(s) >= 0) && (mindist(s) < wmax)) {
          proceede = true;
          weight(s, w0, w1, (time + t) % Period);

          s0 = next_state(s, 0);
          for (d = mindist(s); d < (wmax - w0); d++) {
            Ad_temp(s0, d + w0) += Ad_states(s, d);
            Cd_temp(s0, d + w0) += Cd_states(s, d);
          }
          if (visited_states_temp(s0) == 0) { mindist_temp(s0) = mindist(s) + w0; }
          else { mindist_temp(s0) = ((mindist(s) + w0) < mindist_temp(s0)) ? mindist(s) + w0 :  mindist_temp(s0); }
          visited_states_temp(s0) = 1; //true

          expand_s1 = false;
          if ((block_length == 0) || (t < (block_length - (K - 1)))) {
            expand_s1 = true;
          }

          if (expand_s1) {
            s1 = next_state(s, 1);
            for (d = mindist(s); d < (wmax - w1); d++) {
              Ad_temp(s1, d + w1) += Ad_states(s, d);
              Cd_temp(s1, d + w1) += Cd_states(s, d) + Ad_states(s, d);
            }
            if (visited_states_temp(s1) == 0) { mindist_temp(s1) = mindist(s) + w1; }
            else { mindist_temp(s1) = ((mindist(s) + w1) < mindist_temp(s1)) ? mindist(s) + w1 :  mindist_temp(s1); }
            visited_states_temp(s1) = 1; //true
          }
        }
      }
    }

    Ad_states = Ad_temp;
    Cd_states = Cd_temp;
    if (block_length == 0) {
      spectrum(0) += Ad_temp.get_row(0);
      spectrum(1) += Cd_temp.get_row(0);
    }
    visited_states = visited_states_temp;
    mindist = elem_mult(mindist_temp, visited_states);
    t++;
    if ((t > block_length) && (block_length > 0)) { proceede = false; }
  }

  if (block_length > 0) {
    spectrum(0) = Ad_states.get_row(0);
    spectrum(1) = Cd_states.get_row(0);
  }

}

} // namespace itpp
