/*!
 * \file
 * \brief Implementation of a binary convolutional encoder class
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

#include <itpp/comm/convcode.h>
#include <itpp/base/binary.h>
#include <itpp/base/matfunc.h>
#include <limits>

namespace itpp
{

// ----------------- Protected functions -----------------------------

/*
  The weight of the transition from given state with the input given
*/
int Convolutional_Code::weight(const int state, const int input)
{
  int i, j, temp, out, w = 0, shiftreg = state;

  shiftreg = shiftreg | (input << m);
  for (j = 0; j < n; j++) {
    out = 0;
    temp = shiftreg & gen_pol(j);
    for (i = 0; i < K; i++) {
      out ^= (temp & 1);
      temp = temp >> 1;
    }
    w += out;
    //w += weight_int_table(temp);
  }
  return w;
}

/*
  The weight (of the reverse code) of the transition from given state with
  the input given
*/
int Convolutional_Code::weight_reverse(const int state, const int input)
{
  int i, j, temp, out, w = 0, shiftreg = state;

  shiftreg = shiftreg | (input << m);
  for (j = 0; j < n; j++) {
    out = 0;
    temp = shiftreg & gen_pol_rev(j);
    for (i = 0; i < K; i++) {
      out ^= (temp & 1);
      temp = temp >> 1;
    }
    w += out;
    //w += weight_int_table(temp);
  }
  return w;
}

/*
  The weight of the two paths (input 0 or 1) from given state
*/
void Convolutional_Code::weight(const int state, int &w0, int &w1)
{
  int i, j, temp, out, shiftreg = state;
  w0 = 0;
  w1 = 0;

  shiftreg = shiftreg | (1 << m);
  for (j = 0; j < n; j++) {
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

/*
  The weight (of the reverse code) of the two paths (input 0 or 1) from
  given state
*/
void Convolutional_Code::weight_reverse(const int state, int &w0, int &w1)
{
  int i, j, temp, out, shiftreg = state;
  w0 = 0;
  w1 = 0;

  shiftreg = shiftreg | (1 << m);
  for (j = 0; j < n; j++) {
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

/*
  Output on transition (backwards) with input from state
*/
bvec Convolutional_Code::output_reverse(const int state, const int input)
{
  int j, temp, temp_state;

  bvec output(n);

  temp_state = (state << 1) | input;
  for (j = 0; j < n; j++) {
    temp = temp_state & gen_pol(j);
    output(j) = xor_int_table(temp);
  }

  return output;
}

/*
  Output on transition (backwards) with input from state
*/
void Convolutional_Code::output_reverse(const int state, bvec &zero_output,
                                        bvec &one_output)
{
  int j, temp, temp_state;
  bin one_bit;

  temp_state = (state << 1) | 1;
  for (j = 0; j < n; j++) {
    temp = temp_state & gen_pol(j);
    one_bit = temp & 1;
    temp = temp >> 1;
    one_output(j) = xor_int_table(temp) ^ one_bit;
    zero_output(j) = xor_int_table(temp);
  }
}

/*
  Output on transition (backwards) with input from state
*/
void Convolutional_Code::output_reverse(const int state, int &zero_output,
                                        int &one_output)
{
  int j, temp, temp_state;
  bin one_bit;

  zero_output = 0, one_output = 0;
  temp_state = (state << 1) | 1;
  for (j = 0; j < n; j++) {
    temp = temp_state & gen_pol(j);
    one_bit = temp & 1;
    temp = temp >> 1;

    one_output = (one_output << 1) | int(xor_int_table(temp) ^ one_bit);
    zero_output = (zero_output << 1) | int(xor_int_table(temp));
  }
}

/*
  Output on transition (backwards) with input from state
*/
void Convolutional_Code::calc_metric_reverse(int state,
    const vec &rx_codeword,
    double &zero_metric,
    double &one_metric)
{
  int temp;
  bin one_bit;
  one_metric = zero_metric = 0;

  int temp_state = (state << 1) | 1;
  for (int j = 0; j < n; j++) {
    temp = temp_state & gen_pol(j);
    one_bit = temp & 1;
    temp >>= 1;
    one_metric += (2 * static_cast<int>(xor_int_table(temp) ^ one_bit) - 1)
                  * rx_codeword(j);
    zero_metric += (2 * static_cast<int>(xor_int_table(temp)) - 1)
                   * rx_codeword(j);
  }
}


// calculates metrics for all codewords (2^n of them) in natural order
void Convolutional_Code::calc_metric(const vec &rx_codeword,
                                     vec &delta_metrics)
{
  int no_outputs = pow2i(n), no_loop = pow2i(n - 1), mask = no_outputs - 1,
                                       temp;
  delta_metrics.set_size(no_outputs, false);

  if (no_outputs <= no_states) {
    for (int i = 0; i < no_loop; i++) { // all input possibilities
      delta_metrics(i) = 0;
      temp = i;
      for (int j = n - 1; j >= 0; j--) {
        if (temp & 1)
          delta_metrics(i) += rx_codeword(j);
        else
          delta_metrics(i) -= rx_codeword(j);
        temp >>= 1;
      }
      delta_metrics(i ^ mask) = -delta_metrics(i); // the inverse codeword
    }
  }
  else {
    double zero_metric, one_metric;
    int zero_output, one_output, temp_state;
    bin one_bit;
    for (int s = 0; s < no_states; s++) { // all states
      zero_metric = 0, one_metric = 0;
      zero_output = 0, one_output = 0;
      temp_state = (s << 1) | 1;
      for (int j = 0; j < n; j++) {
        temp = temp_state & gen_pol(j);
        one_bit = temp & 1;
        temp >>= 1;
        if (xor_int_table(temp)) {
          zero_metric += rx_codeword(j);
        }
        else {
          zero_metric -= rx_codeword(j);
        }
        if (static_cast<int>(xor_int_table(temp)^one_bit)) {
        	one_metric += rx_codeword(j);
        } else {
        	one_metric -= rx_codeword(j);
        }
        one_output = (one_output << 1)
                     | static_cast<int>(xor_int_table(temp) ^ one_bit);
        zero_output = (zero_output << 1)
                      | static_cast<int>(xor_int_table(temp));
      }
      delta_metrics(zero_output) = zero_metric;
      delta_metrics(one_output) = one_metric;
    }
  }
}

//! \cond

// MFD codes R=1/2
int Conv_Code_MFD_2[15][2] = {
  {0, 0},
  {0, 0},
  {0, 0},
  {05, 07},
  {015, 017},
  {023, 035},
  {053, 075},
  {0133, 0171},
  {0247, 0371},
  {0561, 0753},
  {01167, 01545},
  {02335, 03661},
  {04335, 05723},
  {010533, 017661},
  {021675, 027123},
};

// MFD codes R=1/3
int Conv_Code_MFD_3[15][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0},
  {05, 07, 07},
  {013, 015, 017},
  {025, 033, 037},
  {047, 053, 075},
  {0133, 0145, 0175},
  {0225, 0331, 0367},
  {0557, 0663, 0711},
  {0117, 01365, 01633},
  {02353, 02671, 03175},
  {04767, 05723, 06265},
  {010533, 010675, 017661},
  {021645, 035661, 037133},
};

// MFD codes R=1/4
int Conv_Code_MFD_4[15][4] = {
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {05, 07, 07, 07},
  {013, 015, 015, 017},
  {025, 027, 033, 037},
  {053, 067, 071, 075},
  {0135, 0135, 0147, 0163},
  {0235, 0275, 0313, 0357},
  {0463, 0535, 0733, 0745},
  {0117, 01365, 01633, 01653},
  {02327, 02353, 02671, 03175},
  {04767, 05723, 06265, 07455},
  {011145, 012477, 015573, 016727},
  {021113, 023175, 035527, 035537},
};


// MFD codes R=1/5
int Conv_Code_MFD_5[9][5] = {
  {0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0},
  {07, 07, 07, 05, 05},
  {017, 017, 013, 015, 015},
  {037, 027, 033, 025, 035},
  {075, 071, 073, 065, 057},
  {0175, 0131, 0135, 0135, 0147},
  {0257, 0233, 0323, 0271, 0357},
};

// MFD codes R=1/6
int Conv_Code_MFD_6[9][6] = {
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0},
  {07, 07, 07, 07, 05, 05},
  {017, 017, 013, 013, 015, 015},
  {037, 035, 027, 033, 025, 035},
  {073, 075, 055, 065, 047, 057},
  {0173, 0151, 0135, 0135, 0163, 0137},
  {0253, 0375, 0331, 0235, 0313, 0357},
};

// MFD codes R=1/7
int Conv_Code_MFD_7[9][7] = {
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {07, 07, 07, 07, 05, 05, 05},
  {017, 017, 013, 013, 013, 015, 015},
  {035, 027, 025, 027, 033, 035, 037},
  {053, 075, 065, 075, 047, 067, 057},
  {0165, 0145, 0173, 0135, 0135, 0147, 0137},
  {0275, 0253, 0375, 0331, 0235, 0313, 0357},
};

// MFD codes R=1/8
int Conv_Code_MFD_8[9][8] = {
  {0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0},
  {07, 07, 05, 05, 05, 07, 07, 07},
  {017, 017, 013, 013, 013, 015, 015, 017},
  {037, 033, 025, 025, 035, 033, 027, 037},
  {057, 073, 051, 065, 075, 047, 067, 057},
  {0153, 0111, 0165, 0173, 0135, 0135, 0147, 0137},
  {0275, 0275, 0253, 0371, 0331, 0235, 0313, 0357},
};

int maxK_Conv_Code_MFD[9] = {0, 0, 14, 14, 14, 8, 8, 8, 8};

void get_MFD_gen_pol(int n, int K, ivec & gen)
{
  gen.set_size(n);

  switch (n) {
  case 2: // Rate 1/2
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[2], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_2[K][0];
    gen(1) = Conv_Code_MFD_2[K][1];
    break;
  case 3: // Rate 1/3
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[3], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_3[K][0];
    gen(1) = Conv_Code_MFD_3[K][1];
    gen(2) = Conv_Code_MFD_3[K][2];
    break;
  case 4: // Rate 1/4
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[4], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_4[K][0];
    gen(1) = Conv_Code_MFD_4[K][1];
    gen(2) = Conv_Code_MFD_4[K][2];
    gen(3) = Conv_Code_MFD_4[K][3];
    break;
  case 5: // Rate 1/5
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[5], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_5[K][0];
    gen(1) = Conv_Code_MFD_5[K][1];
    gen(2) = Conv_Code_MFD_5[K][2];
    gen(3) = Conv_Code_MFD_5[K][3];
    gen(4) = Conv_Code_MFD_5[K][4];
    break;
  case 6: // Rate 1/6
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[6], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_6[K][0];
    gen(1) = Conv_Code_MFD_6[K][1];
    gen(2) = Conv_Code_MFD_6[K][2];
    gen(3) = Conv_Code_MFD_6[K][3];
    gen(4) = Conv_Code_MFD_6[K][4];
    gen(5) = Conv_Code_MFD_6[K][5];
    break;
  case 7: // Rate 1/7
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[7], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_7[K][0];
    gen(1) = Conv_Code_MFD_7[K][1];
    gen(2) = Conv_Code_MFD_7[K][2];
    gen(3) = Conv_Code_MFD_7[K][3];
    gen(4) = Conv_Code_MFD_7[K][4];
    gen(5) = Conv_Code_MFD_7[K][5];
    gen(6) = Conv_Code_MFD_7[K][6];
    break;
  case 8: // Rate 1/8
    it_assert(K >= 3 && K <= maxK_Conv_Code_MFD[8], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_MFD_8[K][0];
    gen(1) = Conv_Code_MFD_8[K][1];
    gen(2) = Conv_Code_MFD_8[K][2];
    gen(3) = Conv_Code_MFD_8[K][3];
    gen(4) = Conv_Code_MFD_8[K][4];
    gen(5) = Conv_Code_MFD_8[K][5];
    gen(6) = Conv_Code_MFD_8[K][6];
    gen(7) = Conv_Code_MFD_8[K][7];
    break;
  default: // wrong rate
    it_assert(false, "This convolutional code doesn't exist in the tables");
  }
}

// ODS codes R=1/2
int Conv_Code_ODS_2[17][2] = {
  {0, 0},
  {0, 0},
  {0, 0},
  {05, 07},
  {015, 017},
  {023, 035},
  {053, 075},
  {0133, 0171},
  {0247, 0371},
  {0561, 0753},
  {01151, 01753},
  {03345, 03613},
  {05261, 07173},
  {012767, 016461},
  {027251, 037363},
  {063057, 044735},
  {0126723, 0152711},
};

// ODS codes R=1/3
int Conv_Code_ODS_3[14][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0},
  {05, 07, 07},
  {013, 015, 017},
  {025, 033, 037},
  {047, 053, 075},
  {0133, 0165, 0171},
  {0225, 0331, 0367},
  {0575, 0623, 0727},
  {01233, 01375, 01671},
  {02335, 02531, 03477},
  {05745, 06471, 07553},
  {013261, 015167, 017451},
};

// ODS codes R=1/4
int Conv_Code_ODS_4[12][4] = {
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {05, 05, 07, 07},
  {013, 015, 015, 017},
  {025, 027, 033, 037},
  {051, 055, 067, 077},
  {0117, 0127, 0155, 0171},
  {0231, 0273, 0327, 0375},
  {0473, 0513, 0671, 0765},
  {01173, 01325, 01467, 01751},
  {02565, 02747, 03311, 03723},
};

int maxK_Conv_Code_ODS[5] = {0, 0, 16, 13, 11};

void get_ODS_gen_pol(int n, int K, ivec & gen)
{
  gen.set_size(n);

  switch (n) {
  case 2: // Rate 1/2
    it_assert(K >= 3 && K <= maxK_Conv_Code_ODS[2], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_ODS_2[K][0];
    gen(1) = Conv_Code_ODS_2[K][1];
    break;
  case 3: // Rate 1/3
    it_assert(K >= 3 && K <= maxK_Conv_Code_ODS[3], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_ODS_3[K][0];
    gen(1) = Conv_Code_ODS_3[K][1];
    gen(2) = Conv_Code_ODS_3[K][2];
    break;
  case 4: // Rate 1/4
    it_assert(K >= 3 && K <= maxK_Conv_Code_ODS[4], "This convolutional code doesn't exist in the tables");
    gen(0) = Conv_Code_ODS_4[K][0];
    gen(1) = Conv_Code_ODS_4[K][1];
    gen(2) = Conv_Code_ODS_4[K][2];
    gen(3) = Conv_Code_ODS_4[K][3];
    break;
  default: // wrong rate
    it_assert(false, "This convolutional code doesn't exist in the tables");
  }
}

//! \endcond

// --------------- Public functions -------------------------

void Convolutional_Code::set_code(const CONVOLUTIONAL_CODE_TYPE type_of_code,
                                  int inverse_rate, int constraint_length)
{
  ivec gen;

  if (type_of_code == MFD)
    get_MFD_gen_pol(inverse_rate, constraint_length, gen);
  else if (type_of_code == ODS)
    get_ODS_gen_pol(inverse_rate, constraint_length, gen);
  else
    it_assert(false, "This convolutional code doesn't exist in the tables");

  set_generator_polynomials(gen, constraint_length);
}

/*
  Set generator polynomials. Given in Proakis integer form
*/
void Convolutional_Code::set_generator_polynomials(const ivec &gen,
    int constraint_length)
{
  it_error_if(constraint_length <= 0, "Convolutional_Code::set_generator_polynomials(): Constraint length out of range");
  gen_pol = gen;
  n = gen.size();
  it_error_if(n <= 0, "Convolutional_Code::set_generator_polynomials(): Invalid code rate");

  // Generate table look-up of weight of integers of size K bits
  if (constraint_length != K) {
    K = constraint_length;
    xor_int_table.set_size(pow2i(K), false);
    for (int i = 0; i < pow2i(K); i++) {
      xor_int_table(i) = (weight_int(K, i) & 1);
    }
  }

  trunc_length = 5 * K;
  m = K - 1;
  no_states = pow2i(m);
  gen_pol_rev.set_size(n, false);
  rate = 1.0 / n;

  for (int i = 0; i < n; i++) {
    gen_pol_rev(i) = reverse_int(K, gen_pol(i));
  }

  int zero_output, one_output;
  output_reverse_int.set_size(no_states, 2, false);

  for (int i = 0; i < no_states; i++) {
    output_reverse(i, zero_output, one_output);
    output_reverse_int(i, 0) = zero_output;
    output_reverse_int(i, 1) = one_output;
  }

  // initialise memory structures
  visited_state.set_size(no_states);
  visited_state = false;
  visited_state(start_state) = true; // mark start state

  sum_metric.set_size(no_states);
  sum_metric.clear();

  trunc_ptr = 0;
  trunc_state = 0;

}

// Reset encoder and decoder states
void Convolutional_Code::reset()
{
  init_encoder();

  visited_state = false;
  visited_state(start_state) = true; // mark start state

  sum_metric.clear();

  trunc_ptr = 0;
  trunc_state = 0;
}


/*
  Encode a binary vector of inputs using specified method
*/
void Convolutional_Code::encode(const bvec &input, bvec &output)
{
  switch (cc_method) {
  case Trunc:
    encode_trunc(input, output);
    break;
  case Tailbite:
    encode_tailbite(input, output);
    break;
  case Tail:
  default:
    encode_tail(input, output);
    break;
  };
}

/*
  Encode a binary vector of inputs starting from state set by the
  set_state function
*/
void Convolutional_Code::encode_trunc(const bvec &input, bvec &output)
{
  int temp;
  output.set_size(input.size() * n, false);

  for (int i = 0; i < input.size(); i++) {
    encoder_state |=  static_cast<int>(input(i)) << m;
    for (int j = 0; j < n; j++) {
      temp = encoder_state & gen_pol(j);
      output(i * n + j) = xor_int_table(temp);
    }
    encoder_state >>= 1;
  }
}

/*
  Encode a binary vector of inputs starting from zero state and also adds
  a tail of K-1 zeros to force the encoder into the zero state. Well
  suited for packet transmission.
*/
void Convolutional_Code::encode_tail(const bvec &input, bvec &output)
{
  int temp;
  output.set_size((input.size() + m) * n, false);

  // always start from state 0
  encoder_state = 0;

  for (int i = 0; i < input.size(); i++) {
    encoder_state |=  static_cast<int>(input(i)) << m;
    for (int j = 0; j < n; j++) {
      temp = encoder_state & gen_pol(j);
      output(i * n + j) = xor_int_table(temp);
    }
    encoder_state >>= 1;
  }

  // add tail of m = K-1 zeros
  for (int i = input.size(); i < input.size() + m; i++) {
    for (int j = 0; j < n; j++) {
      temp = encoder_state & gen_pol(j);
      output(i * n + j) = xor_int_table(temp);
    }
    encoder_state >>= 1;
  }
}

/*
  Encode a binary vector of inputs starting using tail biting
*/
void Convolutional_Code::encode_tailbite(const bvec &input, bvec &output)
{
  int temp;
  output.set_size(input.size() * n, false);

  // Set the start state equal to the end state:
  encoder_state = 0;
  bvec last_bits = input.right(m);
  for (int i = 0; i < m; i++) {
    encoder_state |= static_cast<int>(last_bits(i)) << m;
    encoder_state >>= 1;
  }

  for (int i = 0; i < input.size(); i++) {
    encoder_state |= static_cast<int>(input(i)) << m;
    for (int j = 0; j < n; j++) {
      temp = encoder_state & gen_pol(j);
      output(i * n + j) = xor_int_table(temp);
    }
    encoder_state >>= 1;
  }
}

/*
  Encode a binary bit starting from the internal encoder state.
  To initialize the encoder state use set_start_state() and init_encoder()
*/
void Convolutional_Code::encode_bit(const bin &input, bvec &output)
{
  int temp;
  output.set_size(n, false);

  encoder_state |= static_cast<int>(input) << m;
  for (int j = 0; j < n; j++) {
    temp = encoder_state & gen_pol(j);
    output(j) = xor_int_table(temp);
  }
  encoder_state >>= 1;
}


// --------------- Hard-decision decoding is not implemented -----------------

void Convolutional_Code::decode(const bvec &, bvec &)
{
  it_error("Convolutional_Code::decode(): Hard-decision decoding not implemented");
}

bvec Convolutional_Code::decode(const bvec &)
{
  it_error("Convolutional_Code::decode(): Hard-decision decoding not implemented");
  return bvec();
}


/*
  Decode a block of encoded data using specified method
*/
void Convolutional_Code::decode(const vec &received_signal, bvec &output)
{
  switch (cc_method) {
  case Trunc:
    decode_trunc(received_signal, output);
    break;
  case Tailbite:
    decode_tailbite(received_signal, output);
    break;
  case Tail:
  default:
    decode_tail(received_signal, output);
    break;
  };
}

/*
  Decode a block of encoded data where encode_tail has been used.
  Thus is assumes a decoder start state of zero and that a tail of
  K-1 zeros has been added. No memory truncation.
*/
void Convolutional_Code::decode_tail(const vec &received_signal, bvec &output)
{
  int block_length = received_signal.size() / n; // no input symbols
  it_error_if(block_length - m <= 0,
              "Convolutional_Code::decode_tail(): Input sequence to short");
  int S0, S1;
  vec temp_sum_metric(no_states), temp_rec(n), delta_metrics;
  Array<bool> temp_visited_state(no_states);
  double temp_metric_zero, temp_metric_one;

  path_memory.set_size(no_states, block_length, false);
  output.set_size(block_length - m, false);    // no tail in the output

  // clear visited states
  visited_state = false;
  temp_visited_state = visited_state;
  visited_state(0) = true; // starts in the zero state

  // clear accumulated metrics
  sum_metric.clear();

  // evolve up to m where not all states visited
  for (int l = 0; l < m; l++) { // all transitions including the tail
    temp_rec = received_signal.mid(l * n, n);

    // calculate all metrics for all codewords at the same time
    calc_metric(temp_rec, delta_metrics);

    for (int s = 0; s < no_states; s++) { // all states
      // S0 and S1 are the states that expanded end at state s
      previous_state(s, S0, S1);
      if (visited_state(S0)) { // expand trellis if state S0 is visited
        temp_metric_zero = sum_metric(S0)
                           + delta_metrics(output_reverse_int(s, 0));
        temp_visited_state(s) = true;
      }
      else {
        temp_metric_zero = std::numeric_limits<double>::max();
      }
      if (visited_state(S1)) { // expand trellis if state S0 is visited
        temp_metric_one = sum_metric(S1)
                          + delta_metrics(output_reverse_int(s, 1));
        temp_visited_state(s) = true;
      }
      else {
        temp_metric_one = std::numeric_limits<double>::max();
      }
      if (temp_metric_zero < temp_metric_one) { // path zero survives
        temp_sum_metric(s) = temp_metric_zero;
        path_memory(s, l) = 0;
      }
      else { // path one survives
        temp_sum_metric(s) = temp_metric_one;
        path_memory(s, l) = 1;
      }
    } // all states, s
    sum_metric = temp_sum_metric;
    visited_state = temp_visited_state;
  } // all transitions, l

  // evolve from m and to the end of the block
  for (int l = m; l < block_length; l++) { // all transitions except the tail
    temp_rec = received_signal.mid(l * n, n);

    // calculate all metrics for all codewords at the same time
    calc_metric(temp_rec, delta_metrics);

    for (int s = 0; s < no_states; s++) { // all states
      // S0 and S1 are the states that expanded end at state s
      previous_state(s, S0, S1);
      temp_metric_zero = sum_metric(S0)
                         + delta_metrics(output_reverse_int(s, 0));
      temp_metric_one = sum_metric(S1)
                        + delta_metrics(output_reverse_int(s, 1));
      if (temp_metric_zero < temp_metric_one) { // path zero survives
        temp_sum_metric(s) = temp_metric_zero;
        path_memory(s, l) = 0;
      }
      else { // path one survives
        temp_sum_metric(s) = temp_metric_one;
        path_memory(s, l) = 1;
      }
    } // all states, s
    sum_metric = temp_sum_metric;
  } // all transitions, l

  // minimum metric is the zeroth state due to the tail
  int min_metric_state = 0;
  // trace back to remove tail of zeros
  for (int l = block_length - 1; l > block_length - 1 - m; l--) {
    // previous state calculation
    min_metric_state = previous_state(min_metric_state,
                                      path_memory(min_metric_state, l));
  }

  // trace back to calculate output sequence
  for (int l = block_length - 1 - m; l >= 0; l--) {
    output(l) = get_input(min_metric_state);
    // previous state calculation
    min_metric_state = previous_state(min_metric_state,
                                      path_memory(min_metric_state, l));
  }
}


/*
  Decode a block of encoded data where encode_tailbite has been used.
*/
void Convolutional_Code::decode_tailbite(const vec &received_signal,
    bvec &output)
{
  int block_length = received_signal.size() / n; // no input symbols
  it_error_if(block_length <= 0,
              "Convolutional_Code::decode_tailbite(): Input sequence to short");
  int S0, S1;
  vec temp_sum_metric(no_states), temp_rec(n), delta_metrics;
  Array<bool> temp_visited_state(no_states);
  double temp_metric_zero, temp_metric_one;
  double best_metric = std::numeric_limits<double>::max();
  bvec best_output(block_length), temp_output(block_length);

  path_memory.set_size(no_states, block_length, false);
  output.set_size(block_length, false);


  // Try all start states ss
  for (int ss = 0; ss < no_states; ss++) {
    //Clear the visited_state vector:
    visited_state = false;
    temp_visited_state = visited_state;
    visited_state(ss) = true; // starts in the state s

    // clear accumulated metrics
    sum_metric.zeros();

    for (int l = 0; l < block_length; l++) { // all transitions
      temp_rec = received_signal.mid(l * n, n);
      // calculate all metrics for all codewords at the same time
      calc_metric(temp_rec, delta_metrics);

      for (int s = 0; s < no_states; s++) { // all states
        // S0 and S1 are the states that expanded end at state s
        previous_state(s, S0, S1);
        if (visited_state(S0)) { // expand trellis if state S0 is visited
          temp_metric_zero = sum_metric(S0)
                             + delta_metrics(output_reverse_int(s, 0));
          temp_visited_state(s) = true;
        }
        else {
          temp_metric_zero = std::numeric_limits<double>::max();
        }
        if (visited_state(S1)) { // expand trellis if state S0 is visited
          temp_metric_one = sum_metric(S1)
                            + delta_metrics(output_reverse_int(s, 1));
          temp_visited_state(s) = true;
        }
        else {
          temp_metric_one = std::numeric_limits<double>::max();
        }
        if (temp_metric_zero < temp_metric_one) { // path zero survives
          temp_sum_metric(s) = temp_metric_zero;
          path_memory(s, l) = 0;
        }
        else { // path one survives
          temp_sum_metric(s) = temp_metric_one;
          path_memory(s, l) = 1;
        }
      } // all states, s
      sum_metric = temp_sum_metric;
      visited_state = temp_visited_state;
    } // all transitions, l

    // minimum metric is the ss state due to the tailbite
    int min_metric_state = ss;

    // trace back to calculate output sequence
    for (int l = block_length - 1; l >= 0; l--) {
      temp_output(l) = get_input(min_metric_state);
      // previous state calculation
      min_metric_state = previous_state(min_metric_state,
                                        path_memory(min_metric_state, l));
    }
    if (sum_metric(ss) < best_metric) {
      best_metric = sum_metric(ss);
      best_output = temp_output;
    }
  } // all start states ss
  output = best_output;
}


/*
  Viterbi decoding using truncation of memory (default = 5*K)
*/
void Convolutional_Code::decode_trunc(const vec &received_signal,
                                      bvec &output)
{
  int block_length = received_signal.size() / n; // no input symbols
  it_error_if(block_length <= 0,
              "Convolutional_Code::decode_trunc(): Input sequence to short");
  int S0, S1;
  vec temp_sum_metric(no_states), temp_rec(n), delta_metrics;
  Array<bool> temp_visited_state(no_states);
  double temp_metric_zero, temp_metric_one;

  path_memory.set_size(no_states, trunc_length, false);
  output.set_size(0);

  // copy visited states
  temp_visited_state = visited_state;

  for (int i = 0; i < block_length; i++) {
    // update path memory pointer
    trunc_ptr = (trunc_ptr + 1) % trunc_length;

    temp_rec = received_signal.mid(i * n, n);
    // calculate all metrics for all codewords at the same time
    calc_metric(temp_rec, delta_metrics);

    for (int s = 0; s < no_states; s++) { // all states
      // the states that expanded end at state s
      previous_state(s, S0, S1);
      if (visited_state(S0)) { // expand trellis if state S0 is visited
        temp_metric_zero = sum_metric(S0)
                           + delta_metrics(output_reverse_int(s, 0));
        temp_visited_state(s) = true;
      }
      else {
        temp_metric_zero = std::numeric_limits<double>::max();
      }
      if (visited_state(S1)) { // expand trellis if state S0 is visited
        temp_metric_one = sum_metric(S1)
                          + delta_metrics(output_reverse_int(s, 1));
        temp_visited_state(s) = true;
      }
      else {
        temp_metric_one = std::numeric_limits<double>::max();
      }
      if (temp_metric_zero < temp_metric_one) { // path zero survives
        temp_sum_metric(s) = temp_metric_zero;
        path_memory(s, trunc_ptr) = 0;
      }
      else { // path one survives
        temp_sum_metric(s) = temp_metric_one;
        path_memory(s, trunc_ptr) = 1;
      }
    } // all states, s
    sum_metric = temp_sum_metric;
    visited_state = temp_visited_state;

    // find minimum metric
    int min_metric_state = min_index(sum_metric);

    // normalise accumulated metrics
    sum_metric -= sum_metric(min_metric_state);

    // check if we had enough metrics to generate output
    if (trunc_state >= trunc_length) {
      // trace back to calculate the output symbol
      for (int j = trunc_length; j > 0; j--) {
        // previous state calculation
        min_metric_state =
          previous_state(min_metric_state,
                         path_memory(min_metric_state,
                                     (trunc_ptr + j) % trunc_length));
      }
      output.ins(output.size(), get_input(min_metric_state));
    }
    else { // if not increment trunc_state counter
      trunc_state++;
    }
  } // end for (int i = 0; i < block_length; l++)
}


/*
  Calculate the inverse sequence

  Assumes that encode_tail is used in the encoding process. Returns false
  if there is an error in the coded sequence (not a valid codeword). Do
  not check that the tail forces the encoder into the zeroth state.
*/
bool Convolutional_Code::inverse_tail(const bvec coded_sequence, bvec &input)
{
  int state = 0, zero_state, one_state, zero_temp, one_temp, i, j;
  bvec zero_output(n), one_output(n);

  int block_length = coded_sequence.size() / n - m; // no input symbols
  it_error_if(block_length <= 0, "The input sequence is to short");
  input.set_length(block_length, false); // not include the tail in the output


  for (i = 0; i < block_length; i++) {
    zero_state = state;
    one_state = state | (1 << m);
    for (j = 0; j < n; j++) {
      zero_temp = zero_state & gen_pol(j);
      one_temp = one_state & gen_pol(j);
      zero_output(j) = xor_int_table(zero_temp);
      one_output(j) = xor_int_table(one_temp);
    }
    if (coded_sequence.mid(i*n, n) == zero_output) {
      input(i) = bin(0);
      state = zero_state >> 1;
    }
    else if (coded_sequence.mid(i*n, n) == one_output) {
      input(i) = bin(1);
      state = one_state >> 1;
    }
    else
      return false;
  }

  return true;
}

/*
  Check if catastrophic. Returns true if catastrophic
*/
bool Convolutional_Code::catastrophic(void)
{
  int start, S, W0, W1, S0, S1;
  bvec visited(1 << m);

  for (start = 1; start < no_states; start++) {
    visited.zeros();
    S = start;
    visited(S) = 1;

  node1:
    S0 = next_state(S, 0);
    S1 = next_state(S, 1);
    weight(S, W0, W1);
    if (S1 < start) goto node0;
    if (W1 > 0) goto node0;

    if (visited(S1) == bin(1))
      return true;
    S = S1; // goto node1
    visited(S) = 1;
    goto node1;

  node0:
    if (S0 < start) continue; // goto end;
    if (W0 > 0) continue; // goto end;

    if (visited(S0) == bin(1))
      return true;
    S = S0;
    visited(S) = 1;
    goto node1;

    //end:
    //void;
  }

  return false;
}

/*
  Calculate distance profile. If reverse = true calculate for the reverse code instead.
*/
void Convolutional_Code::distance_profile(ivec &dist_prof, int dmax, bool reverse)
{
  int max_stack_size = 50000;
  ivec S_stack(max_stack_size), W_stack(max_stack_size), t_stack(max_stack_size);

  int stack_pos = -1, t, S, W, W0, w0, w1;

  dist_prof.set_size(K, false);
  dist_prof.zeros();
  dist_prof += dmax; // just select a big number!
  if (reverse)
    W = weight_reverse(0, 1);
  else
    W = weight(0, 1);

  S = next_state(0, 1); // first state 0 and one as input, S = 1<<(m-1);
  t = 0;
  dist_prof(0) = W;

node1:
  if (reverse)
    weight_reverse(S, w0, w1);
  else
    weight(S, w0, w1);

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

      S_stack(stack_pos) = next_state(S, 0); //S>>1;
      W_stack(stack_pos) = W0;
      t_stack(stack_pos) = t + 1;
    }
  }
  else goto stack;

  W += w1;
  if (W > dist_prof(m))
    goto stack;

  t++;
  S = next_state(S, 1); //S = (S>>1)|(1<<(m-1));

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

/*
  Calculate spectrum

  Calculates both the weight spectrum (Ad) and the information weight spectrum (Cd) and
  returns it as ivec:s in the 0:th and 1:st component of spectrum, respectively. Suitable
  for calculating many terms in the spectra (uses an breadth first algorithm). It is assumed
  that the code is non-catastrophic or else it is a possibility for an eternal loop.
  dmax = an upper bound on the free distance
  no_terms = no_terms including the dmax term that should be calculated
*/
void Convolutional_Code::calculate_spectrum(Array<ivec> &spectrum, int dmax, int no_terms)
{
  imat Ad_states(no_states, dmax + no_terms), Cd_states(no_states, dmax + no_terms);
  imat Ad_temp(no_states, dmax + no_terms), Cd_temp(no_states, dmax + no_terms);
  ivec mindist(no_states),  mindist_temp(1 << m);

  spectrum.set_size(2);
  spectrum(0).set_size(dmax + no_terms, false);
  spectrum(1).set_size(dmax + no_terms, false);
  spectrum(0).zeros();
  spectrum(1).zeros();
  Ad_states.zeros();
  Cd_states.zeros();
  mindist.zeros();
  int wmax = dmax + no_terms;
  ivec visited_states(no_states), visited_states_temp(no_states);
  bool proceede;
  int d, w0, w1, s, s0, s1;

  visited_states.zeros(); // 0 = false
  s = next_state(0, 1); // Start in state from 0 with an one input.
  visited_states(s) = 1;  // 1 = true. Start in state 2^(m-1).
  w1 = weight(0, 1);
  Ad_states(s, w1) = 1;
  Cd_states(s, w1) = 1;
  mindist(s) = w1;
  proceede = true;

  while (proceede) {
    proceede = false;
    Ad_temp.zeros();
    Cd_temp.zeros();
    mindist_temp.zeros();
    visited_states_temp.zeros(); //false
    for (s = 1; s < no_states; s++) {
      if ((mindist(s) > 0) && (mindist(s) < wmax)) {
        proceede = true;
        weight(s, w0, w1);
        s0 = next_state(s, 0);
        for (d = mindist(s); d < (wmax - w0); d++) {
          Ad_temp(s0, d + w0) += Ad_states(s, d);
          Cd_temp(s0, d + w0) += Cd_states(s, d);
          visited_states_temp(s0) = 1; //true
        }

        s1 = next_state(s, 1);
        for (d = mindist(s); d < (wmax - w1); d++) {
          Ad_temp(s1, d + w1) += Ad_states(s, d);
          Cd_temp(s1, d + w1) += Cd_states(s, d) + Ad_states(s, d);
          visited_states_temp(s1) = 1; //true
        }
        if (mindist_temp(s0) > 0) { mindist_temp(s0) = (mindist(s) + w0) < mindist_temp(s0) ? mindist(s) + w0 :  mindist_temp(s0); }
        else { mindist_temp(s0) = mindist(s) + w0; }
        if (mindist_temp(s1) > 0) { mindist_temp(s1) = (mindist(s) + w1) < mindist_temp(s1) ? mindist(s) + w1 :  mindist_temp(s1); }
        else { mindist_temp(s1) = mindist(s) + w1; }

      }
    }
    Ad_states = Ad_temp;
    Cd_states = Cd_temp;
    spectrum(0) += Ad_temp.get_row(0);
    spectrum(1) += Cd_temp.get_row(0);
    visited_states = visited_states_temp;
    mindist = elem_mult(mindist_temp, visited_states);
  }
}

/*
  Cederwall's fast algorithm

  See IT No. 6, pp. 1146-1159, Nov. 1989 for details.
  Calculates both the weight spectrum (Ad) and the information weight spectrum (Cd) and
  returns it as ivec:s in the 0:th and 1:st component of spectrum, respectively. The FAST algorithm
  is good for calculating only a few terms in the spectrum. If many terms are desired, use calc_spectrum instead.
  The algorithm returns -1 if the code tested is worse that the input dfree and Cdfree.
  It returns 0 if the code MAY be catastrophic (assuming that test_catastrophic is true),
  and returns 1 if everything went right.

  \arg \c dfree the free distance of the code (or an upper bound)
  \arg \c no_terms including the dfree term that should be calculated
  \ar \c Cdfree is the best value of information weight spectrum found so far
*/
int Convolutional_Code::fast(Array<ivec> &spectrum, const int dfree, const int no_terms, const int Cdfree, const bool test_catastrophic)
{
  int cat_treshold = 7 * K; // just a big number, but not to big!
  int i;
  ivec dist_prof(K), dist_prof_rev(K);
  distance_profile(dist_prof, dfree);
  distance_profile(dist_prof_rev, dfree, true); // for the reverse code

  int dist_prof_rev0 = dist_prof_rev(0);

  bool reverse = false; // false = use dist_prof

  // is the reverse distance profile better?
  for (i = 0; i < K; i++) {
    if (dist_prof_rev(i) > dist_prof(i)) {
      reverse = true;
      dist_prof_rev0 = dist_prof(0);
      dist_prof = dist_prof_rev;
      break;
    }
    else if (dist_prof_rev(i) < dist_prof(i)) {
      break;
    }
  }

  int d = dfree + no_terms - 1;
  int max_stack_size = 50000;
  ivec S_stack(max_stack_size), W_stack(max_stack_size), b_stack(max_stack_size), c_stack(max_stack_size);
  int stack_pos = -1;

  // F1.
  int mf = 1, b = 1;
  spectrum.set_size(2);
  spectrum(0).set_size(dfree + no_terms, false); // ad
  spectrum(1).set_size(dfree + no_terms, false); // cd
  spectrum(0).zeros();
  spectrum(1).zeros();
  int S, S0, S1, W0, W1, w0, w1, c = 0;
  S = next_state(0, 1); //first state zero and one as input
  int W = d - dist_prof_rev0;


F2:
  S0 = next_state(S, 0);
  S1 = next_state(S, 1);

  if (reverse) {
    weight(S, w0, w1);
  }
  else {
    weight_reverse(S, w0, w1);
  }
  W0 = W - w0;
  W1 = W - w1;
  if (mf < m) goto F6;

  //F3:
  if (W0 >= 0) {
    spectrum(0)(d - W0)++;
    spectrum(1)(d - W0) += b;
    // The code is worse than the best found.
    if (((d - W0) < dfree) || (((d - W0) == dfree) && (spectrum(1)(d - W0) > Cdfree)))
      return -1;
  }


F4:
  if ((W1 < dist_prof(m - 1)) || (W < dist_prof(m))) goto F5;
  // select node 1
  if (test_catastrophic && W == W1) {
    c++;
    if (c > cat_treshold)
      return 0;
  }
  else {
    c = 0;
    W = W1;
  }

  S = S1;
  mf = 1;
  b++;
  goto F2;

F5:
  //if (stack_pos == -1) return n_values;
  if (stack_pos == -1) goto end;
  // get node from stack
  S = S_stack(stack_pos);
  W = W_stack(stack_pos);
  b = b_stack(stack_pos);
  c = c_stack(stack_pos);
  stack_pos--;
  mf = 1;
  goto F2;

F6:
  if (W0 < dist_prof(m - mf - 1)) goto F4;

  //F7:
  if ((W1 >= dist_prof(m - 1)) && (W >= dist_prof(m))) {
    // save node 1
    if (test_catastrophic && stack_pos > 10000)
      return 0;

    stack_pos++;
    if (stack_pos >= max_stack_size) {
      max_stack_size = 2 * max_stack_size;
      S_stack.set_size(max_stack_size, true);
      W_stack.set_size(max_stack_size, true);
      b_stack.set_size(max_stack_size, true);
      c_stack.set_size(max_stack_size, true);
    }
    S_stack(stack_pos) = S1;
    W_stack(stack_pos) = W1;
    b_stack(stack_pos) = b + 1; // because gone to node 1
    c_stack(stack_pos) = c; // because gone to node 1
  }
  // select node 0
  S = S0;
  if (test_catastrophic && W == W0) {
    c++;
    if (c > cat_treshold)
      return 0;
  }
  else {
    c = 0;
    W = W0;
  }


  mf++;
  goto F2;

end:
  return 1;
}

//----------- These functions should be moved into some other place -------

/*!
  Reverses the bitrepresentation of in (of size length) and converts to an integer
*/
int reverse_int(int length, int in)
{
  int i, j, out = 0;

  for (i = 0; i < (length >> 1); i++) {
    out = out | ((in & (1 << i)) << (length - 1 - (i << 1)));
  }
  for (j = 0; j < (length - i); j++) {
    out = out | ((in & (1 << (j + i))) >> ((j << 1) - (length & 1) + 1));
    //out = out | ( (in & (1<<j+i)) >> ((j<<1)-(length&1)+1) ); old version with preecedence problems in MSC

  }
  return out;
}

/*!
  Calculate the Hamming weight of the binary representation of in of size length
*/
int weight_int(int length, int in)
{
  int i, w = 0;
  for (i = 0; i < length; i++) {
    w += (in & (1 << i)) >> i;
  }
  return w;
}

/*!
  Compare two distance spectra. Return 1 if v1 is less, 0 if v2 less, and -1 if equal.
*/
int compare_spectra(ivec v1, ivec v2)
{
  it_assert_debug(v1.size() == v2.size(), "compare_spectra: wrong sizes");

  for (int i = 0; i < v1.size(); i++) {
    if (v1(i) < v2(i)) {
      return 1;
    }
    else if (v1(i) > v2(i)) {
      return 0;
    }
  }
  return -1;
}

/*!
  Compare two distance spectra using a weight profile.

  Return 1 if v1 is less, 0 if v2 less, and -1 if equal.
*/
int compare_spectra(ivec v1, ivec v2, vec weight_profile)
{
  double t1 = 0, t2 = 0;
  for (int i = 0; i < v1.size(); i++) {
    t1 += (double)v1(i) * weight_profile(i);
    t2 += (double)v2(i) * weight_profile(i);
  }
  if (t1 < t2) return 1;
  else if (t1 > t2) return 0;
  else return -1;
}

} // namespace itpp
