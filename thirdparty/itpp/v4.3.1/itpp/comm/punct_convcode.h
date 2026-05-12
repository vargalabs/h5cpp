/*!
 * \file
 * \brief Definitions of a Binary Punctured Convolutional Encoder class
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

#ifndef PUNCT_CONVCODE_H
#define PUNCT_CONVCODE_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/array.h>
#include <itpp/comm/convcode.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief Binary Punctured Convolutional Code Class.

  The codes are given as feedforward encoders an given in the Proakis form. That is the binary
  generators (K-tuples) are converted to octal integers. Observe that the constraint length (K)
  is defined as the number of meomory cells plus one (as in Proakis). The puncture matrix should
  be of size \a n * \a Period, where \a Period is the puncturing period.

  Encoding is performed with the \c encode function. By default the \c encode_tail function
  is called which automatically add a tail of K-1 zeros and also assume that the encoder starts in
  the zero state. Observe that \c decode_tail is used for data
  encoded with \c encode_tail, and decode_trunc assumes that the memory truncation length either is the
  default (5*K) or set using the \c set_truncation_length function. Encoding and decoding method can
  be changed by calling the set_method() function.

  Example of use: (rate 1/3 constraint length K=7 ODS code using BPSK over AWGN)
  \code
  BPSK bpsk;
  Punctured_Convolutional_Code code;
  ivec generator(3);
  generator(0)=0133;
  generator(1)=0165;
  generator(2)=0171;
  code.set_generator_polynomials(generator, 7);

  bmat puncture_matrix = "1 1;0 1";
  code.set_puncture_matrix(puncture_matrix);
  code.set_truncation_length(30);

  bvec bits=randb(100), encoded_bits, decoded_bits;
  vec tx_signal, rx_signal;

  code.encode(bits, encoded_bits);
  tx_signal = bpsk.modulate_bits(encoded_bits);
  rx_signal = tx_signal + sqrt(0.5)*randn(tx_signal.size());
  code.decode(rx_signal, decoded_bits);
  \endcode
*/
class ITPP_EXPORT Punctured_Convolutional_Code : public Convolutional_Code
{
public:
  //! Constructor
  Punctured_Convolutional_Code(void) : Convolutional_Code() {}
  //! Destructor
  virtual ~Punctured_Convolutional_Code(void) {}

  /*!
    \brief Set the code according to built-in tables

    The \a type_of_code can be either \a MFD or \a ODS for maximum free distance codes (according to Proakis)
    or Optimum Distance Spectrum Codes according to Frenger, Orten and Ottosson.
  */
  void set_code(const CONVOLUTIONAL_CODE_TYPE type_of_code, int inverse_rate, int constraint_length)
  { Convolutional_Code::set_code(type_of_code, inverse_rate, constraint_length); }
  //! Set generator polynomials. Given in Proakis integer form
  void set_generator_polynomials(const ivec &gen, int constraint_length)
  { Convolutional_Code::set_generator_polynomials(gen, constraint_length); }
  //! Get generator polynomials
  ivec get_generator_polynomials() const { return gen_pol; }

  //! Return rate of code
  virtual double get_rate() const { return rate; }

  //! Set encoding and decoding method (Trunc, Tail, or Tailbite)
  void set_method(const CONVOLUTIONAL_CODE_METHOD method) { Convolutional_Code::set_method(method); }

  //! Set puncture matrix (size n*Period)
  void set_puncture_matrix(const bmat &pmatrix); // add test of matrix size
  //! Get puncture matrix
  bmat get_puncture_matrix() const { return puncture_matrix; }
  //! Get puncturing period
  int get_puncture_period() const { return Period; }

  //! Set the encoder internal state in start_state (set by set_start_state()).
  void init_encoder() { encoder_state = start_state; }

  //! Encode a binary vector of inputs using specified method
  void encode(const bvec &input, bvec &output);
  //! Encode a binary vector of inputs using specified method
  bvec encode(const bvec &input) { bvec output; encode(input, output); return output; }

  //! Encode a binary vector of inputs starting from state set by the set_state function.
  void encode_trunc(const bvec &input, bvec &output);
  //! Encode a binary vector of inputs starting from state set by the set_state function.
  bvec encode_trunc(const bvec &input) { bvec output; encode_trunc(input, output); return output; }

  /*!
    \brief Encoding that begins and ends in the zero state

    Encode a binary vector of inputs starting from zero state and also adds a tail
    of K-1 zeros to force the encoder into the zero state. Well suited for packet
    transmission.
  */
  void encode_tail(const bvec &input, bvec &output);
  /*!
    \brief Encoding that begins and ends in the zero state

    Encode a binary vector of inputs starting from zero state and also adds a tail
    of K-1 zeros to force the encoder into the zero state. Well suited for packet
    transmission.
  */
  bvec encode_tail(const bvec &input) { bvec output; encode_tail(input, output); return output; }

  //! Encode a binary vector of inputs using tailbiting.
  void encode_tailbite(const bvec &input, bvec &output);
  //! Encode a binary vector of inputs using tailbiting.
  bvec encode_tailbite(const bvec &input)
  { bvec output; encode_tailbite(input, output); return output; }


  //! Viterbi decoding using specified method
  virtual void decode(const vec &received_signal, bvec &output);
  //! Viterbi decoding using specified method
  virtual bvec decode(const vec &received_signal) { bvec output; decode(received_signal, output); return output; }

  // ------------ Hard-decision decoding is not implemented -------------------
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
  virtual bvec decode(const bvec &coded_bits);

  //! Viterbi decoding using truncation of memory (default = 5*K)
  void decode_trunc(const vec &received_signal, bvec &output);
  //! Viterbi decoding using truncation of memory (default = 5*K)
  bvec decode_trunc(const vec &received_signal) { bvec output; decode_trunc(received_signal, output); return output; }

  /*!
    \brief Decode a block of encoded data where encode_tail has been used.

    Thus is assumes a decoder start state of zero and that a tail of
    K-1 zeros has been added. No memory truncation.
  */
  void decode_tail(const vec &received_signal, bvec &output);
  /*!
    \brief Decode a block of encoded data where encode_tail has been used.

    Thus is assumes a decoder start state of zero and that a tail of
    K-1 zeros has been added. No memory truncation.
  */
  bvec decode_tail(const vec &received_signal) { bvec output; decode_tail(received_signal, output); return output; }

  //! Decode a block of encoded data where encode_tailbite has been used. Tries all start states.
  void decode_tailbite(const vec &received_signal, bvec &output);
  //! Decode a block of encoded data where encode_tailbite has been used. Tries all start states.
  bvec decode_tailbite(const vec &received_signal)
  { bvec output; decode_tailbite(received_signal, output); return output; }

  /*
    \brief Calculate the inverse sequence

    Assumes that encode_tail is used in the encoding process. Returns false if there is an error in the coded sequence
    (not a valid codeword).
  */
  bool inverse_tail(const bvec coded_sequence, bvec &input);

  //! Check if the code is catastrophic. Returns true if catastrophic
  bool catastrophic(void);

  //! Calculate distance profile. If reverse = true calculate for the reverse code instead.
  void distance_profile(ivec &dist_prof, int time, int dmax = 100000, bool reverse = false);

  /*!
    \brief Calculate spectrum.

    Calculates both the weight spectrum (Ad) and the information weight spectrum (Cd) and
    returns it as ivec:s in the 0:th and 1:st component of spectrum, respectively. For a
    punctrued code the spectrum is a sum of the spectras of all starting positions.
    Suitable for calculating many terms in the spectra (uses an breadth first algorithm).
    It is assumed that the code is non-catastrophic or else it is a possibility for an eternal loop.

    <ul>
    <li> \a dmax = an upper bound on the free distance </li>
    <li> \a no_terms = number of terms including the \a dmax term that should be calculated </li>
    </ul>

    Observe that there is a risk that some of the integers are overflow if many terms are calculated in the spectrum.
  */
  void calculate_spectrum(Array<ivec> &spectrum, int dmax, int no_terms);

  /*!
    \brief Calculate spectrum. Suitable when calculating many terms in the spectra. Breadth first search.

    Use this function to evaluate the spectum whith a speccific puncturing period, or to calculate the
    spectrum for block transmission. To calculate spectra for block transmission:
    <ul>
    <li> Use \a time = 0 if the puncturing is restarted at each block. </li>
    <li> Use \a block_length = 0 (default value) for infinite blocks. </li>
    </ul>
  */
  void calculate_spectrum(Array<ivec> &spectrum, int time, int dmax, int no_terms, int block_length = 0);

  /*!
    \brief Cederwall's fast algorithm.

    <ul>
    <li> See IEEE Trans. Information Theory No. 6, pp. 1146-1159, Nov. 1989 for details.  </li>
    <li> Calculates both the weight spectrum (Ad) and the information weight spectrum (Cd) and returns it as ivec:s in the 0:th and 1:st component of spectrum, respectively. </li>
    <li> The algorithm returns -1 if the code tested is worse that the input \a dfree and \a Cdfree. </li>
    <li> It returns 0 if the code MAY be catastrophic (assuming that \a test_catastrophic is \c true), and returns 1 if everything went right. </li>
    <li> \a dfree = the free distance of the code (or an upper bound). </li>
    <li> \a no_terms = Number of terms including the \a dfree term that should be calculated. </li>
    <li> \a d_best_so_far = the best value of the free distance found so far. </li>
    <li> The FAST algorithm is good for calculating only a few terms in the spectrum. If many terms are desired, use \c calc_spectrum instead. </li>
    </ul>

    Observe that there is a risk that some of the integers are overflow if many terms are calculated in the spectrum.
  */
  int fast(Array<ivec> &spectrum, int time, int dfree, int no_terms, int d_best_so_far = 0, bool test_catastrophic = false);

protected:
  //! The weight of path from \c state with \c input (0 or 1) at transition \c time
  int weight(const int state, const int input, int time);
  //! The weight of the two paths (input 0 or 1) from given state
  void weight(const int state, int &w0, int &w1, int time);
  //! Weight of the reverse code from \c state with \c input (0 or 1) at transition \c time
  int weight_reverse(const int state, const int input, int time);
  //! The weight of the reverse code of two paths (input 0 or 1) from given state
  void weight_reverse(const int state, int &w0, int &w1, int time);

  //! The puncture period (i.e. the number of columns in the puncture matrix)
  int Period;
  //! The number of "1" in the puncture matrix
  int total;
  //! The puncture matrix (\a n rows and \a Period columns)
  bmat puncture_matrix;
};

} // namespace itpp

#endif // #ifndef PUNCT_CONVCODE_H
