/*!
 * \file
 * \brief Definition of a binary convolutional encoder class
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

#ifndef CONVCODE_H
#define CONVCODE_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/array.h>
#include <itpp/base/binary.h>
#include <itpp/comm/channel_code.h>
#include <itpp/itexports.h>
#include <itpp/base/base_exports.h>

namespace itpp
{

/*!
  \brief Type of Convolutional Code
  \ingroup fec
*/
enum CONVOLUTIONAL_CODE_TYPE {MFD, ODS};

/*!
  \brief Encoding and decoding methods for Convolutional codes
  \ingroup fec
*/
enum CONVOLUTIONAL_CODE_METHOD {Trunc, Tail, Tailbite};


/*!
  \ingroup fec
  \brief Binary Convolutional rate 1/n class

  The codes are given as feedforward encoders and given in the
  Proakis form. That is, the binary generators (K-tuples) are
  converted to octal integers. Observe that the constraint length
  (K) is defined as the number of memory cells plus one (as in
  Proakis).

  Encoding is performed with the encode function. The default method
  for encoding is by adding a tail of K-1 zeros and also assume that
  the encoder starts in the zero state (the encode_tail()
  function). Observe that decode() by default also assumes that a
  tail is added. Both encoding and decoding method can be changed by
  the set_method() function.

  Example of use: (rate 1/3 constraint length K=7 ODS code using BPSK
  over AWGN)
  \code
  BPSK bpsk;
  Convolutional_Code code;
  ivec generator(3);
  generator(0)=0133;
  generator(1)=0165;
  generator(2)=0171;
  code.set_generator_polynomials(generator, 7);

  bvec bits=randb(100), encoded_bits, decoded_bits;
  vec tx_signal, rx_signal;

  code.encode_tail(bits, encoded_bits);
  tx_signal = bpsk.modulate_bits(encoded_bits);
  rx_signal = tx_signal + sqrt(0.5)*randn(tx_signal.size());
  code.decode_tail(rx_signal, decoded_bits);
  \endcode

  Comment: ODS-code stand for Optimum Distance Spectrum Code. For
  details see T. Ottosson, "Coding, Modulation and Multiuser
  Decoding for DS-CDMA Systems," Ph.d. thesis, Department of
  Information Theory, Scool of Electrical and Computer Engineering,
  Chalmers University of Technology, Goteborg 1997.

  It is also possible to set the generator polynomials directly using
  the builtin tables which consists of: Maximum Free Distance (MFD)
  Codes of rates R=1/2 through R=1/8 and Optimum Distance Spectrum
  (ODS) Codes of rates R=1/2 through R=1/4.
*/
class ITPP_EXPORT Convolutional_Code : public Channel_Code
{
public:
  //! Default constructor - sets (0133,0171) code with tail
  Convolutional_Code(void): K(0), start_state(0), cc_method(Tail) {
    set_code(MFD, 2, 7);
    init_encoder();
  }

  //! Destructor
  virtual ~Convolutional_Code(void) {}

  //! Set encoding and decoding method (Trunc, Tail, or Tailbite)
  void set_method(const CONVOLUTIONAL_CODE_METHOD method) {
    cc_method = method;
  }

  /*!
    \brief Set the code according to built-in tables

    The \a type_of_code can be either \a MFD or \a ODS for maximum
    free distance codes (according to Proakis) or Optimum Distance
    Spectrum Codes according to Frenger, Orten and Ottosson.
  */
  void set_code(const CONVOLUTIONAL_CODE_TYPE type_of_code, int inverse_rate,
                int constraint_length);

  //! Set generator polynomials. Given in Proakis integer form
  void set_generator_polynomials(const ivec &gen, int constraint_length);
  //! Get generator polynomials
  ivec get_generator_polynomials(void) const { return gen_pol; }

  //! Reset encoder and decoder states
  void reset();


  //@{
  //! Encode an input binary vector using specified method (Tail by default)
  virtual void encode(const bvec &input, bvec &output);
  virtual bvec encode(const bvec &input) {
    bvec output;
    encode(input, output);
    return output;
  }
  //@}

  //@{
  /*!
   * \brief Encode a binary vector starting from the previous encoder state
   *
   * The initial encoder state can be changed using set_start_state()
   * and init_encoder() functions.
   */
  void encode_trunc(const bvec &input, bvec &output);
  bvec encode_trunc(const bvec &input) {
    bvec output;
    encode_trunc(input, output);
    return output;
  }
  //@}

  //@{
  /*!
   * \brief Encoding that starts and ends in the zero state
   *
   * Encode a binary vector of inputs starting from zero state and also
   * adds a tail of K-1 zeros to force the encoder into the zero state.
   * Well suited for packet transmission.
   *
   * \note The init_encoder() function has no effect on the starting state
   * for this method.
   */
  void encode_tail(const bvec &input, bvec &output);
  bvec encode_tail(const bvec &input) {
    bvec output;
    encode_tail(input, output);
    return output;
  }
  //@}

  //@{
  /*!
   * \brief Encode an input binary vector using tailbiting
   *
   * In the Tailbiting method the starting state of the encoder is
   * initialised with the last K-1 bits of the input vector. This gives an
   * additional information to the decoder that the starting and ending
   * states are identical, although not known a priori.
   *
   * Well suited for packet transmission with small packets, because there
   * is no tail overhead as in the Tail method.
   *
   * \note The init_encoder() function has no effect on the starting state
   * for this method.
   */
  void encode_tailbite(const bvec &input, bvec &output);
  bvec encode_tailbite(const bvec &input) {
    bvec output;
    encode_tailbite(input, output);
    return output;
  }
  //@}

  //@{
  /*!
    \brief Encode a binary bit starting from the internal encoder state.

    To initialize the encoder state use set_start_state() and init_encoder()
  */
  void encode_bit(const bin &input, bvec &output);
  bvec encode_bit(const bin &input) {
    bvec output;
    encode_bit(input, output);
    return output;
  }
  //@}

  // ------------ Hard-decision decoding is not implemented ----------------
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
  virtual bvec decode(const bvec &coded_bits);

  //@{
  //! Decode a block of encoded data using specified method (Tail by default)
  virtual void decode(const vec &received_signal, bvec &output);
  virtual bvec decode(const vec &received_signal) {
    bvec output;
    decode(received_signal, output);
    return output;
  }
  //@}

  //@{
  /*!
    \brief Decode a block of encoded data where encode_tail has been used.

    Thus is assumes a decoder start state of zero and that a tail of
    K-1 zeros has been added. No memory truncation.
  */
  virtual void decode_tail(const vec &received_signal, bvec &output);
  virtual bvec decode_tail(const vec &received_signal) {
    bvec output;
    decode_tail(received_signal, output);
    return output;
  }
  //@}

  //@{
  /*!
   * \brief Decode a block of encoded data where encode_tailbite has been
   * used.
   *
   * The decoding algorithm tries all start states, so the
   * decode_tailbite() is \f$2^{K-1}\f$ times more complex than the
   * decode_tail method.
   */
  virtual void decode_tailbite(const vec &received_signal, bvec &output);
  virtual bvec decode_tailbite(const vec &received_signal) {
    bvec output;
    decode_tailbite(received_signal, output);
    return output;
  }
  //@}

  //@{
  //! Viterbi decoding using truncation of memory (default = 5*K)
  virtual void decode_trunc(const vec &received_signal, bvec &output);
  virtual bvec decode_trunc(const vec &received_signal) {
    bvec output;
    decode_trunc(received_signal, output);
    return output;
  }
  //@}


  //! Return rate of code (not including the rate-loss)
  virtual double get_rate(void) const { return rate; }


  //! Set encoder default start state.
  void set_start_state(int state) {
    it_error_if((state < 0) || ((state >= (1 << m)) && m != 0),
                "Convolutional_Code::set_start_state(): Invalid start state");
    start_state = state;
  }

  /*!
   * \brief Initialise internal encoder state with start state. Has no
   * effect on \c Tail and \c Tailbite methods.
   */
  void init_encoder() { encoder_state = start_state; }

  //!  Get the current encoder state
  int get_encoder_state(void) const { return encoder_state; }


  //! Set memory truncation length. Must be at least K.
  void set_truncation_length(const int length) {
    it_error_if(length < K, "Convolutional_Code::set_truncation_length(): "
                "Truncation length shorter than K");
    trunc_length = length;
  }

  //! Get memory truncation length
  int get_truncation_length(void) const { return trunc_length; }


  //! Check if catastrophic. Returns true if catastrophic
  bool catastrophic(void);


  /*!
    \brief Calculate the inverse sequence

    Assumes that encode_tail is used in the encoding
    process. Returns false if there is an error in the coded
    sequence (not a valid codeword). Do not check that the tail
    forces the encoder into the zeroth state.
  */
  bool inverse_tail(const bvec coded_sequence, bvec &input);


  //! \brief Calculate distance profile. If reverse = true calculate for
  //!        the reverse code instead.
  void distance_profile(ivec &dist_prof, int dmax = 100000,
                        bool reverse = false);

  /*!
    \brief Calculate spectrum

    Calculates both the weight spectrum (Ad) and the information
    weight spectrum (Cd) and returns it as ivec:s in the 0:th and
    1:st component of spectrum, respectively. Suitable for
    calculating many terms in the spectra (uses an breadth first
    algorithm). It is assumed that the code is non-catastrophic or
    else it is a possibility for an eternal loop.
    dmax = an upper bound on the free distance
    no_terms = no_terms including the dmax term that should be calculated

    Observe that there is a risk that some of the integers are
    overflow if many terms are calculated in the spectrum.
  */
  void calculate_spectrum(Array<ivec> &spectrum, int dmax, int no_terms);

  /*!
    \brief Cederwall's fast algorithm

    Calculates both the weight spectrum (Ad) and the information
    weight spectrum (Cd) and returns it as ivec:s in the 0:th and
    1:st component of spectrum, respectively. The FAST algorithm is
    good for calculating only a few terms in the spectrum. If many
    terms are desired, use calc_spectrum instead.  The algorithm
    returns -1 if the code tested is worse that the input dfree and
    Cdfree.  It returns 0 if the code MAY be catastrophic (assuming
    that test_catastrophic is true), and returns 1 if everything
    went right.

    \arg \c dfree the free distance of the code (or an upper bound)
    \arg \c no_terms including the dfree term that should be calculated
    \arg \c Cdfree is the best value of information weight spectrum found so far

    Observe that there is a risk that some of the integers are
    overflow if many terms are calculated in the spectrum.

    See IT No. 6, pp. 1146-1159, Nov. 1989 for details.
  */
  int fast(Array<ivec> &spectrum, const int dfree, const int no_terms,
           const int Cdfree = 1000000, const bool test_catastrophic = false);

protected:
  //! Next state from instate given the input
  int next_state(const int instate, const int input) {
    return ((instate >> 1) | (input << (m - 1)));
  }
  //! The previous state from state given the input
  int previous_state(const int state, const int input) {
    return (((state << 1) | input) & ((1 << m) - 1));
  }
  //! The previous state from state given the input
  void previous_state(const int state, int &S0, int &S1) {
    S0 = (state << 1) & (no_states - 1);
    S1 = S0 | 1;
  }
  //! The weight of the transition from given state with the input given
  int weight(const int state, const int input);
  //! The weight of the two paths (input 0 or 1) from given state
  void weight(const int state, int &w0, int &w1);
  //! \brief The weight (of the reverse code) of the transition from given
  //!        state with the input given
  int weight_reverse(const int state, const int input);
  //! \brief The weight (of the reverse code) of the two paths (input 0
  //!        or 1) from given state
  void weight_reverse(const int state, int &w0, int &w1);
  //! Output on transition (backwards) with input from state
  bvec output_reverse(const int state, const int input);
  //! Output on transition (backwards) with input from state
  void output_reverse(const int state, bvec &zero_output, bvec &one_output);
  //! Output on transition (backwards) with input from state
  void output_reverse(const int state, int &zero_output, int &one_output);
  //! Calculate delta metrics for 0 and 1 input branches reaching state
  void calc_metric_reverse(const int state, const vec &rx_codeword,
                           double &zero_metric, double &one_metric);
  //! Calculate delta metrics for all possible codewords
  void calc_metric(const vec &rx_codeword, vec &delta_metrics);
  //! Returns the input that results in state, that is the MSB of state
  int get_input(const int state) { return (state >> (m - 1)); }

  //! Number of generators
  int n;
  //! Constraint length
  int K;
  //! Memory of the encoder
  int m;
  //! Number of states
  int no_states;
  //! Generator polynomials
  ivec gen_pol;
  //! Generator polynomials for the reverse code
  ivec gen_pol_rev;
  //! The current encoder state
  int encoder_state;
  //! The encoder start state
  int start_state;
  //! The decoder truncation length
  int trunc_length;
  //! The rate of the code
  double rate;
  //! Auxilary table used by the codec
  bvec xor_int_table;
  //! output in int format for a given state and input
  imat output_reverse_int;
  //! encoding and decoding method
  CONVOLUTIONAL_CODE_METHOD cc_method;
  //! Path memory (trellis)
  imat path_memory;
  //! Visited states
  Array<bool> visited_state;
  //! Metrics accumulator
  vec sum_metric;
  //! Truncated path memory pointer
  int trunc_ptr;
  //! Truncated memory fill state
  int trunc_state;
};

// --------------- Some other functions that maybe should be moved -----------
/*!
  \relates Convolutional_Code
  \brief Reverses the bitrepresentation of in (of size length) and converts to an integer
*/
ITPP_EXPORT int reverse_int(int length, int in);

/*!
  \relates Convolutional_Code
  \brief Calculate the Hamming weight of the binary representation of in of size length
*/
ITPP_EXPORT int weight_int(int length, int in);

/*!
  \relates Convolutional_Code
  \brief Compare two distance spectra. Return 1 if v1 is less, 0 if v2 less, and -1 if equal.
*/
ITPP_EXPORT int compare_spectra(ivec v1, ivec v2);

/*!
  \relates Convolutional_Code
  \brief Compare two distance spectra using a weight profile.

  Return 1 if v1 is less, 0 if v2 less, and -1 if equal.
*/
ITPP_EXPORT int compare_spectra(ivec v1, ivec v2, vec weight_profile);

} // namespace itpp

#endif // #ifndef CONVCODE_H
