/*!
 * \file
 * \brief Definitions of a Recursive Systematic Convolutional codec class
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

#ifndef REC_SYST_CONV_CODE_H
#define REC_SYST_CONV_CODE_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/comm/convcode.h>
#include <itpp/comm/llr.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief A Recursive Systematic Convolutional Encoder/Decoder class

  The main purpose of this class is to be used by the Turbo_Codec class which uses two recursive systematic
  convolutional encoders. It can however be used as a stand alone class. The map_decode
  member function implementation follows the paper "A Turbo Code Tutorial" by William
  E. Ryan, New Mexico State University. This paper was found on the web and is probably
  unpublished.
*/
class ITPP_EXPORT Rec_Syst_Conv_Code
{
public:

  //! Class constructor
  Rec_Syst_Conv_Code(): infinity(1e30) {}

  //! Class constructor
  virtual ~Rec_Syst_Conv_Code() {}

  /*!
    \brief Set generator polynomials.

    The generator polynomials are given in Proakis integer form. First generator (gen(0)) is the recursive polynomial.

    \param gen A vector containing the generator polynomials of the RSC Code.
    \param constraint_length The Constraint length of the encoder.
  */
  void set_generator_polynomials(const ivec &gen, int constraint_length);

  /*!
    \brief Sets the channel parameters needed for MAP-decoding.

    \param Ec is the energy per channel symbol
    \param N0 is the single-sided power spectral density of the AWGN on the channel.
  */
  void set_awgn_channel_parameters(double Ec, double N0);

  /*!
    \brief Set scaling factor for the decoder

    \param in_Lc is the channel reliability factor (i.e. \a Lc = 4 x sqrt( \a Ec ) / \a N0)
  */
  void set_scaling_factor(double in_Lc);

  /*!
    \brief Set the lookup table for algebra with quantized LLR values (see \c LLR_calc_unit class)
  */
  void set_llrcalc(LLR_calc_unit in_llrcalc);

  /*!
    \brief Encode a binary vector of inputs and also adds a tail of \a K-1 zeros to force the encoder into the zero state.

    The encoder remembers that the trellis is terminated in the zero state at the end of the input block. This is then
    utilized by the decoder when going through the trellis in the reverse direction. The tailbits used are returned in \a tail.
    Parity bits for both the input part and the tail part of the data are returned in the matrix \a parity bits.
  */
  void encode_tail(const bvec &input, bvec &tail, bmat &parity_bits);

  /*!
    \brief Encode a binary vector of inputs starting from zero state without adding of a tail
  */
  void encode(const bvec &input, bmat &parity_bits);

  /*!
    \brief Maximum A Posteriori (MAP) Probability symbol-by-symbol decoder.

    The extrinsic_input is the a priori information on each systematic bit. If no a priori information is available, this vector should
    contain only zeros. The extrinsic_output term may be passed to a subsequent decoder in a Turbo
    scheme. The decision variable is \code L = Lc*rec_systematic + extrinsic_output + extrinsic_input \endcode where \code Lc = 4*sqrt(Ec)/N0 \endcode

    \param rec_systematic Including both systematic bits and tail bits (if any)
    \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
    \param extrinsic_input For all systematic bits
    \param extrinsic_output For all systematic bits
    \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise

    <b>Note:</b> It is recommended to use the log_decode() decoder instead, as it is much faster and more numerically stable.
  */
  virtual void map_decode(const vec &rec_systematic, const mat &rec_parity, const vec &extrinsic_input, vec &extrinsic_output,
                          bool set_terminated = false);

  /*!
    \brief Log domain implementation of the Maximum A Posteriori (MAP) Probability symbol-by-symbol decoder.

    The extrinsic_input is the a priori information on each systematic bit. If no a priori information is available, this vector should
    contain only zeros. The extrinsic_output term may be passed to a subsequent decoder in a Turbo
    scheme. The decision variable is \code L = Lc*rec_systematic + extrinsic_output + extrinsic_input \endcode where \code Lc = 4*sqrt(Ec)/N0 \endcode

    \param rec_systematic Including both systematic bits and tail bits (if any)
    \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
    \param extrinsic_input For all systematic bits
    \param extrinsic_output For all systematic bits
    \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise
    \param metric May be "LOGMAP", "LOGMAX" (default), or "TABLE"

    <b>Note:</b> Unless LOGMAX decoding is desired, it is
    recommended to use the TABLE metric instead of LOGMAP as the
    table-based decoder is much faster and numerically stable.
  */
  virtual void log_decode(const vec &rec_systematic, const mat &rec_parity, const vec &extrinsic_input,
                          vec &extrinsic_output, bool set_terminated = false, std::string metric = "LOGMAX");

  /*!
    \brief Special Log-MAP/Log-MAX decoder implementation for \a n = 2

    \param rec_systematic Including both systematic bits and tail bits (if any)
    \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
    \param extrinsic_input For all systematic bits
    \param extrinsic_output For all systematic bits
    \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise
    \param metric May be "LOGMAP", "LOGMAX" (default), or "TABLE"

    <b>Note:</b> Unless LOGMAX decoding is desired, it is
    recommended to use the TABLE metric instead of LOGMAP as the
    table-based decoder is much faster and numerically stable.
  */
  virtual void log_decode_n2(const vec &rec_systematic,
                             const vec &rec_parity,
                             const vec &extrinsic_input,
                             vec &extrinsic_output,
                             bool set_terminated = false,
                             std::string metric = "LOGMAX");

  // ===== EGL: ADDED FUNCTIONS NOV 2005 (THESE ARE DERIVATIVES OF EXISTING FUNCTIONS) ======

  /*!  \brief Implementation of the log-map decoder using quantized
    LLR values (the \c QLLR type) and table-lookup (using the \c
    LLR_calc_unit class).

    \param rec_systematic Including both systematic bits and tail bits (if any)
    \param rec_parity Matrix including all parity bits from all polynomials as well
    as parity bits from the tail (if terminated)
    \param extrinsic_input For all systematic bits
    \param extrinsic_output For all systematic bits
    \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise

  */
  virtual void log_decode(const QLLRvec &rec_systematic,
                          const QLLRmat &rec_parity,
                          const QLLRvec &extrinsic_input,
                          QLLRvec &extrinsic_output,
                          bool set_terminated = false);

  /*!  \brief Implementation of the log-map decoder for the n=2 case
    using quantized LLR values (the \c QLLR type) and table-lookup
    (using the \c LLR_calc_unit class).

    \param rec_systematic Including both systematic bits and tail bits (if any)
    \param rec_parity Matrix including all parity bits from all polynomials as well
    as parity bits from the tail (if terminated)
    \param extrinsic_input For all systematic bits
    \param extrinsic_output For all systematic bits
    \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise

   */
  virtual void log_decode_n2(const QLLRvec &rec_systematic,
                             const QLLRvec &rec_parity,
                             const QLLRvec &extrinsic_input,
                             QLLRvec &extrinsic_output,
                             bool set_terminated = false);

  // ========================================================

  //! Dummy assignment operator - MSVC++ warning C4512
  Rec_Syst_Conv_Code & operator=(const Rec_Syst_Conv_Code &) { return *this; }

private:
  //! Used for precalculations of the trellis state transitions
  int calc_state_transition(const int instate, const int input, ivec &parity);

  int n, K, m;
  ivec gen_pol, gen_pol_rev;
  int encoder_state, Nstates;
  double rate, Lc;
  imat state_trans, output_parity, rev_state_trans, rev_output_parity;
  bool terminated;
  double ln2;

  /*!
    This instance of an \c LLR_calc_unit contains the tables used for table lookup
    in the table-based map decoder.
  */
  LLR_calc_unit llrcalc;

  // This const value replaces INT definition used previously
  const double infinity;
};

} // namespace itpp

#endif // #ifndef REC_SYST_CONV_CODE_H
