/*!
 * \file
 * \brief Definition of a (punctured) turbo encoder/decoder class
 * \author Pal Frenger and Zbigniew Dlugaszewski. QLLR support by Erik G. Larsson.
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

#ifndef TURBO_H
#define TURBO_H

#include <itpp/comm/rec_syst_conv_code.h>
#include <itpp/comm/interleave.h>
#include <itpp/comm/llr.h>
#include <itpp/itexports.h>

namespace itpp
{

//! \cond

#if (defined(_MSC_VER) && defined (ITPP_SHARED_LIB))
//MSVC explicitely instantiate required template while building the shared library
template class ITPP_EXPORT Sequence_Interleaver<bin>;
template class ITPP_EXPORT Sequence_Interleaver<double>;
#endif

//! \endcond


/*!
  \brief Turbo encoder/decoder Class
  \ingroup fec
  \author Pal Frenger

  To set up the turbo encoder used in e.g. WCDMA the following code can be used (assuming a code block size of 320 bits):
  \code
  Turbo_Codec turbo;
  ivec gen(2);
  gen(0) = 013; gen(1) = 015;
  int constraint_length = 4;
  ivec interleaver_sequence = wcdma_turbo_interleaver_sequence( 320 );
  turbo.set_parameters(gen, gen, constraint_length, interleaver_sequence);
  \endcode
*/
class ITPP_EXPORT Turbo_Codec
{
public:

  //! Class constructor
  Turbo_Codec(void) {}

  //! Class destructor
  virtual ~Turbo_Codec(void) {}

  /*!
    \brief Set parameters for the turbo encoder/decoder

    \param gen1 A vector with \a n1 elements containing the generator polynomials for the first constituent encoder
    \param gen2 A vector with \a n2 elements containing the generator polynomials for the second constituent encoder
    \param constraint_length The constraint length of the two constituent encoders
    \param interleaver_sequence An ivec defining the internal turbo interleaver.
    \param in_iterations The number of decoding iterations. Default value is 8.
    \param in_metric Determines the decoder metric: "MAP", LOGMAP", "LOGMAX", or "TABLE". The default is "LOGMAX".
    \param in_logmax_scale_factor The extrinsic information from each constituent decoder is to optimistic when
    LOGMAX decoding is used.
    This parameter allows for a down-scaling of the extrinsic information that will be passed on to the next decoder.
    The default value
    is 1.0. This parameter is ignored for other metrics than "LOGMAX".
    \param in_adaptive_stop If this parameter is true, then the iterations will stop if the decoding results after one
    full iteration equals the previous iteration. Default value is false.

    \param lcalc This parameter can be used to provide a specific \c LLR_calc_unit which defines the resolution in
    the table-lookup if decoding with the metric "TABLE" is used.

    \note For issues relating to the accuracy of LLR computations,
    please see the documentation of \c LLR_calc_unit
  */
  void set_parameters(ivec gen1, ivec gen2, int constraint_length,
                      const ivec &interleaver_sequence, int in_iterations = 8,
                      const std::string &in_metric = "LOGMAX", double in_logmax_scale_factor = 1.0,
                      bool in_adaptive_stop = false, LLR_calc_unit lcalc = LLR_calc_unit());

  /*!
    \brief Set a new internal interleaver sequence for the turbo encoder/decoder

    By changing the interleaver sequence it is possible to change the code word size of the turbo codec. Note that you
    still must use the \a set_parameters function first to set the other parameters of the turbo codec.
  */
  void set_interleaver(const ivec &interleaver_sequence);

  /*!
    \brief Set the decoder metric

    \param in_metric Determines the decoder metric: "MAP", LOGMAP", "LOGMAX", or "TABLE". The default is "LOGMAX".
    \param in_logmax_scale_factor The extrinsic information from each constituent decoder is to optimistic when
    LOGMAX decoding is used.
    This parameter allows for a down-scaling of the extrinsic information that will be passed on to the next decoder.
    The default value is 1.0. This parameter is ignored for other metrics than "LOGMAX".

    \param lcalc This parameter can be used to provide a specific \c LLR_calc_unit which defines the resolution in
    the table-lookup if decoding with the metric "TABLE" is used.
  */
  void set_metric(std::string in_metric = "LOGMAX", double in_logmax_scale_factor = 1.0,
                  LLR_calc_unit lcalc = LLR_calc_unit());

  /*!
    \brief Sets the number of decoding iterations. Default value is 8.
  */
  void set_iterations(int in_iterations = 8);

  /*!
    \brief Use and adaptive number of iterations

    When the adaptive stop criterion is used the iterations will stop before the maximum number of iterations is
    executed if the decoding results after one full iteration equals the previous iteration. Default value is \c true.
  */
  void set_adaptive_stop(bool in_adaptive_stop = true);

  /*!
    \brief Set parameters for decoding on an AWGN channel

    \param in_Ec The received energy per channel symbol (i.e. per soft input bit)
    \param in_N0 The single sided spectral density of the AWGN noise
  */
  void set_awgn_channel_parameters(double in_Ec, double in_N0);

  /*!
    \brief Set scaling factor for decoding on e.g. Rayleigh fading channels

    Setting the correct value of the channel reliability function is important for MAP decoder algorithms. However, if
    the Log-MAX decoding algorithm is used, then the value of \a Lc is not important. Assuming that the received soft
    values \f$r_k\f$ from the channel equal

    \f[ r_k = h_k c_k + w_k \f]

    where \f$h_k\f$ is the (complex valued) channel gain, \f$c_k\f$ is the transmitted symbol equal to
    \f$\{-\sqrt{E_c},+\sqrt{E_c}\}\f$, and \f$w_k\f$ is the AWGN (complex valued) noise with total variance
    of the real plus imaginary part equal to \f$N_0\f$. The input to the turbo decoder shall then be

    \f[ z_k = \hat{h}_k^{*} r_k \f]

    where \f$\hat{h}_k^{*}\f$ is the conjugate of the channel estimate. Assuming that the channel estimate is perfect,
    the channel reliability factor shall be set to

    \f[ L_c = 4\sqrt{E_c} / {N_0} \f]

    \param in_Lc the channel reliability factor of the channel.
  */
  void set_scaling_factor(double in_Lc);

  /*!
    \brief Encoder function

    This function can encode several consecutive coding blocks. The output is organized as follows:

    \f[ s(1), p_{1,1}(1), p_{1,2}(1), \ldots , p_{1,n_1}(1), p_{2,1}(1), p_{2,2}(1), \ldots , p_{2,n_2}(1), s(2), \ldots \f]

    In the above expression \f$s(n)\f$ is the n-th systematic bit and \f$p_{l,k}(n)\f$ is the n-th bit from the k-th
    encoder polynomial of the l-th constituent encoder. A tail of both systematic and parity bits is added after each
    coding block to force both encoder to the zero state. The tail of each encoder is structured as follows (using
    encoder one as an example):

    \f[ t_1(1), pt_{1,1}(1), pt_{1,2}(1), \ldots , pt_{1,n_1}(1), \ldots pt_{1,n_1}(m) \f]

    The tailbits from the first encoder are placed before the tailbits from the second encoder.

    \param input The input bits to the encoder. Must contain an integer number of code blocks
    \param output The encoded bits including two tail, one from each constituent encoder, after each coding block.
  */
  void encode(const bvec &input, bvec &output);

  /*!
    \brief Decoder function

    This function can decode several consecutive coding blocks that were encoded with the encode member function

    \param received_signal The vector of received bits
    \param decoded_bits A vector of decoded bits
    \param true_bits If this input vector is provided then the iterations will stop as soon as the decoded bits
    equals the \c true_bits. Note that this feature can not be used if the \c in_adaptive_stop parameter in the
    setup function is set to \c true.
  */
  virtual void decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits = "0");

  /*!
    \brief Decoder function

    This function can decode several consecutive coding blocks that were encoded with the encode member function

    \param received_signal The vector of received bits
    \param decoded_bits A vector of decoded bits
    \param nrof_used_iterations Returns the number of used iterations for each code block.
    \param true_bits If this input vector is provided then the iterations will stop as soon as the decoded bits
    equals the \c true_bits. Note that this feature can not be used if the \c in_adaptive_stop parameter in the
    setup function is set to \c true.
  */
  virtual void decode(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                      const bvec &true_bits = "0");

  /*!
    \brief Encode a single block

    This function is useful if rate matching is to be applied after the turbo encoder. The size of \a in1 and \a in2
    equals the size of \a input plus the tail bits. Tailbits are appended ate the end of \a in1 and \a in2. The number
    of rows in \a parity1 and \a parity2 equals the lengths of \a in1 and \a in2 respectively. The number of columns of
    \a parity1 and \a parity2 equals the number of parity bits from the first and the second constituent encoders
    respectively.

    \param input The input bits to the encoder. Must contain a single code block
    \param in1 The input bits to the first constituent encoder with a tail added at the end
    \param in2 The input bits to the second constituent encoder (i.e. the interleaved data bits) with a tail
    added at the end
    \param parity1 The parity bits from the first constituent encoder (including parity bits for the first tail)
    \param parity2 The parity bits from the second constituent encoder (including parity bits for the second tail)
  */
  void encode_block(const bvec &input, bvec &in1, bvec &in2, bmat &parity1, bmat &parity2);

  /*!
    \brief Decode a single block

    This function can decode a single coding blocks that was encoded with the encode_block member function.
    In order to speed up the decoding process it is possible to input the correct data bits. If this information
    is provided the decoder can stop iterating as soon as the decoded bits match the correct data bits. This
    simulation trick can greatly speed up the simulation time for high SNR cases when only a few iterations are
    required. If errors still exist after the maximum number of iterations have been performed, the decoding
    process stops.

    The matrix \a decoded_bits_i contains the result from decoding iteration \a i on row \a i. Even though both
    \a rec_syst1 and \a rec_syst2 are given as inputs, the systematic bits in \a rec_syst2 will in most cases be
    punctured and only the tailbits at the end of the vector \a rec_syst2 will have values different from zero.

    \note This decoding function assumes that the input is scaled as +-2*SNR + noise. This means that the channel
    reliability factor \a Lc must be equal to 1.0. No additional scaling is performed by this function.

    \param rec_syst1 The received input bits to the first constituent decoder with a tail added at the end
    \param rec_syst2 The received input bits to the second constituent decoder with a tail added at the end
    \param rec_parity1 The received parity bits for the first constituent decoder (including parity bits for the
    first tail)
    \param rec_parity2 The received parity bits for the second constituent decoder (including parity bits for
    the second tail)
    \param decoded_bits_i Contains the result from decoding iteration \a i on row \a i.
    \param nrof_used_iterations_i Returns the number of iterations used for decoding
    of this block.
    \param true_bits Optional input parameter. If given, the iterations will stop as soon as the decoded bits
    match the true bits.
  */
  virtual void decode_block(const vec &rec_syst1, const vec &rec_syst2, const mat &rec_parity1, const mat &rec_parity2,
                            bmat &decoded_bits_i, int &nrof_used_iterations_i, const bvec &true_bits = "0");

  //! Get number of coded bits
  int get_Ncoded() const { return Ncoded; }

  //! Get number of uncoded bits
  int get_Nuncoded() const { return Nuncoded; }

protected:

  /*!
    \brief Special decoder function for \a R = 1/3 i.e. two parity bits for each systematic bit
  */
  void decode_n3(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                 const bvec &true_bits = "0");

  //Scalars:
  int interleaver_size;
  int Ncoded, Nuncoded;
  int m_tail, n1, n2, n_tot, iterations;
  double Ec, N0, Lc, R, logmax_scale_factor;
  bool adaptive_stop;
  struct ITPP_EXPORT Metric
  {
    enum Type {Unknown, LOGMAX, LOGMAP, MAP, TABLE};
    Metric() : _t(Unknown) {}
    Metric(Type t) : _t(t) {}
    operator Type () const {return _t;}
  private:
    Type _t;
    template<typename T> operator T () const;
  };
  Metric metric;

  //Vectors:
  bvec decoded_bits_previous_iteration;

  //Classes:
  Rec_Syst_Conv_Code rscc1, rscc2;
  Sequence_Interleaver<bin> bit_interleaver;
  Sequence_Interleaver<double> float_interleaver;
  static std::string string_from_metric(const Metric& m);
};

/*!
  \brief Punctured turbo encoder/decoder Class
  \ingroup fec
  \author Zbigniew Dlugaszewski

  To set up the turbo encoder rate 1/2 with random interleaver the following code can be used:
  \code
  Punctured_Turbo_Codec turbo;
  ivec gen(2);
  gen(0) = 013; gen(1) = 015;
  int constraint_length = 4;
  int block_length = 400;
  turbo.set_parameters(gen, gen, constraint_length, block_length);
  bmat puncture_matrix = "1 1;1 0;0 1";
  turbo.set_puncture_matrix(puncture_matrix);

  \endcode
*/
class ITPP_EXPORT Punctured_Turbo_Codec : public Turbo_Codec
{
public:

  //! Class constructor
  Punctured_Turbo_Codec(void) { Period = 0; };

  //! Class destructor
  virtual ~Punctured_Turbo_Codec(void) {}

  /*!
    \brief Set parameters for the punctured turbo encoder/decoder

    the same parameters as with \c Turbo_Codec :
    \param gen1 A vector with \a n1 elements containing the generator polynomials for the first constituent encoder
    \param gen2 A vector with \a n2 elements containing the generator polynomials for the second constituent encoder
    \param constraint_length The constraint length of the two constituent encoders
    \param interleaver_sequence An ivec defining the internal turbo interleaver.
    \param in_iterations The number of decoding iterations. Default value is 8.
    \param in_metric Determines the decoder metric: "MAP", LOGMAP", "LOGMAX", or "TABLE". The default is "LOGMAX".
    \param in_logmax_scale_factor The extrinsic information from each constituent decoder is to optimistic when
    LOGMAX decoding is used.
    This parameter allows for a down-scaling of the extrinsic information that will be passed on to the next decoder.
    The default value
    is 1.0. This parameter is ignored for other metrics than "LOGMAX".
    \param in_adaptive_stop If this parameter is true, then the iterations will stop if the decoding results after one
    full iteration equals the previous iteration. Default value is false.

    \param lcalc This parameter can be used to provide a specific \c LLR_calc_unit which defines the resolution in
    the table-lookup if decoding with the metric "TABLE" is used.

    additionally:
    \param pmatrix Puncturing matrix containing "1" for transmitted bits and "0" for punctured ones with \a n_tot rows and \a Period columns.

    \note For issues relating to the accuracy of LLR computations,
    please see the documentation of \c LLR_calc_unit
  */
  void set_parameters(ivec gen1, ivec gen2, int constraint_length, const ivec &interleaver_sequence, bmat &pmatrix, int in_iterations = 8, std::string in_metric = "LOGMAX", double in_logmax_scale_factor = 1.0, bool in_adaptive_stop = false, itpp::LLR_calc_unit lcalc = itpp::LLR_calc_unit());

  /*!
    \brief Sets puncturing matrix for turbo code (size \a n_tot * \a Period)

    \param pmatrix Puncturing matrix containing "1" for transmitted bits and "0" for punctured ones with \a n_tot rows and \a Period columns.
  */
  void set_puncture_matrix(const bmat &pmatrix);

  /*!
    \brief Returns puncturing matrix
  */
  bmat get_puncture_matrix(void) { return puncture_matrix; };

  /*!
    \brief Returns \a Period (number of columns) of the puncturing matrix
  */
  int get_puncture_period(void) { return Period; };

  /*!
    \brief Returns size of a single punctured block

    Returns size of a single encoded block after puncturing.
  */
  int get_punctured_size(void) { it_assert(Period != 0, "Punctured_Turbo_Codec: puncture matrix is not set"); return pNcoded; };

  /*!
    \brief Returns rate of the turbo code

    \param nominal Determines whether the function returns nominal code rate or the actual code rate (taking into account tail bits). Default value is true, i.e. nominal rate is returned.
  */
  double get_rate(bool nominal = true);

  /*!
    \brief Returns size of a single block prior to the puncturing

    Returns size of a single encoded block prior to the puncturing.

  */
  int get_coded_size(void) { return Ncoded; };

  /*!
    \brief Encoder function

    This function can encode (using function from the parent \c Turbo_Codec class) and puncture several consecutive
    coding blocks. Encoded bits are organized as follows:

    \f[ s(1), p_{1,1}(1), p_{1,2}(1), \ldots , p_{1,n_1}(1), p_{2,1}(1), p_{2,2}(1), \ldots , p_{2,n_2}(1), s(2), \ldots \f]

    In the above expression \f$s(n)\f$ is the n-th systematic bit and \f$p_{l,k}(n)\f$ is the n-th bit from the k-th
    encoder polynom of the l-th constituent encoder. A tail of both systematic and parity bits is added after each
    coding block to force both encoder to the zero state. The tail of each encoder is structured as follows (using
    encoder one as an example):

    \f[ t_1(1), pt_{1,1}(1), pt_{1,2}(1), \ldots , pt_{1,n_1}(1), \ldots pt_{1,n_1}(m) \f]

    The tailbits from the first encoder are placed before the tailbits from the second encoder.
    Each encoded block is punctured (independenly from other blocks) according to the pncturing scheme given in the
    \a puncture_matrix

    \param input The input bits to the encoder. Must contain an integer number of code blocks
    \param output The encoded bits including two tail, one from each constituent encoder, after each coding block.
  */
  void encode(const bvec &input, bvec &output);

  /*!
    \brief Encoder function returning vector

    \param input The input bits to the encoder. Must contain an integer number of code blocks

    The same function as above but has only one parameter \a input and returns a vector containing encoded bits.
  */
  bvec encode(const bvec &input);

  /*!
    \brief Decoder function

    This function can decode several consecutive coding blocks that were encoded with the encode member function

    \param received_signal The vector of received bits
    \param decoded_bits A vector of decoded bits
    \param true_bits If this input vector is provided then the iterations will stop as soon as the decoded bits
    equals the \c true_bits. Note that this feature can not be used if the \c in_adaptive_stop parameter in the
    setup function is set to \c true.
  */
  virtual void decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits = "0");

  /*!
    \brief Decoder function returning vector

    This function can decode several consecutive coding blocks that were encoded with the encode member function

    \param received_signal The vector of received bits

    Function returns decoded bits.
  */
  virtual bvec decode(const vec &received_signal);

  /*!
    \brief Decoder function

    This function can decode several consecutive coding blocks that were encoded with the encode member function

    \param received_signal The vector of received bits
    \param decoded_bits A vector of decoded bits
    \param nrof_used_iterations Returns the number of used iterations for each code block.
    \param true_bits If this input vector is provided then the iterations will stop as soon as the decoded bits
    equals the \c true_bits. Note that this feature can not be used if the \c in_adaptive_stop parameter in the
    setup function is set to \c true.
  */
  virtual void decode(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                      const bvec &true_bits = "0");

  /*!
    \brief Calculates length of uncoded block

    Calculates length of the turbo code input block required to produce output block of length \a punctured_size.
    Some puncturing patterns might not allow to create block with such a length and then \a fill_bits of
    dummy bits must be used.

    \param tc Instance of \a Punctured_Turbo_Codec
    \param punctured_size Required size of punctured encoded block
    \param fill_bits Number of dummy bits that must be attached to the encoded block (parameter calculated within
    the function)
  */
  friend int calculate_uncoded_size(Punctured_Turbo_Codec &tc, int punctured_size, int &fill_bits);

protected:
  /*!
    \brief Calculates length of the puncturing codeword

    Calculates length of punctured encoded block based on a \a puncture_matrix and length of the input block
    (size of the turbo code interleaver)
  */
  void calculate_punctured_size(void);

  //Scalars:
  int Period; ///< Number of columns in the puncturing matrix
  long pNcoded;
  int punct_total, punct_total1, punct_total2;
  double rate;
  bmat puncture_matrix;
};


/*!
  \relatesalso Turbo_Codec
  \brief Generates the interleaver sequence for the internal turbo encoder interleaver used in WCDMA
*/
ITPP_EXPORT ivec wcdma_turbo_interleaver_sequence(int interleaver_size);

/*!
  \relatesalso Turbo_Codec
  \author qdelfin and Stephan Ludwig
  \brief Generates the interleaver sequence for the internal turbo encoder interleaver used in LTE
*/
ITPP_EXPORT ivec lte_turbo_interleaver_sequence(int interleaver_size);

} // namespace itpp

#endif // #ifndef TURBO_H
