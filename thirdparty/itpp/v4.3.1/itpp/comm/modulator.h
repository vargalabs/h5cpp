/*!
 * \file
 * \brief One- and two-dimensional modulators - header file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2011  (see AUTHORS file for a list of contributors)
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

#ifndef MODULATOR_H
#define MODULATOR_H

#include <itpp/base/mat.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/min_max.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
 * \ingroup modulators
 * \brief Soft demodulation methods
 */
enum Soft_Method {
  LOGMAP,   //!< Log-MAP full calculation
  APPROX   //!< Approximate faster method
};

/*!
  \ingroup modulators
  \brief General modulator for 1D or 2D signal constellations.

  The Modulator class is designed for modeling any kind of 1D (real) or 2D
  (complex) signal constellations. Therefore it is used as a base class for
  such modulations like PAM, PSK, QAM, etc.

  The constellation of the modulator is described with two vectors. The
  first one contains the real or complex values representing the
  constellation points, whereas the other one includes the corresponding
  bit to symbol mapping (in the decimal from).

  Beside hard demapping, this class can also perform soft demodulation. To
  use it properly the received symbols should be equal to: \f[r_k = c_k
  s_k + n_k,\f] where \f$c_k\f$ is the real or complex channel gain,
  \f$s_k\f$ is the transmitted constellation symbol, and \f$n_k\f$ is the
  AWGN of the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.
*/
template <typename T>
class Modulator
{
public:
  //! Default constructor
  Modulator();
  //! Constructor
  Modulator(const Vec<T>& symbols, const ivec& bits2symbols);
  //! Destructor
  virtual ~Modulator() {}

  //! Set the constellation to use in the modulator
  virtual void set(const Vec<T>& symbols, const ivec& bits2symbols);

  //! Returns number of bits per symbol
  virtual int bits_per_symbol() const { return k; }
  
  //! Returns number of bits per symbol
  virtual int get_k() const { return k; }

  //! Returns number of modulation symbols
  virtual int get_M() const { return M; }

  //! Get the symbol values used in the modulator
  virtual Vec<T> get_symbols() const { return symbols; }

  /*!
   * \brief Get the bitmap, which maps input bits into symbols
   *
   * The mapping is done as follows. An input bit sequence in decimal
   * notation is used for indexing the \c bits2symbols table. The indexing
   * result denotes the symbol to be used from the \c symbols table, e.g.:
   *
   * \code
   * PSK mod(8); // assume 8-PSK modulator
   * cvec sym =  mod.get_symbols();
   * ivec bits2sym = mod.get_bits2symbols();
   * bvec in_bits = "100" // input bits
   * int d = bin2dec(in_bits); // decimal representation of in_bits = 4
   * // mapping of d into PSK symbol using bits2sym and sym tables
   * std::complex<double> out_symbol = sym(bits2sym(d));
   * \endcode
   */
  virtual ivec get_bits2symbols() const { return bits2symbols; }

  //! Modulation of symbols
  virtual void modulate(const ivec& symbolnumbers, Vec<T>& output) const;
  //! Modulation of symbols
  virtual Vec<T> modulate(const ivec& symbolnumbers) const;

  //! Demodulation of symbols
  virtual void demodulate(const Vec<T>& signal, ivec& output) const;
  //! Demodulation of symbols
  virtual ivec demodulate(const Vec<T>& signal) const;

  //! Modulation of bits
  virtual void modulate_bits(const bvec& bits, Vec<T>& output) const;
  //! Modulation of bits
  virtual Vec<T> modulate_bits(const bvec& bits) const;

  //! Hard demodulation of bits
  virtual void demodulate_bits(const Vec<T>& signal, bvec& bits) const;
  //! Hard demodulation of bits
  virtual bvec demodulate_bits(const Vec<T>& signal) const;

  /*!
    \brief Soft demodulator for AWGN channels

    This function calculates the log-likelihood ratio (LLR) of the
    received signal from AWGN channels. Depending on the soft demodulation
    method chosen, either full log-MAP calculation is performed (default
    method), according to the following equation: \f[\log \left(
    \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left( \frac{\sum_{s_i
    \in S_0} \exp \left( -\frac{|r_k - s_i|^2}{N_0} \right)} {\sum_{s_i
    \in S_1} \exp \left( -\frac{|r_k - s_i|^2}{N_0} \right)} \right) \f]
    or approximate, but faster calculation is performed.

    The approximate method finds for each bit the closest constellation
    points that have zero and one in the corresponding position. Let
    \f$d_0 = |r_k - s_0|\f$ denote the distance to the closest zero point
    and \f$d_1 = |r_k - s_1|\f$ denote the distance to the closest one
    point for the corresponding bit respectively. The approximate
    algorithm then computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]

    This function can be used on channels where the channel gain
    \f$c_k = 1\f$.

    When this function is to be used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of
    the turbo decoder shall be set to 1. The output from this function can
    also be used by a Viterbi decoder using an AWGN based metric
    calculation function.

    \param rx_symbols The received noisy constellation symbols
    \param N0 The spectral density of the AWGN noise
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (\c Modulator_ND) instead, which is
    based on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const Vec<T>& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channels
  virtual vec demodulate_soft_bits(const Vec<T>& rx_symbols, double N0,
                                   Soft_Method method = LOGMAP) const;

  /*!
    \brief Soft demodulator for fading channels

    This function calculates the log-likelihood ratio (LLR) of the
    received signal from fading channels. Depending on the soft
    demodulation method chosen, either full log-MAP calculation is
    performed (default method), according to the following equation:
    \f[\log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
    \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r_k - c_k s_i|^2}{N_0}
    \right)} {\sum_{s_i \in S_1} \exp \left( -\frac{|r_k - c_k
    s_i|^2}{N_0} \right)} \right) \f] or approximate, but faster
    calculation is performed.

    The approximate method finds for each bit the closest constellation
    points that have zero and one in the corresponding position. Let
    \f$d_0 = |r_k - c_k s_0|\f$ denote the distance to the closest zero
    point and \f$d_1 = |r_k - c_k s_1|\f$ denote the distance to the
    closest one point for the corresponding bit respectively. The
    approximate algorithm then computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]

    When this function is to be used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of
    the turbo decoder shall be set to 1. The output from this function can
    also be used by a Viterbi decoder using an AWGN based metric
    calculation function.

    \param rx_symbols The received noisy constellation symbols \f$r_k\f$
    \param channel The channel values \f$c_k\f$
    \param N0 The spectral density of the AWGN noise
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const Vec<T>& rx_symbols,
                                    const Vec<T>& channel,
                                    double N0, vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for fading channels
  virtual vec demodulate_soft_bits(const Vec<T>& rx_symbols,
                                   const Vec<T>& channel,
                                   double N0,
                                   Soft_Method method = LOGMAP) const;

protected:
  //! Setup indicator
  bool setup_done;
  //! Number of bits per modulation symbol
  int k;
  //! Number of modulation symbols
  int M;
  //! Bit to symbol mapping table (size: M x k)
  bmat bitmap;
  //! Bit to symbol mapping in decimal form (size: M)
  ivec bits2symbols;
  //! Corresponding modulation symbols (size: M)
  Vec<T> symbols;
  /*! \brief Matrix where row k contains the constellation points with '0'
    in bit position k */
  imat S0;
  /*! \brief Matrix where row k contains the constellation points with '1'
    in bit position k */
  imat S1;

  //! This function calculates the soft bit mapping matrices S0 and S1
  void calculate_softbit_matrices();
};


// ----------------------------------------------------------------------
// Type definitions of Modulator_1D and Modulator_2D
// ----------------------------------------------------------------------

/*!
 * \relates Modulator
 * \brief Definition of 1D Modulator (with real symbols)
 */
typedef Modulator<double> Modulator_1D;

/*!
 * \relates Modulator
 * \brief Definition of 2D Modulator (with complex symbols)
 */
typedef Modulator<std::complex<double> > Modulator_2D;


// ----------------------------------------------------------------------
// Implementation of templated Modulator members
// ----------------------------------------------------------------------

template<typename T>
Modulator<T>::Modulator() :
    setup_done(false), k(0), M(0), bitmap(""), bits2symbols(""), symbols(""),
    S0(""), S1("") {}

template<typename T>
Modulator<T>::Modulator(const Vec<T> &symbols, const ivec &bits2symbols)
{
  set(symbols, bits2symbols);
}

template<typename T>
void Modulator<T>::set(const Vec<T> &in_symbols, const ivec &in_bits2symbols)
{
  it_assert(in_symbols.size() == in_bits2symbols.size(),
            "Modulator<T>::set(): Number of symbols and bits2symbols does not match");
  it_assert(is_even(in_symbols.size()) && (in_symbols.size() > 0),
            "Modulator<T>::set(): Number of symbols needs to be even and non-zero");
  it_assert((max(in_bits2symbols) == in_bits2symbols.size() - 1)
            && (min(in_bits2symbols) == 0), "Modulator<T>::set(): Improper bits2symbol vector");
  symbols = in_symbols;
  bits2symbols = in_bits2symbols;
  M = bits2symbols.size();
  k = levels2bits(M);
  bitmap.set_size(M, k);
  for (int m = 0; m < M; m++) {
    bitmap.set_row(bits2symbols(m), dec2bin(k, m));
  }
  calculate_softbit_matrices();
  setup_done = true;
}


template<typename T>
void Modulator<T>::modulate(const ivec &symbolnumbers, Vec<T>& output) const
{
  it_assert_debug(setup_done, "Modulator<T>::modulate(): Modulator not ready.");
  output.set_size(symbolnumbers.length());
  for (int i = 0; i < symbolnumbers.length(); i++)
    output(i) = symbols(symbolnumbers(i));
}

template<typename T>
Vec<T> Modulator<T>::modulate(const ivec &symbolnumbers) const
{
  Vec<T> output(symbolnumbers.length());
  modulate(symbolnumbers, output);
  return output;
}


template<typename T>
void Modulator<T>::demodulate(const Vec<T> &signal, ivec& output) const
{
  it_assert_debug(setup_done, "Modulator<T>::demodulate(): Modulator not ready.");
  double dist, mindist;
  int closest;
  output.set_size(signal.size());

  for (int i = 0; i < signal.size(); i++) {
    mindist = std::abs(symbols(0) - signal(i));
    closest = 0;
    for (int j = 1; j < M; j++) {
      dist = std::abs(symbols(j) - signal(i));
      if (dist < mindist) {
        mindist = dist;
        closest = j;
      }
    }
    output(i) = closest;
  }
}

template<typename T>
ivec Modulator<T>::demodulate(const Vec<T>& signal) const
{
  ivec output(signal.length());
  demodulate(signal, output);
  return output;
}


template<typename T>
void Modulator<T>::modulate_bits(const bvec &bits, Vec<T> &output) const
{
  it_assert_debug(setup_done, "Modulator<T>::modulate_bits(): Modulator not ready.");
  // Check if some bits have to be cut and print warning message in such
  // case.
  if (bits.length() % k) {
    it_warning("Modulator<T>::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
  }
  int no_symbols = bits.length() / k;
  output.set_size(no_symbols);
  for (int i = 0; i < no_symbols; i++) {
    output(i) = symbols(bits2symbols(bin2dec(bits.mid(i * k, k))));
  }
}

template<typename T>
Vec<T> Modulator<T>::modulate_bits(const bvec &bits) const
{
  Vec<T> output;
  modulate_bits(bits, output);
  return output;
}

template<typename T>
void Modulator<T>::demodulate_bits(const Vec<T> &signal, bvec &bits) const
{
  it_assert_debug(setup_done, "Modulator<T>::demodulate_bist(): Modulator not ready.");
  double dist, mindist;
  int closest;
  bits.set_size(k*signal.size());

  for (int i = 0; i < signal.size(); i++) {
    mindist = std::abs(symbols(0) - signal(i));
    closest = 0;
    for (int j = 1; j < M; j++) {
      dist = std::abs(symbols(j) - signal(i));
      if (dist < mindist) {
        mindist = dist;
        closest = j;
      }
    }
    bits.replace_mid(i*k, bitmap.get_row(closest));
  }
}

template<typename T>
bvec Modulator<T>::demodulate_bits(const Vec<T> &signal) const
{
  bvec bits;
  demodulate_bits(signal, bits);
  return bits;
}


template<typename T>
void Modulator<T>::demodulate_soft_bits(const Vec<T> &rx_symbols, double N0,
                                        vec &soft_bits,
                                        Soft_Method method) const
{
  it_assert_debug(setup_done, "Modulator<T>::demodulate_soft_bits(): Modulator not ready.");
  double P0, P1, d0min, d1min, temp;
  vec metric(M);

  soft_bits.set_size(k * rx_symbols.size());

  if (method == LOGMAP) {
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = std::exp(-sqr(rx_symbols(l) - symbols(j)) / N0);
      }
      for (int i = 0; i < k; i++) {
        P0 = P1 = 0;
        for (int j = 0; j < (M >> 1); j++) {
          P0 += metric(S0(i, j));
          P1 += metric(S1(i, j));
        }
        soft_bits(l*k + i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }
  else { // method == APPROX
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = sqr(rx_symbols(l) - symbols(j));
      }
      for (int i = 0; i < k; i++) {
        d0min = d1min = std::numeric_limits<double>::max();
        for (int j = 0; j < (M >> 1); j++) {
          temp = metric(S0(i, j));
          if (temp < d0min) { d0min = temp; }
          temp = metric(S1(i, j));
          if (temp < d1min) { d1min = temp; }
        }
        soft_bits(l*k + i) = (-d0min + d1min) / N0;
      }
    }
  }
}

template<typename T>
vec Modulator<T>::demodulate_soft_bits(const Vec<T> &rx_symbols,
                                       double N0,
                                       Soft_Method method) const
{
  vec output;
  demodulate_soft_bits(rx_symbols, N0, output, method);
  return output;
}

template<typename T>
void Modulator<T>::demodulate_soft_bits(const Vec<T> &rx_symbols,
                                        const Vec<T> &channel, double N0,
                                        vec &soft_bits,
                                        Soft_Method method) const
{
  it_assert_debug(setup_done, "Modulator_2D::demodulate_soft_bits(): Modulator not ready.");
  double P0, P1, d0min, d1min, temp;
  vec metric(M);

  soft_bits.set_size(k * rx_symbols.size());

  if (method == LOGMAP) {
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = std::exp(-sqr(rx_symbols(l) - channel(l) * symbols(j))
                             / N0);
      }
      for (int i = 0; i < k; i++) {
        P0 = P1 = 0;
        for (int j = 0; j < (M >> 1); j++) {
          P0 += metric(S0(i, j));
          P1 += metric(S1(i, j));
        }
        soft_bits(l*k + i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }
  else { // method == APPROX
    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
        metric(j) = sqr(rx_symbols(l) - channel(l) * symbols(j));
      }
      for (int i = 0; i < k; i++) {
        d0min = d1min = std::numeric_limits<double>::max();
        for (int j = 0; j < (M >> 1); j++) {
          temp = metric(S0(i, j));
          if (temp < d0min) { d0min = temp; }
          temp = metric(S1(i, j));
          if (temp < d1min) { d1min = temp; }
        }
        soft_bits(l*k + i) = (-d0min + d1min) / N0;
      }
    }
  }
}

template<typename T>
vec Modulator<T>::demodulate_soft_bits(const Vec<T> &rx_symbols,
                                       const Vec<T> &channel,
                                       double N0,
                                       Soft_Method method) const
{
  vec output;
  demodulate_soft_bits(rx_symbols, channel, N0, output, method);
  return output;
}

template<typename T>
void Modulator<T>::calculate_softbit_matrices()
{
  int count0, count1;

  // Allocate storage space for the result matrices:
  S0.set_size(k, M >> 1, false);
  S1.set_size(k, M >> 1, false);

  for (int i = 0; i < k; i++) {
    count0 = 0;
    count1 = 0;
    for (int j = 0; j < M; j++) {
      if (bitmap(j, i) == bin(0)) {
        S0(i, count0++) = j;
      }
      else {
        S1(i, count1++) = j;
      }
    }
  }
}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Modulator<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Modulator<std::complex<double> >;

//! \endcond

// ----------------------------------------------------------------------
// QAM : Modulator_2D
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief M-ary QAM modulator with square lattice.

  The size of the QAM constellation is \f$M = 2^k\f$, where \f$k = 1, 2,
  \ldots \f$. Symbol values in each dimension are: \f$\{-(\sqrt{M}-1),
  \ldots, -3, -1, 1, 3, \ldots, (\sqrt{M}-1)\}\f$. The bitmap is Gray
  encoded. Symbols are normalized so that the average energy is 1. That
  is, normalized with \f$\sqrt{2(M-1)/3}\f$.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.
*/
class ITPP_EXPORT QAM : public Modulator<std::complex<double> >
{
public:
  //! Default Constructor
  QAM() {}
  //! Class Constructor
  QAM(int M) { set_M(M); }
  //! Destructor
  virtual ~QAM() { }
  //! Change the size of the signal constellation
  void set_M(int M);

  //! Hard demodulation of bits
  void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  bvec demodulate_bits(const cvec& signal) const;

protected:
  //! The square-root of M
  int L;
  //! Scaling factor of square QAM constellation (sqrt((M-1)*2/3))
  double scaling_factor;
};


// ----------------------------------------------------------------------
// PSK : Modulator<std::complex<double> >
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief M-ary PSK modulator.

  This class implements the M-ary PSK modulator with \f$M = 2^k\f$
  constellation points, where \f$k = 1, 2, \ldots \f$. The symbol
  numbering is counter clockwise starting from the real axis, i.e. symbol
  \f$(1, 0)\f$. The bitmap is Gray encoded. The symbol energy is
  normalized to 1.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.
*/
class ITPP_EXPORT PSK : public Modulator<std::complex<double> >
{
public:
  //! Default Constructor
  PSK() {}
  //! Class constructor
  PSK(int M) { set_M(M); }
  //! Destructor
  virtual ~PSK() { }
  //! Change the size of the signal constellation
  void set_M(int M);

  //! Hard demodulation of bits
  void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  bvec demodulate_bits(const cvec& signal) const;
};


// ----------------------------------------------------------------------
// QPSK : PSK : Modulator<std::complex<double> >
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief QPSK modulator.

  This is a special version of the PSK modulator with \f$M = 4\f$
  constellation points. Symbol numbering is counter clockwise starting
  from the real axis. Bits are Gray coded onto symbols. Symbol energy is
  normalized to 1.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.
*/
class ITPP_EXPORT QPSK : public PSK
{
public:
  //! Class Constructor
  QPSK(): PSK(4) {}
  //! Destructor
  virtual ~QPSK() {}

  /*!
    \brief Soft demodulator for AWGN channel

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{2 \sqrt{2}}{N_0} \Im\{r_k \exp \left(j \frac{\Pi}{4} \right)
    \}\f] and \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) = \frac{2
    \sqrt{2}}{N_0} \Re\{r_k \exp \left(j \frac{\Pi}{4} \right) \}\f]
    depending on the bit positon in the QPSK symbol.

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channel
  vec demodulate_soft_bits(const cvec& rx_symbols, double N0,
                           Soft_Method method = LOGMAP) const;


  /*!
    \brief Soft demodulator for a known channel in AWGN

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{2 \sqrt{2}}{N_0} \Im\{r_k c_k \exp \left(j \frac{\Pi}{4} \right)
    \}\f] and \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) = \frac{2
    \sqrt{2}}{N_0} \Re\{r_k c_k \exp \left(j \frac{\Pi}{4} \right) \}\f]
    depending on the bit positon in the QPSK symbol.

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    \param channel The channel coefficients, \f$c\f$
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols,
                                    const cvec& channel, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for a known channel in AWGN
  vec demodulate_soft_bits(const cvec& rx_symbols, const cvec& channel,
                           double N0, Soft_Method method = LOGMAP) const;
};


// ----------------------------------------------------------------------
// BPSK_c : PSK : Modulator<std::complex<double> >
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief BPSK modulator with complex symbols.

  This is a special version of the PSK modulator with \f$M = 2\f$
  constellation points. The following bit to symbol mapping is used:
  - \f$0 \rightarrow 1+0i\f$
  - \f$1 \rightarrow -1+0i\f$.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.

  \note Although constellation points of the BPSK modulator can be
  represented in the real domain only, this class uses complex signals to
  be compatible with other PSK and QAM based modulators.

  \sa BPSK
*/
class ITPP_EXPORT BPSK_c : public PSK
{
public:
  //! Constructor
  BPSK_c(): PSK(2) {}
  //! Destructor
  virtual ~BPSK_c() {}

  //! Modulate bits into BPSK symbols in complex domain
  void modulate_bits(const bvec& bits, cvec& output) const;
  //! Modulate bits into BPSK symbols  in complex domain
  cvec modulate_bits(const bvec& bits) const;
  //! Demodulate noisy BPSK symbols in complex domain into bits
  void demodulate_bits(const cvec& signal, bvec& output) const;
  //! Demodulate noisy BPSK symbols in complex domain into bits
  bvec demodulate_bits(const cvec& signal) const;

  /*!
    \brief Soft demodulator for AWGN channel

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{4 \Re\{r\}} {N_0}\f]

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    (complex but symbols in real part)
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channel
  vec demodulate_soft_bits(const cvec& rx_symbols, double N0,
                           Soft_Method method = LOGMAP) const;

  /*!
    \brief Soft demodulator for a known channel in AWGN

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{4 \Re\{r c^{*}\}}{N_0}\f]

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    (complex but symbols in real part)
    \param channel The channel coefficients, \f$c\f$ (complex)
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols,
                                    const cvec& channel, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for a known channel in AWGN
  vec demodulate_soft_bits(const cvec& rx_symbols, const cvec& channel,
                           double N0, Soft_Method method = LOGMAP) const;
};



// ----------------------------------------------------------------------
// BPSK : Modulator<double>
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief BPSK modulator with real symbols.

  This is a special version of the PSK modulator with \f$M = 2\f$
  constellation points. The following bit to symbol mapping is used:
  - \f$0 \rightarrow 1\f$
  - \f$1 \rightarrow -1\f$.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.

  \note This class uses real values for representing symbols. There is
  a similar class named BPSK_c, which uses complex values for symbols and
  therefore is compatible with other PSK and QAM based modulators.
*/
class ITPP_EXPORT BPSK : public Modulator<double>
{
public:
  //! Constructor
  BPSK(): Modulator<double>("1.0 -1.0", "0 1") {}
  //! Destructor
  virtual ~BPSK() {}

  //! Modulate bits into BPSK symbols in complex domain
  void modulate_bits(const bvec& bits, vec& output) const;
  //! Modulate bits into BPSK symbols  in complex domain
  vec modulate_bits(const bvec& bits) const;
  //! Demodulate noisy BPSK symbols in complex domain into bits
  void demodulate_bits(const vec& signal, bvec& output) const;
  //! Demodulate noisy BPSK symbols in complex domain into bits
  bvec demodulate_bits(const vec& signal) const;

  /*!
    \brief Soft demodulator for AWGN channel

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{4 r}{N_0}\f]

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the N-dimensional
    modulator (Modulator_ND) instead, which is based on the QLLR
    (quantized) arithmetic and therefore is faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const vec& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channel
  vec demodulate_soft_bits(const vec& rx_symbols, double N0,
                           Soft_Method method = LOGMAP) const;

  /*!
    \brief Soft demodulator for a known channel in AWGN

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{4 \Re\{r c^{*}\}}{N_0}\f]

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    (complex but symbols in real part)
    \param channel The channel coefficients, \f$c\f$ (complex)
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the N-dimensional
    modulator (Modulator_ND) instead, which is based on the QLLR
    (quantized) arithmetic and therefore is faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const vec& rx_symbols,
                                    const vec& channel, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for a known channel in AWGN
  vec demodulate_soft_bits(const vec& rx_symbols, const vec& channel,
                           double N0, Soft_Method method = LOGMAP) const;
};


// ----------------------------------------------------------------------
// PAM_c : Modulator<std::complex<double> >
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief M-ary PAM modulator with complex symbols.

  This class implements an M-ary PAM modulator with the following signal
  values: \f$\{-(M-1), \ldots, -3, -1, 1, 3, \ldots, (M-1)\}\f$. Symbol
  numbering is from right to left in the increasing order. The Gray
  encoding of bits to symbols is used.

  The constellation symbols are normalized so that the average energy is
  equal to 1. That is, normalized with \f$ \sqrt{(M^2-1)/3}\f$.

  \note Although the constellation points can be represented in the real
  domain only, this class uses complex based interface to be compatible
  with other PSK and QAM based modulators.

  \sa PAM
*/
class ITPP_EXPORT PAM_c : public Modulator<std::complex<double> >
{
public:
  //! Default Constructor
  PAM_c() {}
  //! Constructor
  PAM_c(int M) { set_M(M); }
  //! Destructor
  virtual ~PAM_c() {}
  //! Set the size of the signal constellation
  void set_M(int M);

  //! Hard demodulation of PAM symbols in complex domain to bits
  void demodulate_bits(const cvec& signal, bvec& output) const;
  //! Hard demodulation of PAM symbols in complex domain to bits
  bvec demodulate_bits(const cvec& signal) const;

  /*!
    \brief Soft demodulator for AWGN channels.

    This function calculates the log-likelihood ratio (LLR) of the
    received signal from AWGN channels. Depending on the soft demodulation
    method chosen, either full log-MAP calculation is performed (default
    method), according to the following equation: \f[\log \left(
    \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left( \frac{\sum_{s_i
    \in S_0} \exp \left( -\frac{|r_k - s_i|^2}{N_0} \right)} {\sum_{s_i
    \in S_1} \exp \left( -\frac{|r_k - s_i|^2}{N_0} \right)} \right) \f]
    or approximate, but faster calculation is performed.

    The approximate method finds for each bit the closest constellation
    points that have zero and one in the corresponding position. Let
    \f$d_0 = |r_k - s_0|\f$ denote the distance to the closest zero point
    and \f$d_1 = |r_k - s_1|\f$ denote the distance to the closest one
    point for the corresponding bit respectively. The approximate
    algorithm then computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]

    This function can be used on channels where the channel gain
    \f$c = 1\f$.

    When this function is to be used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of
    the turbo decoder shall be set to 1. The output from this function can
    also be used by a Viterbi decoder using an AWGN based metric
    calculation function.

    \param rx_symbols The received noisy constellation symbols \f$r_k\f$
    (complex, but symbols are real)
    \param N0 The spectral density of the AWGN noise
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the N-dimensional
    modulator (Modulator_ND) instead, which is based on the QLLR
    (quantized) arithmetic and therefore is faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channels.
  virtual vec demodulate_soft_bits(const cvec& rx_symbols, double N0,
                                   Soft_Method method = LOGMAP) const;

  /*!
    \brief Soft demodulator for known fading channels.

    This function calculates the log-likelihood ratio (LLR) of the
    received signal from fading channels. Depending on the soft
    demodulation method chosen, either full log-MAP calculation is
    performed (default method), according to the following equation:
    \f[\log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
    \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r_k - c_k s_i|^2}{N_0}
    \right)} {\sum_{s_i \in S_1} \exp \left( -\frac{|r_k - c_k
    s_i|^2}{N_0} \right)} \right) \f] or approximate, but faster
    calculation is performed.

    The approximate method finds for each bit the closest constellation
    points that have zero and one in the corresponding position. Let
    \f$d_0 = |r_k - c_k s_0|\f$ denote the distance to the closest zero
    point and \f$d_1 = |r_k - c_k s_1|\f$ denote the distance to the
    closest one point for the corresponding bit respectively. The
    approximate algorithm then computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]

    When this function is to be used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of
    the turbo decoder shall be set to 1. The output from this function can
    also be used by a Viterbi decoder using an AWGN based metric
    calculation function.

    \param rx_symbols The received noisy constellation symbols \f$r_k\f$
    (complex)
    \param channel The channel values \f$c_k\f$
    \param N0 The spectral density of the AWGN noise
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the N-dimensional
    modulator (Modulator_ND) instead, which is based on the QLLR
    (quantized) arithmetic and therefore is faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols,
                                    const cvec& channel, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for known fading channels.
  virtual vec demodulate_soft_bits(const cvec& rx_symbols,
                                   const cvec& channel, double N0,
                                   Soft_Method method = LOGMAP) const;

protected:
  //! Scaling factor used to normalize the average energy to 1
  double scaling_factor;
};


// ----------------------------------------------------------------------
// PAM : Modulator<double>
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief M-ary PAM modulator with real symbols.

  This class implements an M-ary PAM modulator with the following signal
  values: \f$\{-(M-1), \ldots, -3, -1, 1, 3, \ldots, (M-1)\}\f$. Symbol
  numbering is from right to left in the increasing order. The Gray
  encoding of bits to symbols is used.

  The constellation symbols are normalized so that the average energy is
  equal to 1. That is, normalized with \f$ \sqrt{(M^2-1)/3}\f$.

  \note This class uses real values for representing symbols. There is
  a similar class named PAM_c, which uses complex values for symbols and
  therefore is compatible with other PSK and QAM based modulators.
*/
class ITPP_EXPORT PAM : public Modulator<double>
{
public:
  //! Default Constructor
  PAM() {}
  //! Constructor
  PAM(int M) { set_M(M); }
  //! Destructor
  virtual ~PAM() {}
  //! Set the size of the signal constellation
  void set_M(int M);

  //! Hard demodulation of PAM symbols in complex domain to bits
  void demodulate_bits(const vec& signal, bvec& output) const;
  //! Hard demodulation of PAM symbols in complex domain to bits
  bvec demodulate_bits(const vec& signal) const;

protected:
  //! Scaling factor used to normalize the average energy to 1
  double scaling_factor;
};

} // namespace itpp

#endif // #ifndef MODULATOR_H
