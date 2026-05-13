/*!
 * \file
 * \brief Definition of vector (MIMO) modulator classes
 * \author Mirsad Cirkic, Erik G. Larsson and Adam Piatyszek
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

#ifndef MODULATOR_ND_H
#define MODULATOR_ND_H

#include <itpp/base/vec.h>
#include <itpp/base/array.h>
#include <itpp/comm/llr.h>
#include <itpp/itexports.h>
#include <itpp/base/base_exports.h>

namespace itpp
{

/*!
 * \addtogroup modulators
 */

// ----------------------------------------------------------------------
// Modulator_ND
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * \brief Base class for an N-dimensional (ND) vector (MIMO) modulator.
 *
 * See \c ND_UPAM class for examples.
 *
 * \note Can also be used for scalar modulation/demodulation as an
 * alternative to \c Modulator_1D or \c Modulator_2D. Mixed use of \c
 * Modulator_1D or \c Modulator_2D and \c Modulator_ND is <b>not
 * advised</b>.
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT Modulator_ND
{
public:
  //! Soft demodulation method
  enum Soft_Demod_Method {
    //! Log-MAP demodulation by "brute-force" enumeration of all points
    FULL_ENUM_LOGMAP,
    //! Max-Log demodulation by "brute-force" enumeration of all points
    FULL_ENUM_MAXLOG,
    //! Zero-Forcing Log-MAP approximated demodulation
    ZF_LOGMAP
  };

  //! Default constructor
  Modulator_ND(LLR_calc_unit llrcalc_in = LLR_calc_unit()):
    llrcalc(llrcalc_in), demod_initialized(false) {}
  //! Destructor
  virtual ~Modulator_ND() {}

  //! Set LLR calculation unit
  void set_llrcalc(LLR_calc_unit llrcalc_in) {
    llrcalc = llrcalc_in;
  };

  //! Get LLR calculation unit
  LLR_calc_unit get_llrcalc() const {
    return llrcalc;
  }

  //! Get number of dimensions
  int get_dim() const {
    return nt;
  }

  //! Get number of bits per modulation symbol per dimension
  ivec get_k() const {
    return k;
  }

  //! Get number of bits per modulation symbol per dimension
  ivec bits_per_symbol() const {
    return k;
  }

  //! Get number of modulation symbols per dimension
  ivec get_M() const {
    return M;
  }

  //! Get bit pattern in decimal
  Array<ivec> get_bits2symbols() const {
    return bits2symbols;
  }

  //! Get Bit mapping table
  Array<bmat> get_bitmap() const {
    return bitmap;
  }

protected:
  //! Number of dimensions
  int nt;
  //! Number of bits in the symbol vector
  int nb;
  //! LLR calculation unit
  LLR_calc_unit llrcalc;
  //! Number of bits per modulation symbol
  ivec k;
  //! Number of modulation symbols along each dimension
  ivec M;
  //! Flag indicating whether the demodulator has been initialized
  bool demod_initialized;
  //! Bit mapping table (one table per dimension)
  Array<bmat> bitmap;
  //! Bit pattern in decimal form ordered and the corresponding symbols (one pattern per dimension)
  Array<ivec> bits2symbols;
  //! The normalization factor in the exponent (in front of the square norm) in the Gaussian distribution
  double gaussnorm;  
  //! Norms part dependent on H
  itpp::vec hnorms;
  //! Norms part depending on both H and y
  itpp::QLLRvec Qnorms;
  //! A prioi information
  itpp::QLLRvec llrapr;
  //! The bit to column mapping
  itpp::ivec bpos2cpos;
  //! The cumulative sum of bits in the symbol vector
  itpp::ivec bitcumsum;
  //! The Gray to decimal mapping
  itpp::Array<itpp::Vec<unsigned> > gray2dec;
  //! Convert LLR to log-probabilities
  QLLRvec probabilities(QLLR l); // some abuse of what QLLR stands for...
  //! Convert LLR to log-probabilities, vector version
  Array<QLLRvec> probabilities(const QLLRvec &l);
  //! Marginalize (sum) over the bits
  void marginalize_bits(itpp::QLLRvec& llr, Soft_Demod_Method method) const;
  //! Hardcoded implementation of 1:st bit demodulation
  void demodllrbit0(itpp::QLLR& llr) const;
  //! Hardcoded implementation of 2:nd bit demodulation
  void demodllrbit1(itpp::QLLR& llr) const;
  //! Hardcoded implementation of 3:rd bit demodulation
  void demodllrbit2(itpp::QLLR& llr) const;
  //! Hardcoded implementation of 1:st bit demodulation
  void demodmaxbit0(itpp::QLLR& maxllr) const;
  //! Hardcoded implementation of 2:nd bit demodulation
  void demodmaxbit1(itpp::QLLR& maxllr) const;
  //! Hardcoded implementation of 3:rd bit demodulation
  void demodmaxbit2(itpp::QLLR& maxllr) const;

  /*!
   * \brief Update LLR (for internal use)
   *
   * This function updates the numerator and denominator in the expression
   * \f[
   * \log \left( \frac {\sum_{s:b_k=0} \exp(-x^2) P(s)}
   *                   {\sum_{s:b_k=1} \exp(-x^2) P(s)} \right)
   * \f]
   *
   * \param[in]   logP_apriori  Vector of a priori probabilities per bit
   * \param[in]   s             Symbol vector
   * \param[in]   scaled_norm   Argument of the exponents in the above
   *                            equation
   * \param[out]  num           Logarithm of the numerator in the above
   *                            expression
   * \param[out]  denom         Logarithm of the denominator in the above
   *                            expression
   */
  void update_LLR(const Array<QLLRvec> &logP_apriori, const ivec &s,
                  QLLR scaled_norm, QLLRvec &num, QLLRvec &denom);

  /*!
   * \brief Update LLR, for scalar channel (for internal use)
   *
   * This function updates the numerator and denominator in the expression
   * \f[
   * \log \left( \frac {\sum_{s:b_k=0} \exp (-x^2) P(s)}
   *                   {\sum_{s:b_k=1} \exp (-x^2) P(s)} \right)
   * \f]
   *
   * \param[in]   logP_apriori  Vector of a priori probabilities per bit
   * \param[in]   s             Symbol
   * \param[in]   scaled_norm   Argument of the exponents in the above
   *                            equation
   * \param[in]   j             Channel index (dimension)
   * \param[out]  num           Logarithm of the numerator in the above
   *                            expression
   * \param[out]  denom         Logarithm of the denominator in the above
   *                            expression
  */
  void update_LLR(const Array<QLLRvec> &logP_apriori, int s,
                  QLLR scaled_norm, int j, QLLRvec &num, QLLRvec &denom);
};


// ----------------------------------------------------------------------
// Modulator_NRD
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * \brief Base class for N-dimensional vector (MIMO) channel
 * modulators/demodulators with real-valued components.
 *
 * This class can be used to perform modulation and demodulation for a
 * matrix (MIMO) channel of the form
 * \f[ y = Hx+e \f],
 * where H is the channel matrix of dimension \f$n_r\times n_t\f$, \f$y\f$
 * is a received vector of length \f$n_r\f$, \f$x\f$ is a transmitted vector
 * of length \f$n_t\f$ and \f$e\f$ is a noise vector.
 *
 * The class supports soft-input soft-output demodulation. It can also be
 * used for scalar modulation to take advantage of this feature.
 *
 * Complex MIMO channels can be handled by using the \c Modulator_NCD
 * class. Alternatively, if the signal constellation is separable in I/Q
 * then the complex channel can be first transformed to a real channel
 * \f[
 * G = \left[ \begin{array}{cc} H_r & -H_i \\ H_i & H_r \end{array} \right]
 * \f]
 *
 * See \c ND_UPAM for examples.
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT Modulator_NRD : public Modulator_ND
{
public:
  //! Constructor
  Modulator_NRD() {}
  //! Destructor
  virtual ~Modulator_NRD() {}

  //! Get modulation symbols per dimension
  Array<vec> get_symbols() const;

  //! Modulate \c bits into \c symbols
  void modulate_bits(const bvec &bits, vec &symbols) const;

  //! Modulate \c bits vector. Symbols are returned.
  vec modulate_bits(const bvec &bits) const;

	  /*!
  * \brief Soft MAP demodulation for multidimensional channel, by
  * "brute-force" enumeration of all constellation points.
  *
  * This function precomputes the norms
  * \f[\frac{|Hs|^2}{2\sigma^2}\f]
  * used to compute the LLR values
  * \f[
  * LLR(k) = \log \left( \frac
  * {\sum_{s:b_k=0} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * {\sum_{s:b_k=1} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * \right)
  * \f]
  *
  * without approximations. It is assumed that H is
  * real-valued. Complex-valued channels can be handled using the \c
  * Modulator_NCD class.
  */
  void init_soft_demodulator(const itpp::mat& H, const double& sigma2);

	
  /*!
  * \brief Soft MAP demodulation for multidimensional channel, by
  * "brute-force" enumeration of all constellation points.
  *
  * This function computes the LLR values
  * \f[
  * LLR(k) = \log \left( \frac
  * {\sum_{s:b_k=0} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * {\sum_{s:b_k=1} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * \right)
  * \f]
  *
  * without approximations. It is assumed that H, y and s are
  * real-valued. Complex-valued channels can be handled using the \c
  * Modulator_NCD class. Currently the following two demodulation methods
  * are supported:
  * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
  *   enumeration of all constellation points
  * - FULL_ENUM_MAXLOG - max-log approximate demodulation, which use "brute-force"
  *   enumeration to find the constellation points that give the smallest euclidian
  *   distances
  *
  * \param[in]   y                Received vector
  *                               (typically \f$N_0/2\f$)
  * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
  * \param[out]  LLR_aposteriori  Vector of a posteriori LLR values
  * \param[in]   method           Soft demodulation method
  *
  * The function performs an exhaustive search over all possible points
  * \c s in the n-dimensional constellation. This is only feasible for
  * relatively small constellations. The Jacobian logarithm is used to
  * compute the sum-exp expression.
  */
  void demodulate_soft_bits(const vec &y,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori,
                            Soft_Demod_Method method = FULL_ENUM_LOGMAP);

  /*!
   * \brief Soft demodulation wrapper function for various methods
   *
   * Currently the following three demodulation methods are supported:
   * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
   *   enumeration of all constellation points
   * - FULL_ENUM_MAXLOG - max-log approximate demodulation, which use "brute-force"
   *   enumeration to find the constellation points that give the smallest euclidian
   *   distances
   * - ZF_LOGMAP - approximated methods with Zero-Forcing preprocessing,
   *   which sometimes tends to perform poorly, especially for poorly
   *   conditioned H
   *
   * \param[in]   y                Received vector
   * \param[in]   H                Channel matrix
   * \param[in]   sigma2           Noise variance per real dimension
   *                               (typically \f$N_0/2\f$)
   * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
   * \param[out]  LLR_aposteriori  Vector of a posteriori LLR values
   * \param[in]   method           Soft demodulation method
   */
  void demodulate_soft_bits(const vec &y, const mat &H, double sigma2,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori,
                            Soft_Demod_Method method = FULL_ENUM_LOGMAP);


  /*!
   * \brief Soft demodulation wrapper function for various methods
   *
   * Currently the following two demodulation methods are supported:
   * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
   *   enumeration of all constellation points
   * - ZF_LOGMAP - approximated methods with Zero-Forcing preprocessing,
   *   which sometimes tends to perform poorly, especially for poorly
   *   conditioned H
   *
   * \param[in]   y                Received vector
   * \param[in]   H                Channel matrix
   * \param[in]   sigma2           Noise variance per real dimension
   *                               (typically \f$N_0/2\f$)
   * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
   * \param[in]   method           Soft demodulation method
   * \return                       Vector of a posteriori LLR values
   */
  QLLRvec demodulate_soft_bits(const vec &y, const mat &H, double sigma2,
                               const QLLRvec &LLR_apriori,
                               Soft_Demod_Method method = FULL_ENUM_LOGMAP);


  /*!
   * \brief Soft MAP demodulation for parallelchannels without crosstalk.
   *
   * This function is a much faster equivalent to \c demodulate_soft_bits
   * with \f$H = \mbox{diag}(h)\f$. Its complexity is linear in the number
   * of subchannels.
   */
  void demodulate_soft_bits(const vec &y, const vec &h, double sigma2,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori);

  //! Output some properties of the MIMO modulator (mainly to aid debugging)
  friend ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const Modulator_NRD &m);

protected:
  //! Vectors of modulation symbols (along each dimension)
  Array<vec> symbols;

  /*!
   * \brief Update residual norm (for internal use).
   *
   * Update the residual norm \f$|y-Hs|\f$ when moving from one
   * constellation point to an adjacent point.
   *
   * \param[in,out]  norm  Norm to be updated
   * \param[in]      k     Position where s changed
   * \param[in]      sold  Old value of s[k]
   * \param[in]      snew  New value of s[k]
   * \param[in]      ytH   y'H vector
   * \param[in]      HtH   Grammian matrix H'H
   * \param[in]      s     Symbol vector
   */
  void update_norm(double &norm, int k, int sold, int snew, const vec &ytH,
                   const mat &HtH, const ivec &s);

  //! Calculation of the part of the norms that depends on H
  void hxnormupdate(itpp::vec& Hx, unsigned& bitstring, unsigned& ind, unsigned bit);

  //! Calculation of the remaining part of the norms that depends both on H and y
  void yxnormupdate(double& yx, itpp::QLLR& lapr, unsigned& bitstring, unsigned& ind, unsigned bit);

	//! Real channel matrix
	itpp::mat H;
	//! The spacing between different constellation points multiplied by the different H columns
	itpp::Array<itpp::Array<itpp::vec> > hspacings;
	//! The spacing between different constellation points scaled by different y elements
  itpp::Array<itpp::vec> yspacings;
};

/*!
 * \relatesalso Modulator_NRD
 * \brief Print some properties of the MIMO modulator (mainly to aid debugging)
 */
ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const Modulator_NRD &m);


// ----------------------------------------------------------------------
// Modulator_NCD
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * \brief Base class for vector (MIMO) channel modulator/demodulators
 * with complex valued components.
 *
 * This class is equivalent to \c Modulator_NRD except for that all
 * quantities are complex-valued.
 *
 * See \c ND_UPAM for examples.
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT Modulator_NCD : public Modulator_ND
{
public:
  //! Constructor
  Modulator_NCD() {}
  //! Destructor
  virtual ~Modulator_NCD() {}

  //! Get modulation symbols per dimension
  Array<cvec> get_symbols() const;

  //! Modulate \c bits into \c symbols
  void modulate_bits(const bvec &bits, cvec &symbols) const;

  //! Modulation of bits
  cvec modulate_bits(const bvec &bits) const;

	/*!
  * \brief Soft MAP demodulation for multidimensional channel, by
  * "brute-force" enumeration of all constellation points.
  *
  * This function computes the norms
  * \f[\frac{|y - Hs|^2}{2\sigma^2}\f]
  * used to compute the LLR values
  * \f[
  * LLR(k) = \log \left( \frac
  * {\sum_{s:b_k=0} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * {\sum_{s:b_k=1} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * \right)
  * \f]
  *
  * without approximations. It is assumed that H is
  * real-valued. Complex-valued channels can be handled using the \c
  * Modulator_NCD class.
  */
  void init_soft_demodulator(const itpp::cmat& H, const double& sigma2);

  /*!
  * \brief Soft MAP demodulation for multidimensional channel, by
  * "brute-force" enumeration of all constellation points.
  *
  * This function computes the LLR values
  * \f[
  * LLR(k) = \log \left( \frac
  * {\sum_{s:b_k=0} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * {\sum_{s:b_k=1} \exp \left( -\frac{|y - Hs|^2}{2\sigma^2} \right) P(s)}
  * \right)
  * \f]
  *
  * without approximations. Currently the following two demodulation methods
  * are supported:
  * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
  *   enumeration of all constellation points
  * - FULL_ENUM_MAXLOG - max-log approximate demodulation, which use "brute-force"
  *   enumeration to find the constellation points that give the smallest euclidian
  *   distances
  *
  * \param[in]   y                Received vector
  *                               (typically \f$N_0/2\f$)
  * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
  * \param[out]  LLR_aposteriori  Vector of a posteriori LLR values
  * \param[in]   method           Soft demodulation method
  *
  * The function performs an exhaustive search over all possible points
  * \c s in the n-dimensional constellation. This is only feasible for
  * relatively small constellations. The Jacobian logarithm is used to
  * compute the sum-exp expression.
  */
  void demodulate_soft_bits(const cvec &y,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori,
                            Soft_Demod_Method method = FULL_ENUM_LOGMAP);

  //! Soft demodulation wrapper function for various methods
  /*!
   * \brief Soft demodulation wrapper function for various methods
   *
   * Currently the following three demodulation methods are supported:
   * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
   *   enumeration of all constellation points
	* - FULL_ENUM_MAXLOG - max-log approximate demodulation, which use "brute-force"
	*   enumeration to find the constellation points that give the smallest euclidian
	*   distances
   * - ZF_LOGMAP - approximated methods with Zero-Forcing preprocessing,
   *   which sometimes tends to perform poorly, especially for poorly
   *   conditioned H
   *
   * \param[in]   y                Received vector
   * \param[in]   H                Channel matrix
   * \param[in]   sigma2           Noise variance per complex dimension,
   *                               i.e. the sum of real and imaginary parts
   *                               (typically \f$N_0\f$)
   * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
   * \param[out]  LLR_aposteriori  Vector of a posteriori LLR values
   * \param[in]   method           Soft demodulation method
   */
  void demodulate_soft_bits(const cvec &y, const cmat &H, double sigma2,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori,
                            Soft_Demod_Method method = FULL_ENUM_LOGMAP);

  /*!
   * \brief Soft demodulation wrapper function for various methods
   *
   * Currently the following three demodulation methods are supported:
   * - FULL_ENUM_LOGMAP - exact demodulation, which use "brute-force"
   *   enumeration of all constellation points
	* - FULL_ENUM_MAXLOG - max-log approximate demodulation, which use "brute-force"
	*   enumeration to find the constellation points that give the smallest euclidian
	*   distances
   * - ZF_LOGMAP - approximated methods with Zero-Forcing preprocessing,
   *   which sometimes tends to perform poorly, especially for poorly
   *   conditioned H
   *
   * \param[in]   y                Received vector
   * \param[in]   H                Channel matrix
   * \param[in]   sigma2           Noise variance per complex dimension,
   *                               i.e. the sum of real and imaginary parts
   *                               (typically \f$N_0\f$)
   * \param[in]   LLR_apriori      Vector of a priori LLR values per bit
   * \param[in]   method           Soft demodulation method
   * \return                       Vector of a posteriori LLR values
   */
  QLLRvec demodulate_soft_bits(const cvec &y, const cmat &H, double sigma2,
                               const QLLRvec &LLR_apriori,
                               Soft_Demod_Method method = FULL_ENUM_LOGMAP);


  /*!
   * \brief Soft MAP demodulation for parallelchannels without crosstalk.
   *
   * This function is a much faster equivalent to \c demodulate_soft_bits
   * with \f$H = \mbox{diag}(h)\f$. Its complexity is linear in the number
   * of subchannels.
   */
  void demodulate_soft_bits(const cvec &y, const cvec &h, double sigma2,
                            const QLLRvec &LLR_apriori,
                            QLLRvec &LLR_aposteriori);

  //! Print some properties of the MIMO modulator (mainly to aid debugging)
  friend ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const Modulator_NCD &m);

protected:
  //! Vectors of modulation symbols (along each dimension)
  Array<cvec> symbols;
	//! Complex-valued channel matrix
	itpp::cmat H;
	//! The spacing between different constellation points multiplied by the different H columns
	itpp::Array<itpp::Array<itpp::cvec> > hspacings;
	//! The spacing between different constellation points scaled by different y elements
  itpp::Array<itpp::vec> yspacings;
	void hxnormupdate(itpp::cvec& Hx, unsigned& bitstring, unsigned& ind, unsigned bit);
	void yxnormupdate(double& yx, itpp::QLLR& lapr, unsigned& bitstring, unsigned& ind, unsigned bit);
};

/*!
 * \relatesalso Modulator_NCD
 * \brief Print some properties of the MIMO modulator (mainly to aid debugging)
 */
ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const Modulator_NCD &m);


// ----------------------------------------------------------------------
// ND_UPAM
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * \brief Real-valued MIMO channel with uniform PAM along each dimension.
 *
 * <b>Example: (4 x 3 matrix channel with 4-PAM)</b>
 * \code
 * ND_UPAM chan; // multidimensional channel with uniform PAM
 * chan.set_M(3, 4); // 3-dimensional matrix channel, 4-PAM per dimension
 * cout << chan << endl;
 * bvec b = randb(3*2); // 3*2 bits in total
 * vec x = chan.modulate_bits(b);
 * mat H = randn(4,3); // 4 x 3 real matrix channel
 * double sigma2 = 0.01; // noise variance per real dimension
 * vec y = H*x + sqrt(sigma2)*randn(4); // transmit vector x
 * QLLRvec llr; // log-likelihood ratios
 * QLLRvec llr_ap = zeros_i(3*2);  // a priori equiprobable bits
 * chan.demodulate_soft_bits(y, H, sigma2, llr_ap, llr);
 * cout << "True bits:" << b << endl;
 * cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
 * \endcode
 *
 * <b>Example: (scalar channel with 8-PAM)</b>
 * \code
 * ND_UPAM chan;
 * chan.set_M(1, 8); // scalar channel, 8-PAM (3 bits per symbol)
 * cout << chan << endl;
 * bvec b = randb(3);
 * vec x = chan.modulate_bits(b);
 * mat H = "1.0";      // scalar channel
 * double sigma2 = 0.01;
 * vec y= H*x + sqrt(sigma2)*randn(); // transmit vector x
 * QLLRvec llr;
 * QLLRvec llr_ap = zeros_i(3);
 * chan.demodulate_soft_bits(y, H, sigma2, llr_ap, llr);
 * cout << "True bits:" << b << endl;
 * cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
 * \endcode
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT ND_UPAM : public Modulator_NRD
{
public:
  //! Constructor
  ND_UPAM(int nt = 1, int Mary = 2);
  //! Destructor
  virtual ~ND_UPAM() {}

  //! Set component modulators to M-PAM with Gray mapping
  void set_M(int nt = 1, int Mary = 2);

  //! Set component modulators to M-PAM with Gray mapping, different M per component
  void set_M(int nt = 1, ivec Mary = "2");

  /*!
   * \brief Sphere decoding
   *
   * This function solves the integer-constrained minimization problem
   * \f[
   * \mbox{min} |y - Hs|
   * \f]
   * with respect to \f$s\f$ using a sphere decoding algorithm and the
   * Schnorr-Eucner search strategy (see the source code for further
   * implementation notes). The function starts with an initial search
   * radius and increases it with a factor (\c stepup) until the search
   * succeeds.
   *
   * \param[in]  y              received data vector (\f$n_r\times 1\f$)
   * \param[in]  H              channel matrix (\f$n_r\times n_t\f$)
   * \param[in]  rmax           maximum possible sphere radius to try
   * \param[in]  rmin           sphere radius in the first try
   * \param[in]  stepup         factor with which the sphere radius is
   *                            increased if the search fails
   * \param[out] detected_bits  result of the search (hard decisions only,
   *                            QLLR for a sure "1" is set to 1000)
   * \return status of the decoding: 0 if the search suceeds, -1 otherwise
   */
  int sphere_decoding(const vec &y, const mat &H, double rmin, double rmax,
                      double stepup, QLLRvec &detected_bits);

private:
  // Sphere decoding search with Schnorr Eucner strategy.
  int sphere_search_SE(const vec &y, const mat &H, const imat &zrange,
                       double r, ivec &zhat);

  vec spacing;  // spacing between the constellation points

  inline int sign_nozero_i(int a) {
    return (a > 0 ? 1 : -1);
  }
  inline int sign_nozero_i(double a) {
    return (a > 0.0 ? 1 : -1);
  }
};

// ----------------------------------------------------------------------
// ND_UQAM
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * \brief Complex MIMO channel with uniform QAM per dimension
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT ND_UQAM : public Modulator_NCD
{
public:
  //! Constructor
  ND_UQAM(int nt = 1, int Mary = 4);
  //! Destructor
  virtual ~ND_UQAM() {}

  //! Set component modulators to M-QAM with Gray mapping
  void set_M(int nt = 1, int Mary = 4);

  //! Set component modulators to M-QAM with Gray mapping, different M per component
  void set_M(int nt = 1, ivec Mary = "4");

  /*!
   * \brief Set the constellation points
   * \param[in] nth - number of antenna
   * \param[in] inConstellation - new constellation points
   * \param[in] vector of mapping transmitted data symbols (index in the vector) to constellation points (value in the vector).
   *
   * <h3>Example of use:</h3>
   * \code
   * ND_UQAM qam(1,16); //QAM-16 with Gray mapping
   *
   * //now set the new Gray mapping (a little bit different than hardcoded in it++).
   * qam.set_constellation_points(1, qam.get_symbols()(0), "0 1 5 4 2 3 7 6 10 11 15 14 8 9 13 12");
   *
   * bvec bits = "0 1 1 0 0 0 1 1 1 1 0 0 1 0 0 1";
   * cvec modulated_symbols = qam.modulate_bits(bits);
   * \endcode
   */
  void set_constellation_points(const int nth, const cvec& inConstellation, const ivec& in_bit2symbols);

protected:
  ivec L;  //!< the square root of M
};

// ----------------------------------------------------------------------
// ND_UPSK
// ----------------------------------------------------------------------

/*!
 * \ingroup modulators
 * Complex MIMO channel with uniform PSK per dimension
 *
 * \note For issues relating to the accuracy of LLR computations,
 * please see the documentation of \c LLR_calc_unit
 */
class ITPP_EXPORT ND_UPSK : public Modulator_NCD
{
public:
  //! Constructor
  ND_UPSK(int nt = 1, int Mary = 4);
  //! Destructor
  virtual ~ND_UPSK() {}

  //! Set component modulators to M-QAM with Gray mapping
  void set_M(int nt = 1, int Mary = 4);

  //! Set component modulators to M-QAM with Gray mapping, different M per component
  void set_M(int nt = 1, ivec Mary = "4");
};


} // namespace itpp

#endif // #ifndef MODULATOR_ND_H

