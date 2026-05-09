/*!
 * \file
 * \brief Class for numerically efficient log-likelihood algebra
 * \author Erik G. Larsson and Martin Senst
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

#ifndef LLR_H
#define LLR_H

#include <limits>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>
#include <limits>
#include <itpp/itexports.h>

namespace itpp
{

/*! \relates LLR_calc_unit
  The quantized log-likelihood ratio (QLLR) representation, scalar form. See \c LLR_calc_unit.
*/
typedef signed int QLLR;

/*!  \relates LLR_calc_unit
  The quantized log-likelihood ratio (QLLR) representation, vector form. See \c LLR_calc_unit.
*/
typedef Vec<QLLR> QLLRvec;

/*!  \relates LLR_calc_unit
  The quantized log-likelihood ratio (QLLR) representation, matrix form. See \c LLR_calc_unit.
*/
typedef Mat<QLLR> QLLRmat;

/*!  \relates LLR_calc_unit
  The largest possible QLLR value
*/
const QLLR QLLR_MAX = (std::numeric_limits<QLLR>::max() >> 4);
// added some margin to make sure the sum of two LLR is still permissible

/*!
  \brief Log-likelihood algebra calculation unit.

  This class contains functions for algebra with log-likelihood ratios
  (LLRs). The class is mainly useful for modules that rely on certain
  nonlinear operations on LLRs (e.g. boxplus for LDPC and turbo
  decoding and Jacobian logarithm for soft demodulation). The routines
  provided are numerically efficient and use only integer
  arithmetics. Additionally, they allow for arbitrary quantization of
  LLRs in order to study effects of limited precision (fixed point
  number representations). As a consequence, computations are
  approximate. With the standard settings, the numerical precision
  should be sufficient (and practically equivalent to floating point
  precision) for all practical applications.  However, one can
  construct cases where small numerical artifacts due to quantization
  effects (e.g., soft demodulation with very high or very low SNR) can
  be observed.

  An LLR for an information bit b is defined according to \f[ L = \log
  \frac{P(b=0)}{P(b=1)} \f] and it is in general a real number.
  LLR values are represented via the special type, "quantized
  LLR" (QLLR).  The relation between the quantized representation
  and the real (floating-point) LLR value is
  \f[ \mbox{QLLR} = \mbox{round} \left(2^{\mbox{Dint1}}\cdot
  \mbox{LLR}\right) \f]
  The user parameter Dint1 determines the
  granularity of the quantization, and it can be set arbitrarily.
  The functions to_double() and to_qllr() can be used to perform
  conversions between the two representations (QLLR to
  floating-point, and vice versa).

  The class provides functions for the computation of the Jacobian
  logarithm and Hagenauer's "boxplus" operator.  These functions are
  based on a table-lookup.  The resolution of the table is
  determined by the parameters Dint2 and Dint3. See the class
  constructor for more detail.  When an object of LLR_calc_unit is
  created, corresponding lookup-tables are also generated. The
  resolution of these tables can be adjusted by providing parameters
  to the constructor.

  The variable table resolution allows one to study complexity
  versus accuracy (i.e., how different table resolutions would
  degrade performance) to some extent. Yet the main purpose of the
  QLLR representation is to provide a tool for writing efficient
  simulation code, rather than to provide for bit-level
  (fixed-point) simulations. For bit-level simulations, a true fixed
  point representation of LLRs would be preferred/required.  With
  the default setting of the table parameters, using the QLLR type
  is practically as accurate (but much faster) as using "double" to
  represent LLRs.  Decoder implementations may then provide
  functions using QLLR, fixed-point, or double (for compatibility
  reasons) representations of LLR values.

  Note: the QLLR type does not check that the correct quantization
  level is used. I.e., in theory it would be possible to add two
  QLLR types with different quantization (Dint) parameters.  This is
  intentionally implemented this way to achieve maximum runtime
  efficiency.

*/
class ITPP_EXPORT LLR_calc_unit
{
public:
  //! Constructor, using the default table resolution
  LLR_calc_unit();

  /*!
   * \brief Constructor, using a specific table resolution.
   *
   * See init_llr_tables() for more details on the parameters.
   */
  LLR_calc_unit(short int Dint1, short int Dint2, short int Dint3);

  /*! \brief Set the quantization and table parameters

    \param Dint1 Determines the relation between LLR represented as
    real number and as integer.  The relation is
    \f[ \mbox{QLLR} = \mbox{round} \left(2^{\mbox{Dint1}}\cdot
    \mbox{LLR}\right) \f]

    \param Dint2 Number of entries in the table. If this is zero,
    then logmap becomes logmax.

    \param Dint3 Determines the table resolution. The spacing between each
    entry is \f[ 2^{-(Dint1-Dint3)} \f]

    The default parameter values are chosen to give a performance
    practically indistinguishable from that of using floating point
    calculations.

    Example: (recommended settings with "exact" computation via high
    resolution lookup table)
    \code
    LLR_calc_unit lcalc(12, 300, 7);
    \endcode

    Example: (recommended settings with logmax, i.e. no table lookup)
    \code
    LLR_calc_unit lcalc(12, 0, 7);
    \endcode
  */
  void init_llr_tables(short int Dint1 = 12, short int Dint2 = 300,
                       short int Dint3 = 7);

  //! Convert a "real" LLR value to an LLR type
  QLLR to_qllr(double l) const;

  //! Convert a vector of "real" LLR values to an LLR type
  QLLRvec to_qllr(const vec &l) const;

  //! Convert a matrix of "real" LLR values to an LLR type
  QLLRmat to_qllr(const mat &l) const;

  //! Convert an LLR type to a "real" LLR
  double to_double(QLLR l) const;

  //! Convert a vector of LLR types to a "real" LLR
  vec to_double(const QLLRvec &l) const;

  //! Convert a matrix of LLR types to a "real" LLR
  mat to_double(const QLLRmat &l) const;

  /*!
   * \brief Jacobian logarithm.
   *
   * This function computes \f[ \log(\exp(a)+\exp(b)) \f]
   */
  inline QLLR jaclog(QLLR a, QLLR b) const;
  // Note: a version of this function taking "double" values as input
  // is deliberately omitted, because this is rather slow.

  /*!
   * \brief Hagenauer's "Boxplus" operator.
   *
   * This function computes:
   * \f[ \mbox{sign}(a) * \mbox{sign}(b) * \mbox{min}(|a|,|b|)
   * + f(|a+b|) - f(|a-b|) \f]
   * where \f[ f(x) = \log(1+\exp(-x))  \f]
   */
  QLLR Boxplus(QLLR a, QLLR b) const;

  /*!
   * \brief Logexp operator.
   *
   * This function computes \f[ f(x) = \log(1+\exp(-x)) \f]
   */
  inline QLLR logexp(QLLR x) const;

  //! Retrieve the table resolution values
  ivec get_Dint();

  //! Print some properties of the LLR calculation unit in plain text
  friend ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &l);

private:
  //! Compute the table for \f[ f(x) = \log(1+\exp(-x)) \f]
  ivec construct_logexp_table();

  //! The lookup tables for the decoder
  ivec logexp_table;

  //! Decoder (lookup-table) parameters
  short int Dint1, Dint2, Dint3;
};

/*!
  \relatesalso LLR_calc_unit
  \brief Print some properties of the LLR calculation unit in plain text.
*/
ITPP_EXPORT std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &lcu);


// ----------------------------------------------------------------------
// implementation of some inline functions
// ----------------------------------------------------------------------

inline double LLR_calc_unit::to_double(QLLR l) const
{
  return static_cast<double>(l) / (1 << Dint1);
}

inline QLLR LLR_calc_unit::to_qllr(double l) const
{
  double QLLR_MAX_double = to_double(QLLR_MAX);
  // Don't abort when overflow occurs, just saturate the QLLR
  if (l > QLLR_MAX_double) {
    it_info_debug("LLR_calc_unit::to_qllr(): LLR overflow");
    return QLLR_MAX;
  }
  if (l < -QLLR_MAX_double) {
    it_info_debug("LLR_calc_unit::to_qllr(): LLR overflow");
    return -QLLR_MAX;
  }
  return static_cast<QLLR>(std::floor(0.5 + (1 << Dint1) * l));
}


inline QLLR LLR_calc_unit::logexp(QLLR x) const
{
  it_assert_debug(x >= 0, "LLR_calc_unit::logexp(): Wrong LLR value");
  int ind = x >> Dint3;
  if (ind >= Dint2) // outside table
    return 0;

  it_assert_debug(ind >= 0, "LLR_calc_unit::logexp(): Internal error");
  it_assert_debug(ind < Dint2, "LLR_calc_unit::logexp(): internal error");

  // With interpolation
  // int delta=x-(ind<<Dint3);
  // return ((delta*logexp_table(ind+1) + ((1<<Dint3)-delta)*logexp_table(ind)) >> Dint3);

  // Without interpolation
  return logexp_table(ind);
}


inline QLLR LLR_calc_unit::jaclog(QLLR a, QLLR b) const
{
  QLLR x, maxab;

  if (a > b) {
    maxab = a;
    x = a - b;
  }
  else {
    maxab = b;
    x = b - a;
  }

  if (maxab >= QLLR_MAX)
    return QLLR_MAX;
  else
    return (maxab + logexp(x));
}

}

#endif
