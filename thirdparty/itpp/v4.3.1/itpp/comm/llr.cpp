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


#include <itpp/comm/llr.h>


namespace itpp
{

LLR_calc_unit::LLR_calc_unit()
{
  init_llr_tables();
}

LLR_calc_unit::LLR_calc_unit(short int d1, short int d2, short int d3)
{
  init_llr_tables(d1, d2, d3);
}

ivec LLR_calc_unit::get_Dint()
{
  ivec r(3);
  r(0) = Dint1;
  r(1) = Dint2;
  r(2) = Dint3;
  return r;
}

void LLR_calc_unit::init_llr_tables(short int d1, short int d2, short int d3)
{
  Dint1 = d1;      // 1<<Dint1 determines how integral LLRs relate to real LLRs (to_double=(1<<Dint)*int_llr)
  Dint2 = d2;      // number of entries in table for LLR operations
  Dint3 = d3;      // table resolution is 2^(-(Dint1-Dint3))
  logexp_table = construct_logexp_table();
}

ivec LLR_calc_unit::construct_logexp_table()
{
  ivec result(Dint2);
  for (int i = 0; i < Dint2; i++) {
    double x = pow2(static_cast<double>(Dint3 - Dint1)) * i;
    result(i) = to_qllr(std::log(1 + std::exp(-x)));
  }
  it_assert(length(result) == Dint2, "Ldpc_codec::construct_logexp_table()");

  return result;
}

QLLRvec LLR_calc_unit::to_qllr(const vec &l) const
{
  int n = length(l);
  ivec result(n);
  for (int i = 0; i < n; i++) {
    result.set(i, to_qllr(l(i)));
  }
  return result;
}

vec LLR_calc_unit::to_double(const QLLRvec &l) const
{
  int n = length(l);
  vec result(n);
  for (int i = 0; i < n; i++) {
    result.set(i, to_double(l(i)));
  }
  return result;
}

QLLRmat LLR_calc_unit::to_qllr(const mat &l)  const
{
  int m = l.rows();
  int n = l.cols();
  imat result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      result.set(i, j, to_qllr(l(i, j)));
    }
  }
  return result;
}

mat LLR_calc_unit::to_double(const QLLRmat &l) const
{
  int m = l.rows();
  int n = l.cols();
  mat result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      result.set(i, j, to_double(l(i, j)));
    }
  }
  return result;
}

// This function used to be inline, but in my experiments,
// the non-inlined version was actually faster /Martin Senst
QLLR LLR_calc_unit::Boxplus(QLLR a, QLLR b) const
{
  QLLR a_abs = (a > 0 ? a : -a);
  QLLR b_abs = (b > 0 ? b : -b);
  QLLR minabs = (a_abs > b_abs ? b_abs : a_abs);
  QLLR term1 = (a > 0 ? (b > 0 ? minabs : -minabs)
                    : (b > 0 ? -minabs : minabs));

  if (Dint2 == 0) {  // logmax approximation - avoid looking into empty table
    // Don't abort when overflowing, just saturate the QLLR
    if (term1 > QLLR_MAX) {
      it_info_debug("LLR_calc_unit::Boxplus(): LLR overflow");
      return QLLR_MAX;
    }
    if (term1 < -QLLR_MAX) {
      it_info_debug("LLR_calc_unit::Boxplus(): LLR overflow");
      return -QLLR_MAX;
    }
    return term1;
  }

  QLLR apb = a + b;
  QLLR term2 = logexp((apb > 0 ? apb : -apb));
  QLLR amb = a - b;
  QLLR term3 = logexp((amb > 0 ? amb : -amb));
  QLLR result = term1 + term2 - term3;

  // Don't abort when overflowing, just saturate the QLLR
  if (result > QLLR_MAX) {
    it_info_debug("LLR_calc_unit::Boxplus() LLR overflow");
    return QLLR_MAX;
  }
  if (result < -QLLR_MAX) {
    it_info_debug("LLR_calc_unit::Boxplus() LLR overflow");
    return -QLLR_MAX;
  }
  return result;
}

std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &lcu)
{
  os << "---------- LLR calculation unit -----------------" << std::endl;
  os << "LLR_calc_unit table properties:" << std::endl;
  os << "The granularity in the LLR representation is "
  << pow2(static_cast<double>(-lcu.Dint1)) << std::endl;
  os << "The LLR scale factor is " << (1 << lcu.Dint1) << std::endl;
  os << "The largest LLR that can be represented is "
  << lcu.to_double(QLLR_MAX) << std::endl;
  os << "The table resolution is "
  << pow2(static_cast<double>(lcu.Dint3 - lcu.Dint1)) << std::endl;
  os << "The number of entries in the table is " << lcu.Dint2 << std::endl;
  os << "The tables truncates at the LLR value "
  << pow2(static_cast<double>(lcu.Dint3 - lcu.Dint1)) * lcu.Dint2
  << std::endl;
  os << "-------------------------------------------------" << std::endl;
  return os;
}

}
