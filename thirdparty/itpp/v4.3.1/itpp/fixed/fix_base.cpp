/*!
 * \file
 * \brief Implementation of a base class for fixed-point data types
 * \author Johan Bergman
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

#include <itpp/fixed/fix_base.h>
#include <itpp/base/itassert.h>
#include <iostream>


namespace itpp
{

// Definition and initialization of static data member
output_mode Fix_Base::outputmode = OUTPUT_FIX_SHIFT;

void Fix_Base::set_output_mode(std::string o)
{
  if (o == "OUTPUT_FIX")
    outputmode = OUTPUT_FIX;
  else if (o == "OUTPUT_FIX_SHIFT")
    outputmode = OUTPUT_FIX_SHIFT;
  else if (o == "OUTPUT_FLOAT")
    outputmode = OUTPUT_FLOAT;
  else if (o == "OUTPUT_FLOAT_SHIFT")
    outputmode = OUTPUT_FLOAT_SHIFT;
  else
    it_error("Fix_Base::set_output_mode: Illegal output mode!");
}

void Fix_Base::print() const
{
  std::cout << "shift = " << shift << std::endl
            << "wordlen = " << wordlen << std::endl
            << "int(emode) = " << int(emode) << std::endl
            << "int(omode) = " << int(omode) << std::endl
            << "int(qmode) = " << int(qmode) << std::endl
            << "stat_ptr = " << stat_ptr << std::endl
            << "min = " << min << std::endl
            << "max = " << max << std::endl
            << "n_unused_bits = " << n_unused_bits << std::endl;
}

void Fix_Base::init()
{
  switch (emode) {
  case TC:
    it_assert_debug(wordlen >= 1 && wordlen <= 64, "Fix_Base::calc_apply_o_modes: Illegal word length!");
    max = fixrep(UINT64_POW2[wordlen - 1] - 1);
    min = -max - 1;
    break;
  case US:
    it_assert_debug(wordlen >= 0 && wordlen <= 63, "Fix_Base::calc_apply_o_modes: Illegal word length!");
    min = 0;
    max = fixrep(UINT64_POW2[wordlen] - 1);
    break;
  default:
    it_error("Fix_Base::init: Illegal sign encoding mode!");
    break;
  }

  n_unused_bits = MAX_WORDLEN - wordlen;
}

fixrep Fix_Base::apply_o_mode(fixrep x) const
{
  fixrep ret = x;
  bool overflow = false;

  if (ret < min) {
    overflow = true;
    switch (omode) {
    case WRAP:
      ret = fixrep((fixrep(ret) << n_unused_bits) >> n_unused_bits);
      break;
    case SAT:
      ret = min;
      break;
    default:
      it_error("Fix_Base::apply_o_mode: Illegal overflow mode!");
      break;
    }
  }
  else if (ret > max) {
    overflow = true;
    switch (omode) {
    case WRAP:
      ret = fixrep((fixrep(ret) << n_unused_bits) >> n_unused_bits);
      break;
    case SAT:
      ret = max;
      break;
    default:
      it_error("Fix_Base::apply_o_mode: Illegal overflow mode!");
      break;
    }
  }

  if (stat_ptr != 0)
    stat_ptr->sample(double(ret), overflow);

  return ret;
}

fixrep Fix_Base::scale_and_apply_modes(double x, q_mode q) const
{
  it_assert_debug(shift >= -64 && shift <= 63, "Fix_Base::scale_and_apply_modes: Illegal shift!");
  fixrep ret = 0;
  double scaled_value = x * DOUBLE_POW2[shift + 64];

  switch (q) {
  case RND:
    ret = apply_o_mode(fixrep(std::floor(scaled_value + 0.5)));
    break;
  case RND_ZERO:
    if (x < 0)
      ret = apply_o_mode(fixrep(std::floor(scaled_value + 0.5)));
    else
      ret = apply_o_mode(fixrep(-std::floor(-scaled_value + 0.5)));
    break;
  case RND_MIN_INF:
    ret = apply_o_mode(fixrep(-std::floor(-scaled_value + 0.5)));
    break;
  case RND_INF:
    if (x < 0)
      ret = apply_o_mode(fixrep(scaled_value - 0.5));
    else
      ret = apply_o_mode(fixrep(scaled_value + 0.5));
    break;
  case RND_CONV:
    if (scaled_value == std::floor(scaled_value) + 0.5)
      ret = apply_o_mode((fixrep(round(scaled_value)) >> 1) << 1);
    else
      ret = apply_o_mode(fixrep(std::floor(scaled_value + 0.5)));
    break;
  case RND_CONV_ODD:
    if (scaled_value == std::floor(scaled_value) + 0.5)
      if (scaled_value < 0)
        ret = apply_o_mode(((fixrep(std::ceil(scaled_value)) >> 1) << 1) - 1);
      else
        ret = apply_o_mode(((fixrep(std::floor(scaled_value)) >> 1) << 1) + 1);
    else
      ret = apply_o_mode(fixrep(std::floor(scaled_value + 0.5)));
    break;
  case TRN:
    ret = apply_o_mode(fixrep(std::floor(scaled_value)));
    break;
  case TRN_ZERO:
    ret = apply_o_mode(fixrep(scaled_value));
    break;
  default:
    it_error("Fix_Base::scale_and_apply_modes: Illegal quantization mode!");
    break;
  }

  return ret;
}

fixrep Fix_Base::rshift_and_apply_q_mode(fixrep x, int n, q_mode q) const
{
  it_assert_debug(n >= 0, "Fix_Base::rshift_and_apply_q_mode: n cannot be negative!");
  fixrep ret = 0;

  if (n == 0) {
    ret = x;
  }
  else {
    switch (q) {
    case RND:
      // Add the most significant deleted bit to the remaining bits
      ret = ((x >> (n - 1)) + 1) >> 1;
      break;
    case RND_ZERO:
      // If the most significant deleted bit is 1,
      // and either the sign bit or at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x & (fixrep(1) << (n - 1))) && ((x < 0) || (x & ((fixrep(1) << (n - 1)) - 1))))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    case RND_MIN_INF:
      // If the most significant deleted bit is 1,
      // and at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x & (fixrep(1) << (n - 1))) && (x & ((fixrep(1) << (n - 1)) - 1)))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    case RND_INF:
      // If the most significant deleted bit is 1,
      // and either the inverted value of the sign bit or at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x & (fixrep(1) << (n - 1))) && ((x >= 0) || (x & ((fixrep(1) << (n - 1)) - 1))))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    case RND_CONV:
      // If the most significant deleted bit is 1,
      // and either the least significant of the remaining bits or at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x & (fixrep(1) << (n - 1))) && ((x & (fixrep(1) << n)) || (x & ((fixrep(1) << (n - 1)) - 1))))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    case RND_CONV_ODD:
      // If the most significant deleted bit is 1,
      // and either the least significant of the remaining bits is 0 or at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x & (fixrep(1) << (n - 1))) && (!(x & (fixrep(1) << n)) || (x & ((fixrep(1) << (n - 1)) - 1))))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    case TRN:
      // Just copy the remaining bits
      ret = x >> n;
      break;
    case TRN_ZERO:
      // If the sign bit is 1,
      // and either the most significant deleted bit or at least one other deleted bit is 1,
      // add 1 to the remaining bits
      if ((x < 0) && (x & ((fixrep(1) << n) - 1)))
        ret = (x >> n) + 1;
      else
        ret = x >> n;
      break;
    default:
      it_error("Fix_Base::rshift_and_apply_q_mode: Illegal quantization mode!");
      break;
    }
  }

  if (stat_ptr != 0)
    stat_ptr->sample(double(ret), false);

  return ret;
}

} // namespace itpp
