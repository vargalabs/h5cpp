/*!
 * \file
 * \brief Implementation of Bit Error Rate Counter (BERC) and
 *        BLock Error Rate Counter (BLERC) classes
 * \author Pal Frenger and Adam Piatyszek
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

#include <itpp/comm/error_counters.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/converters.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>


namespace itpp
{

//-----------------------------------------------------------
// The Bit error rate counter class (BERC)
//-----------------------------------------------------------

BERC::BERC(int indelay, int inignorefirst, int inignorelast):
    delay(indelay), ignorefirst(inignorefirst), ignorelast(inignorelast),
    errors(0), corrects(0) {}

void BERC::count(const bvec &in1, const bvec &in2)
{
  int countlength = std::min(in1.length(), in2.length()) - std::abs(delay)
                    - ignorefirst - ignorelast;

  if (delay >= 0) {
    for (int i = 0; i < countlength; i++) {
      if (in1(i + ignorefirst) == in2(i + ignorefirst + delay)) {
        corrects++;
      }
      else {
        errors++;
      }
    }
  }
  else {
    for (int i = 0; i < countlength; i++) {
      if (in1(i + ignorefirst - delay) == in2(i + ignorefirst)) {
        corrects++;
      }
      else {
        errors++;
      }
    }
  }
}

void BERC::count(const bool x)
{
  if (x) {
    errors++; 
  } else {
    corrects++;
  }
}

void BERC::estimate_delay(const bvec &in1, const bvec &in2, int mindelay,
                          int maxdelay)
{
  int num, start1, start2;
  int min_input_length = std::min(in1.length(), in2.length());
  int bestdelay = mindelay;
  double correlation;
  double bestcorr = 0;
  for (int i = mindelay; i < maxdelay; i++) {
    num = min_input_length - std::abs(i) - ignorefirst - ignorelast;
    start1 = (i < 0) ? -i : 0;
    start2 = (i > 0) ?  i : 0;
    correlation = fabs(sum(to_vec(elem_mult(in1.mid(start1, num),
                                            in2.mid(start2, num)))));
    if (correlation > bestcorr) {
      bestdelay = i;
      bestcorr  = correlation;
    }
  }
  delay = bestdelay;
}

void BERC::report() const
{
  std::cout.setf(std::ios::fixed);
  std::cout << std::endl
            << "==================================" << std::endl
            << "     Bit Error Counter Report     " << std::endl
            << "==================================" << std::endl
            << " Ignore First           = " << ignorefirst << std::endl
            << " Ignore Last            = " << ignorelast << std::endl
            << " Delay                  = " << delay << std::endl
            << " Number of counted bits = " << std::setprecision(0)
            << (errors + corrects) << std::endl
            << " Number of errors       = " << std::setprecision(0)
            << errors << std::endl
            << "==================================" << std::endl
            << " Error rate             = " << std::setprecision(8)
            << (errors / (errors + corrects)) << std::endl
            << "==================================" << std::endl << std::endl;
}

double BERC::count_errors(const bvec &in1, const bvec &in2, int indelay,
                          int inignorefirst, int inignorelast)
{
  int countlength = std::min(in1.length(), in2.length()) - std::abs(indelay)
                    - inignorefirst - inignorelast;
  int local_errors = 0;

  if (indelay >= 0) {
    for (int i = 0; i < countlength; i++) {
      if (in1(i + inignorefirst) != in2(i + inignorefirst + indelay)) {
        local_errors++;
      }
    }
  }
  else {
    for (int i = 0; i < countlength; i++) {
      if (in1(i + inignorefirst - indelay) != in2(i + inignorefirst)) {
        local_errors++;
      }
    }
  }

  return local_errors;
}


//-----------------------------------------------------------
// The Block error rate counter class (BERC)
//-----------------------------------------------------------

BLERC::BLERC(): setup_done(false), blocksize(0), errors(0),
    corrects(0) {}


BLERC::BLERC(int inblocksize): setup_done(true), blocksize(inblocksize),
    errors(0), corrects(0) {}


void BLERC::set_blocksize(int inblocksize, bool clear)
{
  blocksize = inblocksize;
  if (clear) {
    errors = 0;
    corrects = 0;
  }
  setup_done = true;
}

void BLERC::count(const bvec &in1, const bvec &in2)
{
  it_assert(setup_done == true,
            "BLERC::count(): Block size has to be setup before counting errors.");
  int min_input_length = std::min(in1.length(), in2.length());
  it_assert(blocksize <= min_input_length,
            "BLERC::count(): Block size must not be longer than input vectors.");

  for (int i = 0; i < (min_input_length / blocksize); i++) {
    CORR = true;
    for (int j = 0; j < blocksize; j++) {
      if (in1(i * blocksize + j) != in2(i * blocksize + j)) {
        CORR = false;
        break;
      }
    }
    if (CORR) {
      corrects++;
    }
    else {
      errors++;
    }
  }
}

void BLERC::count(const bool x)
{
  if (x) {
    errors++; 
  } else {
    corrects++;
  }
}

} // namespace itpp
