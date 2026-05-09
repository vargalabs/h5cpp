/*!
 * \file
 * \brief Deterministic sources test program
 * \author Tobias Ringstrom, Tony Ottosson and Adam Piatyszek
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

#include <itpp/itsignal.h>
#include <itpp/itstat.h>
#include <iomanip>

using namespace itpp;
using namespace std;

#define LOOP_SIZE 100000
#define THRESHOLD 1e-13

#define REALRUN(name,s)     \
  for (int i = 0; i < LOOP_SIZE; i++)   \
    real_result(i) = s();    \
  show(name, mean(real_result), variance(real_result));


void show(const char *name, double sm, double sv)
{
  cout << setw(18) << name << "  "
       << setw(20) << round_to_zero(sm, THRESHOLD) << "  "
       << setw(20) << round_to_zero(sv, THRESHOLD) << endl;
}


int main()
{
  Sine_Source      s10(20.0 / LOOP_SIZE);
  Square_Source    s11(20.0 / LOOP_SIZE);
  Triangle_Source  s12(20.0 / LOOP_SIZE);
  Sawtooth_Source  s13(20.0 / LOOP_SIZE);
  Impulse_Source   s14(20.0 / LOOP_SIZE);
  Pattern_Source   s15(vec("1 3"));

  RNG_reset(12345);

  cout.setf(ios::fixed);
  cout.precision(8);
  cout << setw(18) << "Source" << "  "
       << setw(20) << "sim mean" << "  "
       << setw(20) << "sim var" << endl
       << "============================================================================" << endl;

  vec real_result(LOOP_SIZE);

  REALRUN("Sine",        s10);
  REALRUN("Square",      s11);
  REALRUN("Triangle",    s12);
  REALRUN("Sawtooth",    s13);
  REALRUN("Impulse",     s14);
  REALRUN("Pattern",     s15);

  return 0;
}
