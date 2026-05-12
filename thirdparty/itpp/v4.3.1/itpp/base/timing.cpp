/*!
 * \file
 * \brief Implementation of Timing classes
 * \author Thomas Eriksson, Tony Ottosson, Tobias Ringstrom and Adam Piatyszek
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#ifdef TIME_WITH_SYS_TIME
#  include <sys/time.h>
#  include <ctime>
#else
#  ifdef HAVE_SYS_TIME_H
#    include <sys/time.h>
#  else
#    include <ctime>
#  endif
#endif

#include <itpp/base/timing.h>
#include <cstdio>
#include <iostream>
#include <cmath>


#if defined(_WIN32) && !defined(__CYGWIN__)
#include <windows.h>

int gettimeofday(struct timeval* p, void*)
{
  union {
    long long ns100; /* time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } _now;

  GetSystemTimeAsFileTime(&(_now.ft));
  p->tv_usec = (long)((_now.ns100 / 10LL) % 1000000LL);
  /* time since 1 Jan 1970 */
  p->tv_sec = (long)((_now.ns100 - 116444736000000000LL) / 10000000LL);
  return 0;
}
#endif


namespace itpp
{

//! Global object for tic and toc functions
Real_Timer __tic_toc_timer;

//----------------------------------------------------------------------------
// class Timer
//----------------------------------------------------------------------------
Timer::Timer()
{
  reset();
}

void Timer::start(void)
{
  if (!running) {
    start_time = get_current_time();
    running = true;
  }
}

double Timer::stop(void)
{
  if (running) {
    stop_time = get_current_time();
    elapsed_time += stop_time - start_time;
    running = false;
  }

  return elapsed_time;
}

void Timer::reset(double t)
{
  elapsed_time = t;
  start_time = 0;
  stop_time = 0;
  running = false;
}

double Timer::get_time() const
{
  return running ?
         elapsed_time + get_current_time() - start_time :
         elapsed_time;
}

void Timer::tic(void)
{
  reset();
  start();
}

double Timer::toc(void)
{
  return get_time() ;
}

void Timer::toc_print(void)
{
  std::cout << "Elapsed time = " << get_time() << " seconds" << std::endl;
}

//----------------------------------------------------------------------------
// class CPU_Timer
//----------------------------------------------------------------------------
double CPU_Timer::get_current_time() const
{
  return static_cast<double>(clock()) / CLOCKS_PER_SEC;
}

//----------------------------------------------------------------------------
// class Real_Timer
//----------------------------------------------------------------------------
double Real_Timer::get_current_time() const
{
  struct timeval t;
  gettimeofday(&t, 0);
  return t.tv_sec + t.tv_usec * 1.0e-6;
}


void tic()
{
  __tic_toc_timer.tic();
}

double toc()
{
  return __tic_toc_timer.toc();
}

void toc_print()
{
  __tic_toc_timer.toc_print();
}

void pause(double t)
{
  if (t == -1) {
    std::cout << "(Press enter to continue)" << std::endl;
    getchar();
  }
  else {
    Real_Timer T;
    T.start();
    while (T.get_time() < t);
  }
}

} // namespace itpp
