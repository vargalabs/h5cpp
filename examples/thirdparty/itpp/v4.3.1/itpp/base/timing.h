/*!
 * \file
 * \brief Definitions of Timing classes
 * \author Thomas Eriksson, Tony Ottosson and Tobias Ringstrom
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

#ifndef TIMING_H
#define TIMING_H

#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup timers
*/

/*!
  \brief A virtual base class for timers.
  \ingroup timers
*/
class ITPP_EXPORT Timer
{
public:
  //! Create a new timer.  Sets the time to zero.
  Timer();
  //! Virtual destructor
  virtual ~Timer() { }
  //! Start the timer.  This does not set the time to zero.
  void start(void);
  //! Stop the timer.  Returns the elapsed time in seconds.
  double stop(void);
  //! Sets the time to time t, which is zero by default. Stops the timer if it is running.
  void reset(double t = 0.0);
  //! Resets the timer and starts it.
  void tic(void);
  //! Returns the elapsed time since last tic()
  double toc(void);
  //! Prints the elapsed time since last tic()
  void toc_print(void);
  //! Returns the elapsed time.
  double get_time() const;

protected:
  //! Vitrual function that returns teh current time
  virtual double get_current_time() const = 0;
  //! The start time of the timer
  double start_time;
  //! The stop time of the timer
  double stop_time;
  //! The ellapsed time from start to stop
  double elapsed_time;
  //! A bool that indicates if the timer is running or not
  bool running;
};

/*!
  \brief A CPU time timer class
  \ingroup timers

  Measures the time spent by the CPU on the current process.  If two processes
  are running concurrently, one real seconds equal 5 CPU seconds per process.
  The resolution is not very good (in the order of 0.01 seconds).

  Usage: Define a time object:
  \code
  CPU_Timer timer;
  \endcode

  Actions:
  Reset: \code timer.reset(); \endcode
  Start: \code timer.start(); \endcode
  Stop: \code timer.stop(); \endcode
  Get time: \code elapsedtime = timer.get_time(); \endcode

  It is possible to get elapsed time without stopping the timer.
  Observe that it is also possible to use the macros "time.tic();"
  to reset and start clock and "time.toc();" stop and print the elapsed time.
  \warning May give an negative answer if the measured time is too long.
*/
class ITPP_EXPORT CPU_Timer : public Timer
{
public:
  //! Create a new timer. Sets the time to zero.
  CPU_Timer() { }

protected:
  //!
  double get_current_time() const;
};

/*!
  \brief A real time timer class
  \ingroup timers

  Measures real time.

  Usage: Define a time object:
  \code
  Real_Timer timer;
  \endcode

  Actions:
  Reset: \code timer.reset(); \endcode
  Start: \code timer.start(); \endcode
  Stop: \code timer.stop(); \endcode
  Get time: \code elapsedtime = timer.get_time(); \endcode

  It is possible to get elapsed time without stopping the timer.
  Observe that it is also possible to use the macros "time.tic();"
  to reset and start clock and "time.toc_print();" to print the elapsed time.

  \warning May give an negative answer if the measured time is too long.
*/
class ITPP_EXPORT Real_Timer : public Timer
{
public:
  //! Create a new timer. Sets the time to zero.
  Real_Timer() { }

protected:
  //!
  double get_current_time() const;
};

/*!
  \brief Reset and start timer
  \ingroup timers
*/
ITPP_EXPORT void tic();

/*!
  \brief Returns the elapsed time since last tic()
  \ingroup timers
*/
ITPP_EXPORT double toc();

/*!
  \brief Prints the elapsed time since last tic()
  \ingroup timers
*/
ITPP_EXPORT void toc_print();

/*!
  \brief pause
  \ingroup timers

  \code pause(n); \endcode  Pauses for n seconds before continuing
  \code pause(); \endcode  Pauses until a key is pressed

*/
ITPP_EXPORT void pause(double t = -1);

} // namespace itpp

#endif // #ifndef TIMING_H
