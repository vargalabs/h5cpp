/*!
 * \file
 * \brief Error handling functions - header file
 * \author Tobias Ringstrom and Adam Piatyszek
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

#ifndef ITASSERT_H
#define ITASSERT_H

#include <sstream>
#include <string>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup errorhandlingfunc

  For the following macros, the argument \c s is a string that is displayed.

  \code
  it_assert(t,s);           // Abort if t is not true
  it_assert_debug(t,s);     // Abort if t is not true and NDEBUG is not defined
  it_error_if(t,s);         // Abort if t is true
  it_error(s);              // Abort
  it_info(s);               // Display a message
  it_info_debug(s);         // Display a message if NDEBUG is not defined
  it_info_no_endl(s);       // Display a message without appended "std::endl"
  it_info_no_endl_debug(s); // Display a message without appended "std::endl" if NDEBUG is not defined
  it_warning(s);            // Display a warning
  \endcode

  \c it_assert(), \c it_error(), \c it_error_if(), \c it_info(), \c
  it_info_no_endl() and \c it_warning() are always active, whereas \c
  it_assert_debug(), \c it_info_debug() and \c it_info_no_endl_debug()
  depends on the \c NDEBUG compile time definition. If \c NDEBUG is
  defined, then none of these macros is executed.

  \note \c it_assert0() and \c it_assert1() macros are still defined for
  backward compatibility, but \c it_assert_debug() should be used instead
  of them.
*/
//!@{

//! Helper function for the \c it_assert and \c it_assert_debug macros
ITPP_EXPORT void it_assert_f(std::string ass, std::string msg, std::string file, int line);
//! Helper function for the \c it_error and \c it_error_if macros
ITPP_EXPORT void it_error_f(std::string msg, std::string file, int line);
//! Helper function for the \c it_info and \c it_info_debug macros
ITPP_EXPORT void it_info_f(std::string msg);
//! Helper function for the \c it_warning macro
ITPP_EXPORT void it_warning_f(std::string msg, std::string file, int line);

//! Enable/disable using exceptions for error handling.
ITPP_EXPORT void it_enable_exceptions(bool on);
//! Enable warnings
ITPP_EXPORT void it_enable_warnings();
//! Disable warnings
ITPP_EXPORT void it_disable_warnings();
//! Redirect warnings to the ostream warn_stream
ITPP_EXPORT void it_redirect_warnings(std::ostream *warn_stream);

//! Style of assert, error and warning messages.
enum error_msg_style { Full, Minimum };

//! Set preferred style of assert, error and warning messages
ITPP_EXPORT void it_error_msg_style(error_msg_style style);


//! Abort if \c t is not true
#define it_assert(t,s)      \
  if (!(t)) {       \
    std::ostringstream m_sout;     \
    m_sout << s;      \
    itpp::it_assert_f(#t,m_sout.str(),__FILE__,__LINE__); \
  } else       \
    ((void) 0)

#if defined(NDEBUG)
//! Abort if \c t is not true and NDEBUG is not defined
#  define it_assert_debug(t,s) ((void) (t))
#else
//! Abort if \c t is not true and NDEBUG is not defined
#  define it_assert_debug(t,s) it_assert(t,s)
#endif // if defined(NDEBUG)

//! Deprecated macro. Please use \c it_assert_debug() instead.
#define it_assert0(t,s) it_assert_debug(t,s)
//! Deprecated macro. Please use \c it_assert_debug() instead.
#define it_assert1(t,s) it_assert_debug(t,s)


//! Abort if \c t is true
#define it_error_if(t,s)    \
  if((t)) {      \
    std::ostringstream m_sout;    \
    m_sout << s;     \
    itpp::it_error_f(m_sout.str(),__FILE__,__LINE__); \
  } else      \
    ((void) 0)

//! Abort unconditionally
#define it_error(s)     \
  if (true) {      \
    std::ostringstream m_sout;    \
    m_sout << s;     \
    itpp::it_error_f(m_sout.str(),__FILE__,__LINE__); \
  } else      \
    ((void) 0)


//! Print information message
#define it_info(s)    \
  if (true) {     \
    std::ostringstream m_sout;   \
    m_sout << s << std::endl;   \
    itpp::it_info_f(m_sout.str());  \
  } else     \
    ((void) 0)

//! Print information message withot \c std::endl at the end
#define it_info_no_endl(s)   \
  if (true) {     \
    std::ostringstream m_sout;   \
    m_sout << s;    \
    itpp::it_info_f(m_sout.str());  \
  } else     \
    ((void) 0)

#if defined(NDEBUG)
//! Print information message if NDEBUG is not defined
#  define it_info_debug(s) ((void) 0)
/*!
  \brief Print information message without \c std::endl at the end if
  NDEBUG is not defined
*/
#  define it_info_no_endl_debug(s) ((void) 0)
#else
//! Print information message if NDEBUG is not defined
#  define it_info_debug(s) it_info(s)
/*!
  \brief Print information message withot \c std::endl at the end if
  NDEBUG is not defined
*/
#  define it_info_no_endl_debug(s) it_info_no_endl(s)
#endif // if defined(NDEBUG)


//! Display a warning message
#define it_warning(s)     \
  if (true) {      \
    std::ostringstream m_sout;    \
    m_sout << s;     \
    itpp::it_warning_f(m_sout.str(),__FILE__,__LINE__); \
  } else      \
    ((void) 0)

//!@}

} // namespace itpp

#endif // #ifndef ITASSERT_H
