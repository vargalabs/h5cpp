/*!
 * \file
 * \brief Include file for exporting/importing IT++ library symbols
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#ifndef ITEXPORTS_H
#define ITEXPORTS_H

/*defined from cmake when using the shared version of IT++ library*/
#cmakedefine ITPP_SHARED_LIB

/*needed to export shared library symbols on Windows*/
#if defined(ITPP_SHARED_LIB) && defined(_MSC_VER)
  #ifndef ITPP_EXPORT
    #if defined(itpp_EXPORTS) || defined(itpp_debug_EXPORTS) /*automatically defined by cmake*/
      #define ITPP_EXPORT __declspec(dllexport)
    #else
      #define ITPP_EXPORT __declspec(dllimport)
    #endif
  #endif
#endif

#if (__GNUC__ >= 4) /*UNIX*/
  #ifndef ITPP_EXPORT_TEMPLATE
    #define ITPP_EXPORT_TEMPLATE extern
  #endif
#endif

#ifndef ITPP_EXPORT
  #define ITPP_EXPORT
#endif
#ifndef ITPP_EXPORT_TEMPLATE
  #define ITPP_EXPORT_TEMPLATE
#endif

#endif
