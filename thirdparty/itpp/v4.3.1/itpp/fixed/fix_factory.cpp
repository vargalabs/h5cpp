/*!
 * \file
 * \brief Implementation of a class factory for fixed-point data types Fix
 *        and CFix
 * \author Johan Bergman and Adam Piatyszek
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

#include <itpp/fixed/fix_factory.h>
#include <itpp/fixed/cfix.h>


namespace itpp
{

void Fix_Factory::create(Fix* &ptr, const int n) const
{
  void *p = ::operator new(sizeof(Fix) * n);
  ptr = reinterpret_cast<Fix*>(p);
  // Set fixed-point restrictions
  for (int i = 0; i < n; ++i) {
    new(ptr + i) Fix(0.0, 0, wordlen, emode, omode, qmode, stat_ptr);
  }
}

void Fix_Factory::create(CFix* &ptr, const int n) const
{
  void *p = ::operator new(sizeof(CFix) * n);
  ptr = reinterpret_cast<CFix*>(p);
  // Set fixed-point restrictions
  for (int i = 0; i < n; ++i) {
    new(ptr + i) CFix(0.0, 0, wordlen, emode, omode, qmode, stat_ptr);
  }
}

template<>
void create_elements<Fix>(Fix* &ptr, const int n, const Factory &f)
{
  if (const Fix_Factory *fix_factory_ptr = dynamic_cast<const Fix_Factory*>(&f)) {
    // Yes, f seems to be a Fix_Factory. Now call the Fix_Factory::create method
    fix_factory_ptr->create(ptr, n);
  }
  else {
    // No, f does not seem to be a Fix_Factory. As a fallback solution,
    // assume that f is DEFAULT_FACTORY and use the default constructor
    void *p = ::operator new(sizeof(Fix) * n);
    ptr = reinterpret_cast<Fix*>(p);
    for (int i = 0; i < n; i++) {
      new(ptr + i) Fix();
    }
  }
}

template<>
void create_elements<CFix>(CFix* &ptr, const int n, const Factory &f)
{
  if (const Fix_Factory *fix_factory_ptr = dynamic_cast<const Fix_Factory*>(&f)) {
    // Yes, f seems to be a Fix_Factory. Now call the Fix_Factory::create method
    fix_factory_ptr->create(ptr, n);
  }
  else {
    // No, f does not seem to be a Fix_Factory. As a fallback solution,
    // assume that f is DEFAULT_FACTORY and use the default constructor
    void *p = ::operator new(sizeof(CFix) * n);
    ptr = reinterpret_cast<CFix*>(p);
    for (int i = 0; i < n; i++) {
      new(ptr + i) CFix();
    }
  }
}

} // namespace itpp
