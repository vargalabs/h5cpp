/*!
 * \file
 * \brief Implementation of spread spectrum classes and functions
 * \author Tony Ottosson
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

#include <itpp/comm/spread.h>
#include <itpp/base/converters.h>
#include <itpp/stat/misc_stat.h>


namespace itpp
{

//------------- Spread1d -------------------

Spread_1d::Spread_1d(const vec &incode)
{
  set_code(incode);
}

void Spread_1d::set_code(const vec &incode)
{
  N = incode.size();
  code = incode;
  code /= norm(code);
}

vec Spread_1d::get_code()
{
  return code;
}

void Spread_1d::spread(const vec &symbols, vec &out)
{
  out.set_size(symbols.length()*N, false);

  for (int i = 0;i < symbols.length();i++)
    out.replace_mid(i*N, symbols(i)*code);
}

void Spread_1d::despread(const vec &rec_signal, vec &out, int timing)
{
  int nosymbols = (int)std::floor(double((rec_signal.length() - timing)) / N);
  out.set_size(nosymbols);

  for (int i = 0;i < nosymbols;i++)
    out(i) = rec_signal.mid(i * N + timing, N) * code;
}


//---------------- Spread2d ----------------------

vec Spread_2d::get_codeI()
{
  return spreadI.get_code();
}

vec Spread_2d::get_codeQ()
{
  return spreadQ.get_code();
}

Spread_2d::Spread_2d(const vec &incodeI, const vec &incodeQ)
{
  set_code(incodeI, incodeQ);
}

void Spread_2d::set_code(const vec &incodeI, const vec &incodeQ)
{
  it_assert(incodeI.length() == incodeQ.length(), "Size of I and Q codes doesn't match");
  spreadI.set_code(incodeI);
  spreadQ.set_code(incodeQ);
}

void Spread_2d::spread(const cvec &symbols, cvec &out)
{
  out = to_cvec(spreadI.spread(real(symbols)), spreadQ.spread(imag(symbols)));
}

void Spread_2d::despread(const cvec &rec_signal, cvec &out, int timing)
{
  out = to_cvec(spreadI.despread(real(rec_signal), timing), spreadQ.despread(imag(rec_signal), timing));
}



//------------- Multicode_Spread_1d ----------------


Multicode_Spread_1d::Multicode_Spread_1d(const mat &incodes)
{
  set_codes(incodes);
}

void Multicode_Spread_1d::set_codes(const mat &incodes)
{
  codes = incodes;
  N = incodes.cols();
  L = incodes.rows();
  for (int i = 0; i < L; i++)
    codes.set_row(i, codes.get_row(i) / norm(codes.get_row(i)));
}

mat Multicode_Spread_1d::get_codes()
{
  return codes;
}

vec Multicode_Spread_1d::spread(const vec &symbols)
{
  int i;
  int nomcsymbols = (int)std::floor(double(symbols.length() / L));
  vec temp(nomcsymbols*N);

  for (i = 0;i < nomcsymbols;i++) {
    temp.replace_mid(i*N, codes.T() * symbols.mid(i*L, L)); // TODO: this is now very slow
  }

  return temp;
}

vec Multicode_Spread_1d::despread(const vec &receivedsignal, int timing)
{
  int i;
  int nosymbols = (int)std::floor(double((receivedsignal.length() - timing)) / N);
  vec temp(nosymbols*L);

  for (i = 0;i < nosymbols;i++) {
    temp.replace_mid(i*L, codes*receivedsignal.mid(i*N + timing, N));
  }
  return temp;
}


//----------------- Multicode_Spread_2d -------------------


Multicode_Spread_2d::Multicode_Spread_2d(const mat &incodesI, const mat &incodesQ)
{
  set_codes(incodesI, incodesQ);
}

mat Multicode_Spread_2d::get_codesI()
{
  return mcspreadI.get_codes();
}

mat Multicode_Spread_2d::get_codesQ()
{
  return mcspreadQ.get_codes();
}

void Multicode_Spread_2d::set_codes(const mat &incodesI, const mat &incodesQ)
{
  it_assert(incodesI.rows() == incodesQ.rows() && incodesI.cols() == incodesQ.cols(),
            "Multicode_Spread_2d::set_codes(): dimension mismatch");
  mcspreadI.set_codes(incodesI);
  mcspreadQ.set_codes(incodesQ);
}

cvec Multicode_Spread_2d::spread(const cvec &symbols)
{
  return to_cvec(mcspreadI.spread(real(symbols)), mcspreadQ.spread(imag(symbols)));
}

cvec Multicode_Spread_2d::despread(const cvec &receivedsignal, int timing)
{
  return to_cvec(mcspreadI.despread(real(receivedsignal), timing), mcspreadQ.despread(imag(receivedsignal), timing));
}

} // namespace itpp
