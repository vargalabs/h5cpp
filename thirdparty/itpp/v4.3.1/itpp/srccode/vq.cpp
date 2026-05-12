/*!
 * \file
 * \brief Implementation of a vector quantizer class (unconstrained)
 * \author Thomas Eriksson
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

#include <itpp/srccode/vq.h>
#include <itpp/base/array.h>
#include <itpp/base/matfunc.h>
#include <fstream>
#include <iostream>
#include <cstdlib>

//! \cond

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

namespace itpp
{

//--------------------------------------------------------------------
//    class VQ
//--------------------------------------------------------------------

Vector_Quantizer::Vector_Quantizer() : CodeBook()
{
  LatestDist = 0;
  Size = 0;
  Dim = 0;
}

Vector_Quantizer::Vector_Quantizer(const char *Name) : CodeBook()
{
  LatestDist = 0;
  Size = 0;
  Dim = 0;
  load(Name);
}


int Vector_Quantizer::encode(const vec &x)
{
  int i;
  double S, MinS = 1.0E30F;
  int MinIndex = 0;
  int j, pos = 0;
  double a;

  for (i = 0;i < Size;i++) {
    S = 0;
    for (j = 0;j < Dim;j++) {
      a = x._elem(j) - CodeBook._elem(pos + j);
      S += a * a;
      if (S >= MinS) goto sune;
    }
    MinS = S;
    MinIndex = i;
  sune:
    pos += Dim;
  }
  LatestDist = MinS;
  return MinIndex;
}

ivec Vector_Quantizer::encode(const vec &x, int num)
{
  double S, a;
  vec  MinS(num);
  ivec MinIndex(num);
  int i, j, index, pos = 0;

  MinS.clear();
  MinS += 1.0E30F;
  MinIndex.clear();
  for (i = 0;i < Size;i++) {
    S = 0;
    for (j = 0;j < Dim;j++) {
      a = x._elem(j) - CodeBook._elem(pos + j);
      S += a * a;
      if (S >= MinS[num-1]) goto sune;
    }
    for (index = num - 2;(index >= 0) && (S < MinS[index]);index--);
    for (j = MinS.length() - 2;j > index;j--) {
      MinS[j+1] = MinS[j];// memcpy, memmov
      MinIndex[j+1] = MinIndex[j];
    }
    MinS[index+1] = S;
    MinIndex[index+1] = i;
  sune:
    pos += Dim;
  }
  LatestDist = MinS[0];
  return MinIndex;
}

Array<vec> Vector_Quantizer::decode(const ivec &Index) const
{
  Array<vec> Temp(Index.length());

  for (int i = 0;i < Temp.length();i++) {
    Temp(i) = get_codevector(Index(i));
  }
  return Temp;
}


ifstream &operator>>(ifstream &ifs, vec &v)
{
  int    i;
  char    str[2000];
  char    *ptr, *ptr_old;
  bool flag;
  if (length(v) != 0) {
    for (i = 0;i < length(v);i++) {
      ifs.operator >> (v[i]) ;
    }
  }
  else {
    v.set_length(50);
    ifs.getline(str, 2000);
    if (strlen(str) == 0) ifs.getline(str, 2000);
    i = 0;
    v[i++] = atof(str);
    ptr = str;
    ptr_old = ptr;
    ptr = strchr(ptr, ' ');
    while (ptr == ptr_old) {
      ptr++;
      ptr_old = ptr;
      ptr = strchr(ptr, ' ');
    }
    while (ptr) {
      if (i >= v.length()) v.set_length(2*v.length(), true);
      v[i++] = atof(ptr);

      ptr_old = ptr;
      ptr = strchr(ptr, ' ');
      while (ptr == ptr_old) {
        ptr++;
        ptr_old = ptr;
        ptr = strchr(ptr, ' ');
      }
    }
    flag = true;
    flag = false;
    v.set_length(i, true);
  }
  return ifs;
}


void Vector_Quantizer::load(const char *Name)
{
  vec   Temp;
  ifstream CodeBookFile(Name);
  vec   v;
  int   n;
  int   d;

  it_error_if(!CodeBookFile, std::string("Vector_Quantizer::load : cannot open file ") + Name);
  cout << "Reading the codebook " << Name ;
  cout.flush() ;
  CodeBookFile >> v ;
  d = length(v);
  Temp.set_length(d*16);
  n = 0;
  while (!CodeBookFile.eof()) {
    if (n*d >= Temp.length()) Temp.set_length(2*Temp.length(), true);
    Temp.replace_mid(n*d, v);
    n++;
    CodeBookFile >> v ;
  }
  Size = n;
  Dim = d;
  CodeBook.set_length(Size*Dim);
  for (n = 0;n < CodeBook.length();n++) CodeBook(n) = Temp(n);
  cout << "  size:" << size() << "  dim:" << dim() << endl ;
}

void Vector_Quantizer::save(const char *Name) const
{
  ofstream CodeBookFile(Name);

  cout << "Saving the codebook " << Name << endl ;
  for (int i = 0;i < Size;i++) {
    vec v = CodeBook.mid(i * Dim, Dim);
    for (int j = 0;j < v.length();j++) {
      CodeBookFile.operator << (v[j]);
      if (j < v.length() - 1) CodeBookFile.put(' ') ;
    }
    CodeBookFile << endl ;
  }
  CodeBookFile.close();
}

void Vector_Quantizer::modify_codevector(int no, double mul, const vec &add)
{
  int    pos = Dim * no;

  for (int i = 0;i < Dim;i++) {
    CodeBook._elem(pos + i) *= mul;
    CodeBook._elem(pos + i) += add[i];
  }
}

vec Vector_Quantizer::get_codevector(int Index) const
{
  return CodeBook.mid(Index*Dim, Dim);
}

void Vector_Quantizer::set_codevector(int Index, const vec &v)
{
  it_error_if(Dim != length(v), "Vector_Quantizer::set_codevector : Wrong dimension");
  for (int i = 0;i < length(v);i++) {
    CodeBook._elem(Index*Dim + i) = v._elem(i);
  }
}

void Vector_Quantizer::set_codebook(const mat &CB)
{
  Size = CB.cols();
  Dim = CB.rows();
  CodeBook.set_length(Size*Dim);
  for (int i = 0;i < Size;i++) {
    for (int j = 0;j < Dim;j++) {
      CodeBook(i*Dim + j) = CB(j, i);
    }
  }
}

mat Vector_Quantizer::get_codebook() const
{
  mat CB(Dim, Size);

  for (int i = 0;i < Size;i++) {
    for (int j = 0;i < Dim;i++) {
      CB(j, i) = CodeBook(i * Dim + j);
    }
  }
  return CB;
}

//--------------------------------------------------------------------
//    class SQ
//--------------------------------------------------------------------

Scalar_Quantizer::Scalar_Quantizer()
{
}

// SQ(const char *Name);

int Scalar_Quantizer::encode(double x) const
{
  int il = 0, ih = Levels.length() - 1, im;

  while (il < ih - 1) {
    im = (il + ih) / 2;
    if (x < Levels(im)) ih = im;
    else il = im;
  }
  if (Levels(ih) - x < x - Levels(il)) return ih;
  else return il;
}

ivec Scalar_Quantizer::encode(const vec &x) const
{
  int  i;
  ivec Index(x.length());

  for (i = 0;i < x.length();i++) {
    Index(i) = encode(x(i));
  }
  return Index;
}

vec Scalar_Quantizer::decode(const ivec &Index) const
{
  int i;
  vec y(Index.length());

  for (i = 0;i < Index.length();i++) {
    y(i) = decode(Index(i));
  }
  return y;
}

vec Scalar_Quantizer::Q(const vec &x) const
{
  int i;
  vec y(x.length());

  for (i = 0;i < x.length();i++) {
    y(i) = Q(x(i));
  }
  return y;
}

// void load(const char *Name);
// void save(const char *Name) const;


//-------------------------------------------------------------------------


int scalar_encode(double x, vec &Levels)
{
  int il = 0, ih = Levels.length() - 1, im;

  while (il < ih - 1) {
    im = (il + ih) / 2;
    if (x < Levels(im)) ih = im;
    else il = im;
  }
  if (Levels(ih) - x < x - Levels(il)) return ih;
  else return il;
}

ivec scalar_encode(vec &x, vec &Levels)
{
  ivec ind(x.length());
  for (int i = 0;i < x.length();i++) ind(i) = scalar_encode(x(i), Levels);
  return ind;
}

} // namespace itpp

//! \endcond
