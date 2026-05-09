/*!
 * \file
 * \brief Implementation of a vector quantizer training functions
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

#include <itpp/srccode/vqtrain.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/random.h>
#include <itpp/base/sort.h>
#include <itpp/base/specmat.h>
#include <itpp/base/math/min_max.h>
#include <itpp/stat/misc_stat.h>
#include <iostream>

//! \cond

namespace itpp
{

// the cols contains codevectors
double kmeansiter(Array<vec> &DB, mat &codebook)
{
  int    DIM = DB(0).length(), SIZE = codebook.cols(), T = DB.length();
  vec    x, xnum(SIZE);
  mat    xsum(DIM, SIZE);
  int    MinIndex, i, j, k;
  double   MinS, S, D, *xp, *cp;

  xsum.clear();
  xnum.clear();

  D = 1E20;
  D = 0;
  for (k = 0;k < T;k++) {
    x = DB(k);
    xp = x._data();
    MinS = 1E20;
    MinIndex = 0;
    for (i = 0;i < SIZE;i++) {
      cp = &codebook(0, i);
      S = sqr(xp[0] - cp[0]);
      for (j = 1;j < DIM;j++) {
        S += sqr(xp[j] - cp[j]);
        if (S >= MinS) goto sune;
      }
      MinS = S;
      MinIndex = i;
    sune:
    void();
    }
    D += MinS;
    cp = &xsum(0, MinIndex);
    for (j = 0;j < DIM;j++) {
      cp[j] += xp[j];
    }
    xnum(MinIndex)++;
  }
  for (i = 0;i < SIZE;i++) {
    for (j = 0;j < DIM;j++) {
      codebook(j, i) = xsum(j, i) / xnum(i);
    }
  }
  return D;
}

mat kmeans(Array<vec> &DB, int SIZE, int NOITER, bool VERBOSE)
{
  int    DIM = DB(0).length(), T = DB.length();
  mat    codebook(DIM, SIZE);
  int    n, i, j;
  double   D, Dold;
  ivec   ind(SIZE);

  for (i = 0;i < SIZE;i++) {
    ind(i) = randi(0, T - 1);
    j = 0;
    while (j < i) {
      if (ind(j) == ind(i)) {
        ind(i) = randi(0, T - 1);
        j = 0;
      }
      j++;
    }
    codebook.set_col(i, DB(ind(i)));
  }


  if (VERBOSE) std::cout << "Training VQ..." << std::endl ;

  D = 1E20;
  for (n = 0;n < NOITER;n++) {
    Dold = D;
    D = kmeansiter(DB, codebook);
    if (VERBOSE) std::cout << n << ": " << D / T << " ";
    if (std::abs((D - Dold) / D) < 1e-4) break;
  }
  return codebook;
}

mat lbg(Array<vec> &DB, int SIZE, int NOITER, bool VERBOSE)
{
  int  S = 1, DIM = DB(0).length(), T = DB.length(), i, n;
  mat  cb;
  vec  delta = 0.001 * randn(DIM), x;
  double D, Dold;

  x = zeros(DIM);
  for (i = 0;i < T;i++) {
    x += DB(i);
  }
  x /= T;
  cb.set_size(DIM, 1);
  cb.set_col(0, x);
  while (cb.cols() < SIZE) {
    S = cb.cols();
    if (2*S <= SIZE) cb.set_size(DIM, 2*S, true);
    else cb.set_size(DIM, 2*S, true);
    for (i = S;i < cb.cols();i++) {
      cb.set_col(i, cb.get_col(i - S) + delta);
    }
    D = 1E20;
    for (n = 0;n < NOITER;n++) {
      Dold = D;
      D = kmeansiter(DB, cb);
      if (VERBOSE) std::cout << n << ": " << D / T << " ";
      if (std::abs((D - Dold) / D) < 1e-4) break;
    }
  }

  return cb;
}

mat vqtrain(Array<vec> &DB, int SIZE, int NOITER, double STARTSTEP, bool VERBOSE)
{
  int    DIM = DB(0).length();
  vec    x;
  vec    codebook(DIM*SIZE);
  int    n, MinIndex, i, j;
  double   MinS, S, D, step, *xp, *cp;

  for (i = 0;i < SIZE;i++) {
    codebook.replace_mid(i*DIM, DB(randi(0, DB.length() - 1)));
  }
  if (VERBOSE) std::cout << "Training VQ..." << std::endl ;

res:
  D = 0;
  for (n = 0;n < NOITER;n++) {
    step = STARTSTEP * (1.0 - double(n) / NOITER);
    if (step < 0) step = 0;
    x = DB(randi(0, DB.length() - 1)); // seems unnecessary! Check it up.
    xp = x._data();

    MinS = 1E20;
    MinIndex = 0;
    for (i = 0;i < SIZE;i++) {
      cp = &codebook(i * DIM);
      S = sqr(xp[0] - cp[0]);
      for (j = 1;j < DIM;j++) {
        S += sqr(xp[j] - cp[j]);
        if (S >= MinS) goto sune;
      }
      MinS = S;
      MinIndex = i;
    sune:
      i = i;
    }
    D += MinS;
    cp = &codebook(MinIndex * DIM);
    for (j = 0;j < DIM;j++) {
      cp[j] += step * (xp[j] - cp[j]);
    }
    if ((n % 20000 == 0) && (n > 1)) {
      if (VERBOSE) std::cout << n << ": " << D / 20000 << " ";
      D = 0;
    }
  }

  // checking training result
  vec dist(SIZE), num(SIZE);

  dist.clear();
  num.clear();
  for (n = 0;n < DB.length();n++) {
    x = DB(n);
    xp = x._data();
    MinS = 1E20;
    MinIndex = 0;
    for (i = 0;i < SIZE;i++) {
      cp = &codebook(i * DIM);
      S = sqr(xp[0] - cp[0]);
      for (j = 1;j < DIM;j++) {
        S += sqr(xp[j] - cp[j]);
        if (S >= MinS) goto sune2;
      }
      MinS = S;
      MinIndex = i;
    sune2:
      i = i;
    }
    dist(MinIndex) += MinS;
    num(MinIndex) += 1;
  }
  dist = 10 * log10(dist * dist.length() / sum(dist));
  if (VERBOSE) std::cout << std::endl << "Distortion contribution: " << dist << std::endl ;
  if (VERBOSE) std::cout << "Num spread: " << num / DB.length()*100 << " %" << std::endl << std::endl ;
  if (min(dist) < -30) {
    std::cout << "Points without entries! Retraining" << std::endl ;
    j = min_index(dist);
    i = max_index(num);
    codebook.replace_mid(j*DIM, codebook.mid(i*DIM, DIM));
    goto res;
  }

  mat cb(DIM, SIZE);
  for (i = 0;i < SIZE;i++) {
    cb.set_col(i, codebook.mid(i*DIM, DIM));
  }
  return cb;
}

vec sqtrain(const vec &inDB, int SIZE)
{
  vec  DB(inDB);
  vec  Levels, Levels_old;
  ivec indexlist(SIZE + 1);
  int  il, im, ih, i;
  int  SIZEDB = inDB.length();
  double x;

  sort(DB);
  Levels = DB(round_i(linspace(0.01 * SIZEDB, 0.99 * SIZEDB, SIZE)));
  Levels_old = zeros(SIZE);

  while (energy(Levels - Levels_old) > 0.0001) {
    Levels_old = Levels;
    for (i = 0;i < SIZE - 1;i++) {
      x = (Levels(i) + Levels(i + 1)) / 2;
      il = 0;
      ih = SIZEDB - 1;
      while (il < ih - 1) {
        im = (il + ih) / 2;
        if (x < DB(im)) ih = im;
        else il = im;
      }
      indexlist(i + 1) = il;
    }
    indexlist(0) = -1;
    indexlist(SIZE) = SIZEDB - 1;
    for (i = 0;i < SIZE;i++) Levels(i) = mean(DB(indexlist(i) + 1, indexlist(i + 1)));
  }
  return Levels;
}

ivec bitalloc(const vec &variances, int nobits)
{
  ivec bitvec(variances.length());
  bitvec.clear();
  int  i, bits = nobits;
  vec  var = variances;

  while (bits > 0) {
    i = max_index(var);
    var(i) /= 4;
    bitvec(i)++;
    bits--;
  }
  return bitvec;
}

} // namespace itpp

//! \endcond
