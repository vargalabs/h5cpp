/*!
 * \file
 * \brief Implementations of linear prediction functions, and conversion
 * between common representations of linear predictive parameters
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

#include <itpp/srccode/lpcfunc.h>
#include <itpp/base/matfunc.h>
#include <itpp/signal/sigfun.h>
#include <itpp/stat/misc_stat.h>
#include <iostream>

//! \cond

using std::cout;
using std::endl;

namespace itpp
{

// Autocorrelation sequence to reflection coefficients conversion.
vec ac2rc(const vec &ac);
// Autocorrelation sequence to prediction polynomial conversion.
vec ac2poly(const vec &ac);
// Inverse sine parameters to reflection coefficients conversion.
vec is2rc(const vec &is);
// Reflection coefficients to autocorrelation sequence conversion.
vec rc2ac(const vec &rc);
// Reflection coefficients to inverse sine parameters conversion.
vec rc2is(const vec &rc);

vec autocorr(const vec &x, int order)
{
  if (order < 0) order = x.size();

  vec R(order + 1);
  double sum;
  int i, j;

  for (i = 0;i < order + 1;i++) {
    sum = 0;
    for (j = 0;j < x.size() - i;j++) {
      sum += x[j] * x[j+i];
    }
    R[i] = sum;
  }
  return R;
}

vec levinson(const vec &R2, int order)
{
  vec R = R2;
  R[0] = R[0] * (1. + 1.e-9);

  if (order < 0) order = R.length() - 1;
  double k, alfa, s;
  double *any = new double[order+1];
  double *a = new double[order+1];
  int j, m;
  vec out(order + 1);

  a[0] = 1;
  alfa = R[0];
  if (alfa <= 0) {
    out.clear();
    out[0] = 1;
    return out;
  }
  for (m = 1;m <= order;m++) {
    s = 0;
    for (j = 1;j < m;j++) {
      s = s + a[j] * R[m-j];
    }

    k = -(R[m] + s) / alfa;
    if (fabs(k) >= 1.0) {
      cout << "levinson : panic! abs(k)>=1, order " << m << ". Aborting..." << endl ;
      for (j = m;j <= order;j++) {
        a[j] = 0;
      }
      break;
    }
    for (j = 1;j < m;j++) {
      any[j] = a[j] + k * a[m-j];
    }
    for (j = 1;j < m;j++) {
      a[j] = any[j];
    }
    a[m] = k;
    alfa = alfa * (1 - k * k);
  }
  for (j = 0;j < out.length();j++) {
    out[j] = a[j];
  }
  delete any;
  delete a;
  return out;
}

vec lpc(const vec &x, int order)
{
  return levinson(autocorr(x, order), order);
}

vec poly2ac(const vec &poly)
{
  vec  a = poly;
  int order = a.length() - 1;
  double alfa, s, *any = new double[order+1];
  int j, m;
  vec  r(order + 1);
  vec  k = poly2rc(a);

  it_error_if(a[0] != 1, "poly2ac : not an lpc filter");
  r[0] = 1;
  alfa = 1;
  for (m = 1;m <= order;m++) {
    s = 0;
    for (j = 1;j < m;j++) {
      s = s + a[j] * r[m-j];
    }
    r[m] = -s - alfa * k[m-1];
    for (j = 1;j < m;j++) {
      any[j] = a[j] + k[m-1] * a[m-j];
    }
    for (j = 1;j < m;j++) {
      a[j] = any[j];
    }
    a[m] = k[m-1];
    alfa = alfa * (1 - sqr(k[m-1]));
  }
  delete any;
  return r;
}

vec poly2rc(const vec &a)
{
  // a is [1 xx xx xx], a.size()=order+1
  int   m, i;
  int    order = a.size() - 1;
  vec k(order);
  vec any(order + 1), aold(a);

  for (m = order - 1;m > 0;m--) {
    k[m] = aold[m+1] ;
    if (fabs(k[m]) > 1) k[m] = 1.0 / k[m];
    for (i = 0;i < m;i++) {
      any[i+1] = (aold[i+1] - aold[m-i] * k[m]) / (1 - k[m] * k[m]);
    }
    aold = any;
  }
  k[0] = any[1];
  if (fabs(k[0]) > 1) k[0] = 1.0 / k[0];
  return k;
}

vec rc2poly(const vec &k)
{
  int  m, i;
  vec a(k.length() + 1), any(k.length() + 1);

  a[0] = 1;
  any[0] = 1;
  a[1] = k[0];
  for (m = 1;m < k.size();m++) {
    any[m+1] = k[m];
    for (i = 0;i < m;i++) {
      any[i+1] = a[i+1] + a[m-i] * k[m];
    }
    a = any;
  }
  return a;
}

vec rc2lar(const vec &k)
{
  short m;
  vec LAR(k.size());

  for (m = 0;m < k.size();m++) {
    LAR[m] = std::log((1 + k[m]) / (1 - k[m]));
  }
  return LAR;
}

vec lar2rc(const vec &LAR)
{
  short m;
  vec k(LAR.size());

  for (m = 0;m < LAR.size();m++) {
    k[m] = (std::exp(LAR[m]) - 1) / (std::exp(LAR[m]) + 1);
  }
  return k;
}

double FNevChebP_double(double  x, const double c[], int n)
{
  int i;
  double b0 = 0.0, b1 = 0.0, b2 = 0.0;

  for (i = n - 1; i >= 0; --i) {
    b2 = b1;
    b1 = b0;
    b0 = 2.0 * x * b1 - b2 + c[i];
  }
  return (0.5 * (b0 - b2 + c[0]));
}

double FNevChebP(double  x, const double c[], int n)
{
  int i;
  double b0 = 0.0, b1 = 0.0, b2 = 0.0;

  for (i = n - 1; i >= 0; --i) {
    b2 = b1;
    b1 = b0;
    b0 = 2.0 * x * b1 - b2 + c[i];
  }
  return (0.5 * (b0 - b2 + c[0]));
}

vec poly2lsf(const vec &pc)
{
  int np = pc.length() - 1;
  vec lsf(np);

  vec fa((np + 1) / 2 + 1), fb((np + 1) / 2 + 1);
  vec ta((np + 1) / 2 + 1), tb((np + 1) / 2 + 1);
  double *t;
  double xlow, xmid, xhigh;
  double ylow, ymid, yhigh;
  double xroot;
  double dx;
  int i, j, nf;
  int odd;
  int na, nb, n;
  double ss, aa;
  double DW = (0.02 * pi);
  int  NBIS = 4;

  odd = (np % 2 != 0);
  if (odd) {
    nb = (np + 1) / 2;
    na = nb + 1;
  }
  else {
    nb = np / 2 + 1;
    na = nb;
  }

  fa[0] = 1.0;
  for (i = 1, j = np; i < na; ++i, --j)
    fa[i] = pc[i] + pc[j];

  fb[0] = 1.0;
  for (i = 1, j = np; i < nb; ++i, --j)
    fb[i] = pc[i] - pc[j];

  if (odd) {
    for (i = 2; i < nb; ++i)
      fb[i] = fb[i] + fb[i-2];
  }
  else {
    for (i = 1; i < na; ++i) {
      fa[i] = fa[i] - fa[i-1];
      fb[i] = fb[i] + fb[i-1];
    }
  }

  ta[0] = fa[na-1];
  for (i = 1, j = na - 2; i < na; ++i, --j)
    ta[i] = 2.0 * fa[j];

  tb[0] = fb[nb-1];
  for (i = 1, j = nb - 2; i < nb; ++i, --j)
    tb[i] = 2.0 * fb[j];

  nf = 0;
  t = ta._data();
  n = na;
  xroot = 2.0;
  xlow = 1.0;
  ylow = FNevChebP_double(xlow, t, n);


  ss = std::sin(DW);
  aa = 4.0 - 4.0 * std::cos(DW)  - ss;
  while (xlow > -1.0 && nf < np) {
    xhigh = xlow;
    yhigh = ylow;
    dx = aa * xhigh * xhigh + ss;
    xlow = xhigh - dx;
    if (xlow < -1.0)
      xlow = -1.0;
    ylow = FNevChebP_double(xlow, t, n);
    if (ylow * yhigh <= 0.0) {
      dx = xhigh - xlow;
      for (i = 1; i <= NBIS; ++i) {
        dx = 0.5 * dx;
        xmid = xlow + dx;
        ymid = FNevChebP_double(xmid, t, n);
        if (ylow * ymid <= 0.0) {
          yhigh = ymid;
          xhigh = xmid;
        }
        else {
          ylow = ymid;
          xlow = xmid;
        }
      }
      if (yhigh != ylow)
        xmid = xlow + dx * ylow / (ylow - yhigh);
      else
        xmid = xlow + dx;
      lsf[nf] = std::acos((double) xmid);
      ++nf;
      if (xmid >= xroot) {
        xmid = xlow - dx;
      }
      xroot = xmid;
      if (t == ta._data()) {
        t = tb._data();
        n = nb;
      }
      else {
        t = ta._data();
        n = na;
      }
      xlow = xmid;
      ylow = FNevChebP_double(xlow, t, n);
    }
  }
  if (nf != np) {
    cout << "poly2lsf: WARNING: failed to find all lsfs" << endl ;
  }
  return lsf;
}

vec lsf2poly(const vec &f)
{
  int m = f.length();
  vec  pc(m + 1);
  double c1, c2, *a;
  vec  p(m + 1), q(m + 1);
  int mq, n, i, nor;

  it_error_if(m % 2 != 0, "lsf2poly: THIS ROUTINE WORKS ONLY FOR EVEN m");
  pc[0] = 1.0;
  a = pc._data() + 1;
  mq = m >> 1;
  for (i = 0 ; i <= m ; i++) {
    q[i] = 0.;
    p[i] = 0.;
  }
  p[0] = q[0] = 1.;
  for (n = 1; n <= mq; n++) {
    nor = 2 * n;
    c1 = 2 * std::cos(f[nor-1]);
    c2 = 2 * std::cos(f[nor-2]);
    for (i = nor; i >= 2; i--) {
      q[i] += q[i-2] - c1 * q[i-1];
      p[i] += p[i-2] - c2 * p[i-1];
    }
    q[1] -= c1;
    p[1] -= c2;
  }
  a[0] = 0.5 * (p[1] + q[1]);
  for (i = 1, n = 2; i < m ; i++, n++)
    a[i] = 0.5 * (p[i] + p[n] + q[n] - q[i]);

  return pc;
}

vec poly2cepstrum(const vec &a)
{
  vec c(a.length() - 1);

  for (int n = 1;n <= c.length();n++) {
    c[n-1] = a[n];
    for (int k = 1;k < n;k++) {
      c[n-1] -= double(k) / n * a[n-k] * c[k-1];
    }
  }
  return c;
}

vec poly2cepstrum(const vec &a, int num)
{
  it_error_if(num < a.length(), "a2cepstrum : not allowed cepstrum length");
  vec c(num);
  int n;

  for (n = 1;n < a.length();n++) {
    c[n-1] = a[n];
    for (int k = 1;k < n;k++) {
      c[n-1] -= double(k) / n * a[n-k] * c[k-1];
    }
  }
  for (n = a.length();n <= c.length();n++) {
    c[n-1] = 0;
    for (int k = n - a.length() + 1;k < n;k++) {
      c[n-1] -= double(k) / n * a[n-k] * c[k-1];
    }
  }
  return c;
}

vec cepstrum2poly(const vec &c)
{
  vec a(c.length() + 1);

  a[0] = 1;
  for (int n = 1;n <= c.length();n++) {
    a[n] = c[n-1];
    for (int k = 1;k < n;k++) {
      a[n] += double(k) / n * a[n-k] * c[k-1];
    }
  }
  return a;
}

vec chirp(const vec &a, double factor)
{
  vec    temp(a.length());
  int    i;
  double   f = factor;

  it_error_if(a[0] != 1, "chirp : a[0] should be 1");
  temp[0] = a[0];
  for (i = 1;i < a.length();i++) {
    temp[i] = a[i] * f;
    f *= factor;
  }
  return temp;
}

vec schurrc(const vec &R, int order)
{
  if (order == -1) order = R.length() - 1;

  vec    k(order), scratch(2*order + 2);

  int m;
  int h;
  double ex;
  double *ep;
  double *en;

  ep = scratch._data();
  en = scratch._data() + order + 1;

  m = 0;
  while (m < order) {
    m++;
    ep[m] = R[m];
    en[m] = R[m-1];
  }
  if (en[1] < 1.0) en[1] = 1.0;
  h = -1;
  while (h < order) {
    h++;
    k[h] = -ep[h+1] / en[1];
    en[1] = en[1] + k[h] * ep[h+1];
    if (h == (order - 1)) {
      // cout << "k: " << k << endl ;
      return k;
    }
    ep[order] = ep[order] + k[h] * en[order-h];
    m = h + 1;
    while (m < (order - 1)) {
      m++;
      ex = ep[m] + k[h] * en[m-h];
      en[m-h] = en[m-h] + k[h] * ep[m];
      ep[m] = ex;
    }
  }
  return k;  // can never come here
}

vec lerouxguegenrc(const vec &R, int order)
{
  vec    k(order);

  double  *r, *rny;
  int j, m;
  int M = order;

  r = new double[2*M+1];
  rny = new double[2*M+1];

  for (j = 0;j <= M;j++) {
    r[M-j] = r[M+j] = R[j];
  }
  for (m = 1;m <= M;m++) {
    k[m-1] = -r[M+m] / r[M];
    for (j = -M;j <= M;j++) {
      rny[M+j] = r[M+j] + k[m-1] * r[M+m-j];
    }
    for (j = -M;j <= M;j++) {
      r[M+j] = rny[M+j];
    }
  }
  delete r;
  delete rny;
  return k;
}

double sd(const vec &In1, const vec &In2)
{
  return std::sqrt(37.722339402*energy(poly2cepstrum(In1, 32) - poly2cepstrum(In2, 32)));
}

// highestfreq=1 gives entire band
double sd(const vec &In1, const vec &In2, double highestfreq)
{
  vec Diff = sqr(abs(log10(filter_spectrum(In1, In2))));
  double S = 0;
  for (int i = 0;i < round(highestfreq*129);i++) {
    S = S + Diff(i);
  }
  S = S * 100 / round(highestfreq * 129);
  return std::sqrt(S);
}

} // namespace itpp

//! \endcond
