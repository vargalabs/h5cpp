/*!
 * \file
 * \brief Implementation of a Gaussian Mixture Model Class
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

#include <itpp/srccode/gmm.h>
#include <itpp/srccode/vqtrain.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/random.h>
#include <itpp/base/timing.h>
#include <iostream>
#include <fstream>

//! \cond

namespace itpp
{

GMM::GMM()
{
  d = 0;
  M = 0;
}

GMM::GMM(std::string filename)
{
  load(filename);
}

GMM::GMM(int M_in, int d_in)
{
  M = M_in;
  d = d_in;
  m = zeros(M * d);
  sigma = zeros(M * d);
  w = 1. / M * ones(M);

  for (int i = 0;i < M;i++) {
    w(i) = 1.0 / M;
  }
  compute_internals();
}

void GMM::init_from_vq(const vec &codebook, int dim)
{

  mat  C(dim, dim);
  int  i;
  vec  v;

  d = dim;
  M = codebook.length() / dim;

  m = codebook;
  w = ones(M) / double(M);

  C.clear();
  for (i = 0;i < M;i++) {
    v = codebook.mid(i * d, d);
    C = C + outer_product(v, v);
  }
  C = 1. / M * C;
  sigma.set_length(M*d);
  for (i = 0;i < M;i++) {
    sigma.replace_mid(i*d, diag(C));
  }

  compute_internals();
}

void GMM::init(const vec &w_in, const mat &m_in, const mat &sigma_in)
{
  int  i, j;
  d = m_in.rows();
  M = m_in.cols();

  m.set_length(M*d);
  sigma.set_length(M*d);
  for (i = 0;i < M;i++) {
    for (j = 0;j < d;j++) {
      m(i*d + j) = m_in(j, i);
      sigma(i*d + j) = sigma_in(j, i);
    }
  }
  w = w_in;

  compute_internals();
}

void GMM::set_mean(const mat &m_in)
{
  int  i, j;

  d = m_in.rows();
  M = m_in.cols();

  m.set_length(M*d);
  for (i = 0;i < M;i++) {
    for (j = 0;j < d;j++) {
      m(i*d + j) = m_in(j, i);
    }
  }
  compute_internals();
}

void GMM::set_mean(int i, const vec &means, bool compflag)
{
  m.replace_mid(i*length(means), means);
  if (compflag) compute_internals();
}

void GMM::set_covariance(const mat &sigma_in)
{
  int  i, j;

  d = sigma_in.rows();
  M = sigma_in.cols();

  sigma.set_length(M*d);
  for (i = 0;i < M;i++) {
    for (j = 0;j < d;j++) {
      sigma(i*d + j) = sigma_in(j, i);
    }
  }
  compute_internals();
}

void GMM::set_covariance(int i, const vec &covariances, bool compflag)
{
  sigma.replace_mid(i*length(covariances), covariances);
  if (compflag) compute_internals();
}

void GMM::marginalize(int d_new)
{
  it_error_if(d_new > d, "GMM.marginalize: cannot change to a larger dimension");

  vec  mnew(d_new*M), sigmanew(d_new*M);
  int  i, j;

  for (i = 0;i < M;i++) {
    for (j = 0;j < d_new;j++) {
      mnew(i*d_new + j) = m(i * d + j);
      sigmanew(i*d_new + j) = sigma(i * d + j);
    }
  }
  m = mnew;
  sigma = sigmanew;
  d = d_new;

  compute_internals();
}

void GMM::join(const GMM &newgmm)
{
  if (d == 0) {
    w = newgmm.w;
    m = newgmm.m;
    sigma = newgmm.sigma;
    d = newgmm.d;
    M = newgmm.M;
  }
  else {
    it_error_if(d != newgmm.d, "GMM.join: cannot join GMMs of different dimension");

    w = concat(double(M) / (M + newgmm.M) * w, double(newgmm.M) / (M + newgmm.M) * newgmm.w);
    w = w / sum(w);
    m = concat(m, newgmm.m);
    sigma = concat(sigma, newgmm.sigma);

    M = M + newgmm.M;
  }
  compute_internals();
}

void GMM::clear()
{
  w.set_length(0);
  m.set_length(0);
  sigma.set_length(0);
  d = 0;
  M = 0;
}

void GMM::save(std::string filename)
{
  std::ofstream f(filename.c_str());
  int   i, j;

  f << M << " " << d << std::endl ;
  for (i = 0;i < w.length();i++) {
    f << w(i) << std::endl ;
  }
  for (i = 0;i < M;i++) {
    f << m(i*d) ;
    for (j = 1;j < d;j++) {
      f << " " << m(i*d + j) ;
    }
    f << std::endl ;
  }
  for (i = 0;i < M;i++) {
    f << sigma(i*d) ;
    for (j = 1;j < d;j++) {
      f << " " << sigma(i*d + j) ;
    }
    f << std::endl ;
  }
}

void GMM::load(std::string filename)
{
  std::ifstream GMMFile(filename.c_str());
  int   i, j;

  it_error_if(!GMMFile, std::string("GMM::load : cannot open file ") + filename);

  GMMFile >> M >> d ;


  w.set_length(M);
  for (i = 0;i < M;i++) {
    GMMFile >> w(i) ;
  }
  m.set_length(M*d);
  for (i = 0;i < M;i++) {
    for (j = 0;j < d;j++) {
      GMMFile >> m(i*d + j) ;
    }
  }
  sigma.set_length(M*d);
  for (i = 0;i < M;i++) {
    for (j = 0;j < d;j++) {
      GMMFile >> sigma(i*d + j) ;
    }
  }
  compute_internals();
  std::cout << "  mixtures:" << M << "  dim:" << d << std::endl ;
}

double GMM::likelihood(const vec &x)
{
  double fx = 0;
  int  i;

  for (i = 0;i < M;i++) {
    fx += w(i) * likelihood_aposteriori(x, i);
  }
  return fx;
}

vec GMM::likelihood_aposteriori(const vec &x)
{
  vec  v(M);
  int  i;

  for (i = 0;i < M;i++) {
    v(i) = w(i) * likelihood_aposteriori(x, i);
  }
  return v;
}

double GMM::likelihood_aposteriori(const vec &x, int mixture)
{
  int  j;
  double s;

  it_error_if(d != x.length(), "GMM::likelihood_aposteriori : dimensions does not match");
  s = 0;
  for (j = 0;j < d;j++) {
    s += normexp(mixture * d + j) * sqr(x(j) - m(mixture * d + j));
  }
  return normweight(mixture)*std::exp(s);;
}

void GMM::compute_internals()
{
  int  i, j;
  double s;
  double constant = 1.0 / std::pow(2 * pi, d / 2.0);

  normweight.set_length(M);
  normexp.set_length(M*d);

  for (i = 0;i < M;i++) {
    s = 1;
    for (j = 0;j < d;j++) {
      normexp(i*d + j) = -0.5 / sigma(i * d + j);  // check time
      s *= sigma(i * d + j);
    }
    normweight(i) = constant / std::sqrt(s);
  }

}

vec GMM::draw_sample()
{
  static bool first = true;
  static vec cumweight;
  double u = randu();
  int  k;

  if (first) {
    first = false;
    cumweight = cumsum(w);
    it_error_if(std::abs(cumweight(length(cumweight) - 1) - 1) > 1e-6, "weight does not sum to 0");
    cumweight(length(cumweight) - 1) = 1;
  }
  k = 0;
  while (u > cumweight(k)) k++;

  return elem_mult(sqrt(sigma.mid(k*d, d)), randn(d)) + m.mid(k*d, d);
}

GMM gmmtrain(Array<vec> &TrainingData, int M, int NOITER, bool VERBOSE)
{
  mat   mean;
  int   i, j, d = TrainingData(0).length();
  vec   sig;
  GMM   gmm(M, d);
  vec   m(d*M);
  vec   sigma(d*M);
  vec   w(M);
  vec   normweight(M);
  vec   normexp(d*M);
  double  LL = 0, LLold, fx;
  double  constant = 1.0 / std::pow(2 * pi, d / 2.0);
  int   T = TrainingData.length();
  vec   x1;
  int   t, n;
  vec   msum(d*M);
  vec   sigmasum(d*M);
  vec   wsum(M);
  vec   p_aposteriori(M);
  vec   x2;
  double  s;
  vec   temp1, temp2;
  //double  MINIMUM_VARIANCE=0.03;

  //-----------initialization-----------------------------------

  mean = vqtrain(TrainingData, M, 200000, 0.5, VERBOSE);
  for (i = 0;i < M;i++) gmm.set_mean(i, mean.get_col(i), false);
  // for (i=0;i<M;i++) gmm.set_mean(i,TrainingData(randi(0,TrainingData.length()-1)),false);
  sig = zeros(d);
  for (i = 0;i < TrainingData.length();i++) sig += sqr(TrainingData(i));
  sig /= TrainingData.length();
  for (i = 0;i < M;i++) gmm.set_covariance(i, 0.5*sig, false);

  gmm.set_weight(1.0 / M*ones(M));

  //-----------optimization-----------------------------------

  tic();
  for (i = 0;i < M;i++) {
    temp1 = gmm.get_mean(i);
    temp2 = gmm.get_covariance(i);
    for (j = 0;j < d;j++) {
      m(i*d + j) = temp1(j);
      sigma(i*d + j) = temp2(j);
    }
    w(i) = gmm.get_weight(i);
  }
  for (n = 0;n < NOITER;n++) {
    for (i = 0;i < M;i++) {
      s = 1;
      for (j = 0;j < d;j++) {
        normexp(i*d + j) = -0.5 / sigma(i * d + j);  // check time
        s *= sigma(i * d + j);
      }
      normweight(i) = constant * w(i) / std::sqrt(s);
    }
    LLold = LL;
    wsum.clear();
    msum.clear();
    sigmasum.clear();
    LL = 0;
    for (t = 0;t < T;t++) {
      x1 = TrainingData(t);
      x2 = sqr(x1);
      fx = 0;
      for (i = 0;i < M;i++) {
        s = 0;
        for (j = 0;j < d;j++) {
          s += normexp(i * d + j) * sqr(x1(j) - m(i * d + j));
        }
        p_aposteriori(i) = normweight(i) * std::exp(s);
        fx += p_aposteriori(i);
      }
      p_aposteriori /= fx;
      LL = LL + std::log(fx);

      for (i = 0;i < M;i++) {
        wsum(i) += p_aposteriori(i);
        for (j = 0;j < d;j++) {
          msum(i*d + j) += p_aposteriori(i) * x1(j);
          sigmasum(i*d + j) += p_aposteriori(i) * x2(j);
        }
      }
    }
    for (i = 0;i < M;i++) {
      for (j = 0;j < d;j++) {
        m(i*d + j) = msum(i * d + j) / wsum(i);
        sigma(i*d + j) = sigmasum(i * d + j) / wsum(i) - sqr(m(i * d + j));
      }
      w(i) = wsum(i) / T;
    }
    LL = LL / T;

    if (std::abs((LL - LLold) / LL) < 1e-6) break;
    if (VERBOSE) {
      std::cout << n << ":   " << LL << "   " << std::abs((LL - LLold) / LL) << "   " << toc() <<  std::endl ;
      std::cout << "---------------------------------------" << std::endl ;
      tic();
    }
    else {
      std::cout << n << ": LL =  " << LL << "   " << std::abs((LL - LLold) / LL) << "\r" ;
      std::cout.flush();
    }
  }
  for (i = 0;i < M;i++) {
    gmm.set_mean(i, m.mid(i*d, d), false);
    gmm.set_covariance(i, sigma.mid(i*d, d), false);
  }
  gmm.set_weight(w);
  return gmm;
}

} // namespace itpp

//! \endcond
