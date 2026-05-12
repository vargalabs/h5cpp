/*!
 * \file
 * \brief Implementation of multilateration class for indoor localization
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

#include <itpp/comm/multilateration.h>
#ifdef _MSC_VER
#include <float.h>
#define isfinite _finite
#endif

using namespace std;

namespace itpp
{

struct Point {
  double x;
  double y;
  double z;
};

//class for managing a set of integers (indices)
class IndexSet
{
public:
  IndexSet(unsigned int max_size) :
    size_(0), max_size_(max_size), data_(new unsigned int[max_size]) {}
  unsigned int get_size() const {
    return size_;
  }
  ~IndexSet() {
    delete[] data_;
  }
  unsigned int& operator[](unsigned int n) {
    return data_[n];
  }
  unsigned int operator[](unsigned int n) const {
    return data_[n];
  }
  void set_size(unsigned int size) {
    size_ = size;
  }
  bool exist(unsigned int elem) {
    bool res = false;
    for(unsigned int i = 0; i < size_; ++i) {
      if(data_[i] == elem) {
        res = true;
        break;
      }
    }
    return res;
  }
  bool insert(unsigned int elem) {
    bool res = false;
    if((size_ < max_size_) && (false == exist(elem))) {
      data_[size_++] = elem;
      res = true;
    }
    return res;
  }
  bool remove(unsigned int elem) {
    bool res = false;
    for(unsigned int i = 0; i < size_; ++i) {
      if(data_[i] == elem) {
        if((i + 1) < size_) {
          memcpy(data_ + i, data_ + i + 1, (size_ - i - 1)*sizeof(*data_));
        }
        --size_;
        res = true;
        break;
      }
    }
    return res;
  }
private:
  unsigned int size_;
  unsigned int max_size_;
  unsigned int *data_;
};

//helper class
class PosHelper
{
public:
  static double get_dist(const Point *p0, const Point *p1) {
    double out = 0.0;
    double *coord_p0 = (double*)p0;
    double *coord_p1 = (double*)p1;
    for(int n = 0; n < 3; ++n) {
      out += (coord_p1[n] - coord_p0[n]) * (coord_p1[n] - coord_p0[n]);
    }
    return std::sqrt(out);
  }
  static bool test_uniqueness(const Point *bs_pos, int bs_pos_len, double eps) {
    int i, j;

    for(i = bs_pos_len - 1; i >= 0; --i) {
      for(j = i - 1; j >= 0; --j) {
        if((fabs(bs_pos[i].x - bs_pos[j].x) < eps) && (fabs(bs_pos[i].y - bs_pos[j].y) < eps) &&
            (fabs(bs_pos[i].z - bs_pos[j].z) < eps)) {
          it_warning("found 2 identical BSs");
          return false;
        }
      }
    }
    return true;
  }
  static double det_3by3(const double mat[9]) {
    return mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]) + mat[3] * (mat[2] * mat[7] - mat[1] * mat[8]) + mat[6] * (mat[1] * mat[5] - mat[2] * mat[4]);
  }
  static double det_4by4(const double mat[16]) {
    return mat[0] * mat[5] * mat[10] * mat[15] + mat[0] * mat[9] * mat[14] * mat[7] + mat[0] * mat[13] * mat[6] * mat[11] +
           mat[4] * mat[1] * mat[14] * mat[11] + mat[4] * mat[9] * mat[2] * mat[15] + mat[4] * mat[13] * mat[10] * mat[3] +
           mat[8] * mat[1] * mat[6] * mat[15] + mat[8] * mat[5] * mat[14] * mat[3] + mat[8] * mat[13] * mat[2] * mat[7] +
           mat[12] * mat[1] * mat[10] * mat[7] + mat[12] * mat[5] * mat[2] * mat[11] + mat[12] * mat[9] * mat[6] * mat[3] -
           mat[0] * mat[5] * mat[14] * mat[11] - mat[0] * mat[9] * mat[6] * mat[15] - mat[0] * mat[13] * mat[10] * mat[7] -
           mat[4] * mat[1] * mat[10] * mat[15] - mat[4] * mat[9] * mat[14] * mat[3] - mat[4] * mat[13] * mat[2] * mat[11] -
           mat[8] * mat[1] * mat[14] * mat[7] - mat[8] * mat[5] * mat[2] * mat[15] - mat[8] * mat[13] * mat[6] * mat[3] -
           mat[12] * mat[1] * mat[6] * mat[11] - mat[12] * mat[5] * mat[10] * mat[3] - mat[12] * mat[9] * mat[2] * mat[7];
  }
  static bool inv_3by3(double *inv, const double *mat) {
    double det = det_3by3(mat);
    if((0 == det) || (0 == isfinite(det))) {
      memset(inv, 0, 9 * sizeof(*inv));
      return false;
    }
    inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / det;
    inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / det;
    inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / det;
    inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / det;
    inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / det;
    inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / det;
    inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / det;
    inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / det;
    inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / det;

    return true;
  }
  static bool test_noncoplanar(const Point *bs_pos, int bs_pos_len) {
    double mat[3 * 3];
    double *coord[] = {NULL, NULL};
    int i;
    int j;

    if((NULL == bs_pos) || (4 != bs_pos_len)) {
      it_warning("4 points are needed");
      return false;
    }

    coord[0] = (double*)bs_pos;
    for(i = 0; i < 3; ++i) {
      coord[1] = (double*)(bs_pos + i + 1);
      for(j = 0; j < 3; ++j) {
        mat[i + 3 * j] = (*(coord[1] + j)) - (*(coord[0] + j));
      }
    }
    return (0 != det_3by3(mat)) ? true : false;
  }
  static bool combination(unsigned int *comb_mat, unsigned int *nb_cols, unsigned int nb_rows, unsigned int len) {
    unsigned int idx = 0;
    unsigned int j, k, l, m;

    if((NULL == comb_mat) || (NULL == nb_cols)) {
      it_warning("invalid input");
      return false;
    }
    if(nb_rows > len) {
      //it_warning("number of rows superior to total length");
      return false;
    }

    switch(nb_rows) {
    case 5: {
      for(j = 1; j < len; ++j) {
        for(k = j + 1; k < len; ++k) {
          for(l = k + 1; l < len; ++l) {
            for(m = l + 1; m < len; ++m) {
              comb_mat[idx++] = 0;
              comb_mat[idx++] = j;
              comb_mat[idx++] = k;
              comb_mat[idx++] = l;
              comb_mat[idx++] = m;
            }
          }
        }
      }
    }
    break;
    case 4: {
      for(j = 1; j < len; ++j) {
        for(k = j + 1; k < len; ++k) {
          for(l = k + 1; l < len; ++l) {
            comb_mat[idx++] = 0;
            comb_mat[idx++] = j;
            comb_mat[idx++] = k;
            comb_mat[idx++] = l;
          }
        }
      }
    }
    break;
    case 3: {
      for(j = 1; j < len; ++j) {
        for(k = j + 1; k < len; ++k) {
          comb_mat[idx++] = 0;
          comb_mat[idx++] = j;
          comb_mat[idx++] = k;
        }
      }
    }
    break;
    case 2: {
      for(j = 1; j < len; ++j) {
        comb_mat[idx++] = 0;
        comb_mat[idx++] = j;
      }
    }
    break;
    default:
      it_warning("invalid nb_rows");
      return false;
    }

    *nb_cols = idx / nb_rows;
    return true;
  }
  static unsigned int comb_nb(unsigned int n, unsigned int k) {
    unsigned int i;
    unsigned int nom = 1;
    unsigned int denom = 1;

    for(i = k + 1; i <= n; ++i) {
      nom *= i;
    }

    for(i = 2; i <= (n - k); ++i) {
      denom *= i;
    }
    return nom / denom;
  }
  static bool update_subsets(IndexSet &used_subset_idx, IndexSet &unused_subset_idx,
                             const unsigned int *subset_idx, unsigned int subset_idx_len) {
    unsigned int i;
    bool res = true;
    for(i = 0; i < subset_idx_len; ++i) {
      if((true == used_subset_idx.insert(subset_idx[i])) &&
          (false == unused_subset_idx.remove(subset_idx[i]))) {
        it_warning("cannot update subsets");
        res = false;
        break;
      }
    }
    return res;
  }
  static bool store_subset(unsigned int **subsets, unsigned int *in_subset, unsigned int in_subset_len, unsigned int in_subset_idx) {
    bool rc = false;

    /* realloc memory */
    *subsets = (unsigned int*)realloc(*subsets, (in_subset_idx + 1) * in_subset_len * sizeof(**subsets));
    if(NULL != *subsets) {
      /* copy subset */
      memcpy((*subsets) + in_subset_idx * in_subset_len, in_subset, in_subset_len * sizeof(*in_subset));
      rc = true;
    }
    else {
      it_warning("cannot allocate memory");
    }
    return rc;
  }
};

// interface for multilateration algorithms: spherical or hyperbolic
class Algorithm
{
public:
  enum GeoPos {GEO_POS_NO = 0, GEO_POS_YES, GEO_POS_EXACT};
  virtual GeoPos validate(const Point *bs_pos, unsigned int bs_pos_len, const double *r) = 0;
  virtual bool setup(const Point *bs_pos, unsigned int bs_pos_len) = 0;
  /* compute an unique MS position */
  virtual bool get_pos(Point *ms_pos, const double *range, unsigned int range_len) = 0;
  virtual void set_meas_mat(const double *meas) = 0;
  virtual void get_meas(double *meas, const unsigned int *bs_subset_idx, unsigned int len) = 0;
  virtual bool get_meas_mult_mat(const unsigned int *bs_pos_idx, unsigned int rows, unsigned int cols) = 0;
  virtual void get_meas_mult(unsigned int *meas_mult, const unsigned int *subset_idx, unsigned int meas_mult_len) = 0;
  virtual bool get_grad(double *grad, const Point *bs_pos, unsigned int bs_nb, const Point *ms_pos) = 0;
  virtual double* get_meas_mat() const = 0;
  virtual unsigned int* get_mult_mat() const = 0;
  virtual ~Algorithm() {}
};

class Spherical : public Algorithm
{
public:
  explicit Spherical(unsigned int nb_bs);
  virtual ~Spherical();
  virtual Algorithm::GeoPos validate(const Point *bs_pos, unsigned int bs_pos_len, const double *r);
  virtual bool setup(const Point *bs_pos, unsigned int bs_pos_len);
  virtual bool get_pos(Point *ms_pos, const double *range, unsigned int range_len);
  virtual void set_meas_mat(const double *meas) {
    memcpy(meas_mat_, meas, nb_bs_ * sizeof(*meas));
  }
  virtual void get_meas(double *meas, const unsigned int *bs_subset_idx, unsigned int len);
  virtual bool get_meas_mult_mat(const unsigned int *bs_pos_idx, unsigned int rows, unsigned int cols);
  virtual void get_meas_mult(unsigned int *meas_mult, const unsigned int *subset_idx, unsigned int meas_mult_len);
  virtual bool get_grad(double *grad, const Point *bs_pos, unsigned int bs_nb, const Point *ms_pos);
  virtual double* get_meas_mat() const {
    return meas_mat_;
  }
  virtual unsigned int* get_mult_mat() const {
    return mult_mat_;
  }
private:
  double mult_vTMw(double v[2], double M[4], double w[2]) {
    return w[0] * (M[0] * v[0] + M[1] * v[1]) + w[1] * (M[2] * v[0] + M[3] * v[1]);
  }
  double comp_quadratic(double *u, const double *range) {
    double v[10];
    double scalar;
    int i;

    v[0] = u[0] = 1;
    v[1] = u[1] = range[0] * range[0];
    v[2] = u[2] = range[1] * range[1];
    v[3] = u[3] = range[2] * range[2];

    v[4] = u[1] * u[1];
    v[5] = u[2] * u[2];
    v[6] = u[3] * u[3];
    v[7] = u[1] * u[2];
    v[8] = u[1] * u[3];
    v[9] = u[2] * u[3];

    /*could generate numerical overflows*/
    scalar = xi[0] * v[0];
    for(i = 1; i < 10; ++i) {
      scalar += xi[i] * v[i];
    }
    return scalar;
  }
  //data members
  Point bs_pos_;
  unsigned int bs_pos_len_;
  double xi[10];
  double Lambda[12];
  double mu[6];

  unsigned int nb_bs_;
  double *meas_mat_;
  unsigned int *mult_mat_;
};

Spherical::Spherical(unsigned int nb_bs) :
  nb_bs_(nb_bs),
  meas_mat_(new double[nb_bs]),
  mult_mat_(new unsigned int [nb_bs])
{
}

Spherical::~Spherical()
{
  delete[] meas_mat_;
  delete[] mult_mat_;
}

/*validate a given set of BSs and the associated measurements*/
Algorithm::GeoPos Spherical::validate(const Point *bs_pos, unsigned int bs_pos_len, const double *r)
{
  GeoPos out = GEO_POS_NO;

  if(3 > bs_pos_len) {
    it_warning("invalid input");
    return GEO_POS_NO;
  }
  if(true == Spherical::setup(bs_pos, bs_pos_len)) {
    double u[4];/*not used*/
    double scalar = comp_quadratic(u, r);
    if(0 == scalar) {
      out = GEO_POS_EXACT;
    }
    else if((0 < scalar) && (1 == bs_pos_len_)) {
      out = GEO_POS_YES;
    }
    else {
      out = GEO_POS_NO;
    }
  }
  return out;
}

bool Spherical::setup(const Point *bs_pos, unsigned int bs_pos_len)
{
  double S[3];
  double W[4];
  double detW;
  double invW[4];
  double G[4];/*W^{-T}*W^{-1}*/
  double d[2];
  double rh1[2];
  double a;
  double delta[2];
  double e[2];
  double lambda[4];
  double w1Td;
  double w2Td;
  int i;

  if(3 > bs_pos_len) {
    it_warning("at least 3 BSs are needed");
    return false;
  }

  memset(xi, 0, sizeof(xi));
  memset(Lambda, 0, sizeof(Lambda));
  memset(mu, 0, sizeof(mu));
  bs_pos_len_ = 0;

  for(i = 0; i < 3; ++i) {
    S[i] = bs_pos[i].x * bs_pos[i].x + bs_pos[i].y * bs_pos[i].y + bs_pos[i].z * bs_pos[i].z;
  }
  W[0] = bs_pos[1].x - bs_pos[0].x;
  W[1] = bs_pos[2].x - bs_pos[0].x;
  W[2] = bs_pos[1].y - bs_pos[0].y;
  W[3] = bs_pos[2].y - bs_pos[0].y;
  detW = W[0] * W[3] - W[1] * W[2];
  if(0 == detW) {
    //it_warning("please make sure that the 3 known points are not colinear");
    return false;
  }
  invW[0] = W[3] / detW;
  invW[1] = -W[1] / detW;
  invW[2] = -W[2] / detW;
  invW[3] = W[0] / detW;
  G[0] = invW[0] * invW[0] + invW[1] * invW[1];
  G[1] = invW[0] * invW[2] + invW[1] * invW[3];
  G[2] = G[1];
  G[3] = invW[2] * invW[2] + invW[3] * invW[3];
  rh1[0] = bs_pos[0].x;
  rh1[1] = bs_pos[0].y;
  for(i = 0; i < 2; ++i) {
    d[i] = bs_pos[i + 1].z - bs_pos[0].z;
    delta[i] = S[0] - S[i + 1];
  }
  for(i = 0; i < 2; ++i) {
    e[i] = 0.5 * (delta[0] * G[2 * i] + delta[1] * G[2 * i + 1]) + rh1[0] * invW[2 * i] + rh1[1] * invW[2 * i + 1];
  }
  a = 1.0 + mult_vTMw(d, G, d);

  lambda[0] = -(2 * mult_vTMw(rh1, invW, d) - 2 * bs_pos[0].z + mult_vTMw(d, G, delta)) / (2 * a);
  lambda[1] = ((bs_pos[1].z - bs_pos[0].z) * (G[0] + G[1]) + (bs_pos[2].z - bs_pos[0].z) * (G[1] + G[3])) / (2 * a);
  lambda[2] = -((bs_pos[1].z - bs_pos[0].z) * G[0] + (bs_pos[2].z - bs_pos[0].z) * G[1]) / (2 * a);
  lambda[3] = -((bs_pos[2].z - bs_pos[0].z) * G[3] + (bs_pos[1].z - bs_pos[0].z) * G[1]) / (2 * a);

  xi[0] = lambda[0] * lambda[0] - (0.25 * mult_vTMw(delta, G, delta) + mult_vTMw(rh1, invW, delta) + S[0]) / a;
  xi[1] = 2 * lambda[0] * lambda[1] + (1 + e[0] + e[1]) / a;
  xi[2] = 2 * lambda[0] * lambda[2] - e[0] / a;
  xi[3] = 2 * lambda[0] * lambda[3] - e[1] / a;
  xi[4] = lambda[1] * lambda[1] - (G[0] + 2 * G[1] + G[3]) / (4 * a);
  xi[5] = lambda[2] * lambda[2] - G[0] / (4 * a);
  xi[6] = lambda[3] * lambda[3] - G[3] / (4 * a);
  xi[7] = 2 * lambda[1] * lambda[2] + (G[0] + G[1]) / (2 * a);
  xi[8] = 2 * lambda[1] * lambda[3] + (G[1] + G[3]) / (2 * a);
  xi[9] = 2 * lambda[2] * lambda[3] - G[1] / (2 * a);

  w1Td = invW[0] * d[0] + invW[2] * d[1];
  w2Td = invW[1] * d[0] + invW[3] * d[1];
  /*row 0*/
  Lambda[0] = -w1Td * lambda[0] - (invW[0] * delta[0] + invW[2] * delta[1]) / 2;
  Lambda[3] = -w1Td * lambda[1] + (bs_pos[2].y - bs_pos[1].y) / (2 * detW);
  Lambda[6] = -w1Td * lambda[2] - (bs_pos[2].y - bs_pos[0].y) / (2 * detW);
  Lambda[9] = -w1Td * lambda[3] + (bs_pos[1].y - bs_pos[0].y) / (2 * detW);
  /*row 1*/
  Lambda[1] = -w2Td * lambda[0] - (invW[1] * delta[0] + invW[3] * delta[1]) / 2;
  Lambda[4] = -w2Td * lambda[1] + (bs_pos[1].x - bs_pos[2].x) / (2 * detW);
  Lambda[7] = -w2Td * lambda[2] + (bs_pos[2].x - bs_pos[0].x) / (2 * detW);
  Lambda[10] = -w2Td * lambda[3] - (bs_pos[1].x - bs_pos[0].x) / (2 * detW);
  /*row 2*/
  Lambda[2] = lambda[0];
  Lambda[5] = lambda[1];
  Lambda[8] = lambda[2];
  Lambda[11] = lambda[3];

  mu[0] = -w1Td;
  mu[1] = -w2Td;
  mu[2] = 1;

  mu[3] = w1Td;
  mu[4] = w2Td;
  mu[5] = -1;

  if(4 == bs_pos_len) {
    if(false == PosHelper::test_noncoplanar(bs_pos, bs_pos_len)) {
      //it_warning("4 noncoplanar BSs are needed");
      return false;
    }
    bs_pos_ = bs_pos[3];
    bs_pos_len_ = 1;
  }

  return true;
}

bool Spherical::get_pos(Point *ms_pos, const double *range, unsigned int range_len)
{
  double u[4];
  double scalar;
  int i;
  int j;
  double *coord = NULL;

  if(3 > range_len) {
    it_warning("at least 3 measurements are needed");
    return false;
  }
  if(0 == bs_pos_len_) {
    it_warning("geo_spheric_setup needs to be called first");
    return false;
  }

  scalar = comp_quadratic(u, range);
  if(scalar < 0) {
    it_warning("square root from negative number");
    scalar = std::sqrt(-scalar);
  }
  else {
    scalar = std::sqrt(scalar);
  }

  /* first solution */
  coord = (double*)ms_pos;
  for(j = 0; j < 3; ++j) {
    coord[j] = mu[j] * scalar;
    for(i = 0; i < 4; ++i) {
      coord[j] += Lambda[3 * i + j] * u[i];
    }
  }

  if(0 != scalar) {
    double r1[2];
    Point tmp_pos;

    /* second solution */
    coord = (double*)(&tmp_pos);
    for(j = 0; j < 3; ++j) {
      coord[j] = mu[3 + j] * scalar;
      for(i = 0; i < 4; ++i) {
        coord[j] += Lambda[3 * i + j] * u[i];
      }
    }
    /*smallest distance wrt the fourth measure gives the position*/
    r1[0] = fabs(PosHelper::get_dist(ms_pos, &bs_pos_) - range[3]);
    r1[1] = fabs(PosHelper::get_dist(&tmp_pos, &bs_pos_) - range[3]);
    if(r1[0] > r1[1]) {
      *ms_pos = tmp_pos;
    }
  }

  return true;
}

void Spherical::get_meas(double *meas, const unsigned int *bs_subset_idx, unsigned int len)
{
  unsigned int i;
  for(i = 0; i < len; ++i) {
    meas[i] = meas_mat_[bs_subset_idx[i]];
  }
}

bool Spherical::get_meas_mult_mat(const unsigned int *bs_pos_idx, unsigned int rows, unsigned int cols)
{
  memset(mult_mat_, 0, nb_bs_ * sizeof(*mult_mat_));
  if(1 == rows) {
    it_warning("nothing to do");
    return true;
  }
  for(unsigned int i = 0; i < (rows * cols); ++i) {
    if(nb_bs_ <= bs_pos_idx[i]) {
      return false;
    }
    ++mult_mat_[bs_pos_idx[i]];
  }
  return true;
}

void Spherical::get_meas_mult(unsigned int *meas_mult, const unsigned int *subset_idx, unsigned int meas_mult_len)
{
  for(unsigned int i = 0; i < meas_mult_len; ++i) {
    meas_mult[i] = mult_mat_[subset_idx[i]];
  }
}

bool Spherical::get_grad(double *grad, const Point *bs_pos, unsigned int bs_nb, const Point *ms_pos)
{
  double *coord_ms = NULL;
  double *coord_bs = NULL;
  double denom = 0;
  unsigned int n;
  int i;

  coord_ms = (double*)ms_pos;
  for(n = 0; n < bs_nb; ++n) {
    denom = PosHelper::get_dist(ms_pos, bs_pos + n);
    if(0 == denom) {
      it_warning("division by zero");
      return false;
    }
    coord_bs = (double*)(bs_pos + n); /*BSn*/
    for(i = 0; i < 3; ++i) {
      grad[i + 3 * n] = (coord_ms[i] - coord_bs[i]) / denom;
    }
  }
  return true;
}

class Hyperbolic : public Algorithm
{
public:
  explicit Hyperbolic(unsigned int nb_bs);
  virtual ~Hyperbolic();
  virtual Algorithm::GeoPos validate(const Point *bs_pos, unsigned int bs_pos_len, const double *r);
  virtual bool setup(const Point *bs_pos, unsigned int bs_pos_len);
  virtual bool get_pos(Point *ms_pos, const double *range, unsigned int range_len);
  virtual void set_meas_mat(const double *meas) {
    memcpy(meas_mat_, meas, nb_bs_ * nb_bs_ * sizeof(*meas));
  }
  virtual void get_meas(double *meas, const unsigned int *bs_subset_idx, unsigned int len);
  virtual bool get_meas_mult_mat(const unsigned int *bs_pos_idx, unsigned int rows, unsigned int cols);
  virtual void get_meas_mult(unsigned int *meas_mult, const unsigned int *subset_idx, unsigned int meas_mult_len);
  virtual bool get_grad(double *grad, const Point *bs_pos, unsigned int bs_nb, const Point *ms_pos);
  virtual double* get_meas_mat() const {
    return meas_mat_;
  }
  virtual unsigned int* get_mult_mat() const {
    return mult_mat_;
  }
private:
  double dot_prod(const double *v, const double *w, int len) {
    double out = v[0] * w[0];
    int i;
    for(i = 1; i < len; ++i) {
      out += v[i] * w[i];
    }
    return out;
  }
  char comp_quadratic(double *c, double *r1, double *beta, const double *r) {
    double diff_range2[3];
    double alpha[3];
    double d, e, f;
    char out = 0;
    int i;
    double delta;
    double *coord;

    for(i = 0; i < 3; ++i) {
      diff_range2[i] = r[i] * r[i];
    }
    coord = (double*)(bs_pos_);
    for(i = 0; i < 3; ++i) {
      c[i] = b[i] - 0.5 * dot_prod(A + 3 * i, diff_range2, 3);
      alpha[i] = coord[i] - c[i];
      beta[i] = dot_prod(A + 3 * i, r, 3);
    }

    /* quadratic equation */
    d = -1.0;
    e = f = 0;
    for(i = 0; i < 3; ++i) {
      d += beta[i] * beta[i];
      e += 2 * alpha[i] * beta[i];
      f += alpha[i] * alpha[i];
    }
    delta = e * e - 4 * d * f;
    if(0 > delta) {
      //it_warning("square root of negative number, inverting sign");
      delta = -delta;
      out = -1;
    }
    else if(0 < delta) {
      out = 1;
    }
    /* compute both solutions */
    r1[0] = (-e + std::sqrt(delta)) / (2 * d);
    r1[1] = (-e - std::sqrt(delta)) / (2 * d);

    return out;
  }
  //internal data
  double b[3];
  double A[9];
  Point bs_pos_[2];/*BS 1 and 5*/
  unsigned int bs_pos_len_;

  double *meas_mat_;
  unsigned int *mult_mat_;
  unsigned int nb_bs_;
};

Hyperbolic::Hyperbolic(unsigned int nb_bs) :
  nb_bs_(nb_bs),
  meas_mat_(new double[nb_bs*nb_bs]),
  mult_mat_(new unsigned int[nb_bs*nb_bs])
{
}

Hyperbolic::~Hyperbolic()
{
  delete[] meas_mat_;
  delete[] mult_mat_;
}

bool Hyperbolic::setup(const Point *bs_pos, unsigned int bs_pos_len)
{
  double K[3];
  double K0;
  double tmp[3 * 3];
  int i;

  if((4 > bs_pos_len) || (NULL == bs_pos)) {
    it_warning("at least 4 BSs are needed");
    return false;
  }

  memset(b, 0, sizeof(b));
  memset(A, 0, sizeof(A));
  memset(bs_pos_, 0, sizeof(bs_pos_));
  bs_pos_len_ = 0;

  for(i = 0; i < 3; ++i) {
    tmp[i] = bs_pos[i + 1].x - bs_pos[0].x;
    tmp[i + 3] = bs_pos[i + 1].y - bs_pos[0].y;
    tmp[i + 2 * 3] = bs_pos[i + 1].z - bs_pos[0].z;
  }
  if(false == PosHelper::inv_3by3(A, tmp)) {
    //it_warning("base stations cannot be all in the same plane");
    return false;
  }
  /*transpose as A needs to be stored row-wise*/
  for(i = 0; i < 3; ++i) {
    tmp[3 * i] = A[i];
    tmp[3 * i + 1] = A[3 + i];
    tmp[3 * i + 2] = A[6 + i];
  }
  memcpy(A, tmp, sizeof(A));

  K0 = bs_pos[0].x * bs_pos[0].x + bs_pos[0].y * bs_pos[0].y + bs_pos[0].z * bs_pos[0].z;
  for(i = 0; i < 3; ++i) {
    K[i] = bs_pos[i + 1].x * bs_pos[i + 1].x + bs_pos[i + 1].y * bs_pos[i + 1].y + bs_pos[i + 1].z * bs_pos[i + 1].z - K0;
  }
  for(i = 0; i < 3; ++i) {
    b[i] = 0.5 * dot_prod(A + 3 * i, K, 3);
  }

  memcpy(bs_pos_, bs_pos, sizeof(bs_pos_[0]));
  bs_pos_len_ = 1;
  if(5 == bs_pos_len) {
    /*test that all BSs are unique*/
    /*TODO: condition on 5th BS position*/
    if(false == PosHelper::test_uniqueness(bs_pos, bs_pos_len, 1e-3)) {
      return false;
    }
    /*copy the 5th BS position*/
    memcpy(bs_pos_ + 1, bs_pos + 4, sizeof(bs_pos_[1]));
    ++bs_pos_len_;
  }

  return true;
}

Algorithm::GeoPos Hyperbolic::validate(const Point *bs_pos, unsigned int bs_pos_len, const double *r)
{
  double c[3];
  double r1[2];
  double beta[3];
  Algorithm::GeoPos out = GEO_POS_NO;
  char delta_sign;

  if((NULL == bs_pos) || (0 == bs_pos_len) || (NULL == r)) {
    it_warning("invalid input");
    return out;
  }

  if(true == setup(bs_pos, bs_pos_len)) {

    /*quadratic equation*/
    delta_sign = comp_quadratic(c, r1, beta, r);

    /*decide if a precise location is possible*/
    if(0 == delta_sign) {
      /* very uncommon situation */
      out = GEO_POS_EXACT;
    }
    else if((0 < delta_sign) && (2 == bs_pos_len_)) {
      /*compute the determinant of the 4x4 matrix*/
      double tmp[4 * 4];
      int i;
      double det;
      for(i = 0; i < 4; ++i) {
        tmp[i] = bs_pos[i + 1].x - bs_pos[0].x;
        tmp[4 + i] = bs_pos[i + 1].y - bs_pos[0].y;
        tmp[2 * 4 + i] = bs_pos[i + 1].z - bs_pos[0].z;
        tmp[3 * 4 + i] = r[i];
      }
      det = PosHelper::det_4by4(tmp);
      out = (0 != det) ? GEO_POS_YES : GEO_POS_NO;
    }
  }
  return out;
}

bool Hyperbolic::get_pos(Point *ms_pos, const double *r, unsigned int r_len)
{
  int i, j;
  double c[3];
  double r1[2];
  double beta[3];
  char delta_sign;
  double *coord;

  if((3 > r_len) || (NULL == r)) {
    it_warning("at least 3 measurements are needed");
    return false;
  }
  if(NULL == ms_pos) {
    it_warning("output is NULL");
    return false;
  }
  if(0 == bs_pos_len_) {
    it_warning("geo_hyper_setup needs to be called first");
    return false;
  }

  /* quadratic equation */
  delta_sign = comp_quadratic(c, r1, beta, r);

  /* bs_position */
  if(0 == delta_sign) {
    /* one solution */
    coord = (double*)(ms_pos);
    for(j = 0; j < 3; ++j) {
      coord[j] = c[j] - beta[j] * r1[0];
    }
  }
  else {
    Point out[2];/* two possible solutions */

    if((4 != r_len) || (2 != bs_pos_len_)) {
      it_warning("4 measurements from 5 BSs are needed");
      return false;
    }
    /* compute both positions */
    for(i = 0; i < 2; ++i) {
      coord = (double*)(out + i);
      for(j = 0; j < 3; ++j) {
        coord[j] = c[j] - beta[j] * r1[i];
      }
    }
    /* compute TDOA distance for each solution */
    for(i = 0; i < 2; ++i) {
      r1[i] = PosHelper::get_dist(bs_pos_ + 1, out + i) - PosHelper::get_dist(bs_pos_, out + i);
      r1[i] = fabs(r1[i] - r[3]);
    }
    /*smallest distance wrt the fourth measure gives the position*/
    /*TODO: should put a condition on the 5th BS*/
    j = (r1[0] > r1[1]);
    memcpy(ms_pos, out + j, sizeof(*ms_pos));
  }
  return true;
}

void Hyperbolic::get_meas(double *meas, const unsigned int *bs_subset_idx, unsigned int len)
{
  for(unsigned int i = 0; i < (len - 1); ++i) {
    meas[i] = meas_mat_[bs_subset_idx[i + 1] + nb_bs_ * bs_subset_idx[0]];
  }
}

bool Hyperbolic::get_meas_mult_mat(const unsigned int *bs_pos_idx, unsigned int rows, unsigned int cols)
{
  unsigned int i, j;

  memset(mult_mat_, 0, nb_bs_ * nb_bs_ * sizeof(*mult_mat_));
  if(1 == rows) {
    it_warning("nothing to do");
    return true;
  }
  for(j = 0; j < cols; ++j) {
    if(nb_bs_ <= bs_pos_idx[rows * j]) {
      return false;
    }
    for(i = 0; i < (rows - 1); ++i) {
      if(nb_bs_ <= bs_pos_idx[i + 1 + rows * j]) {
        return false;
      }
      ++mult_mat_[bs_pos_idx[i + 1 + rows * j] + nb_bs_ * bs_pos_idx[rows * j]];
    }
  }
  return true;
}

void Hyperbolic::get_meas_mult(unsigned int *meas_mult, const unsigned int *subset_idx, unsigned int meas_mult_len)
{
  for(unsigned int i = 0; i < meas_mult_len; ++i) {
    meas_mult[i] = mult_mat_[subset_idx[i + 1] + nb_bs_ * subset_idx[0]];
  }
}

/* get gradient for the measurement function generated by the current set of BSs and the MS
 * By convention, the first BS is the reference BS
 * Output
 * grad: matrix 3x(bs_nb-1) stored column-wise*/
bool Hyperbolic::get_grad(double *grad, const Point *bs_pos, unsigned int bs_nb, const Point *ms_pos)
{
  double dr[3];
  double *coord_ms = NULL;
  double *coord_bs = NULL;
  double denom = 0;
  int i;
  unsigned int n;

  denom = PosHelper::get_dist(ms_pos, bs_pos);
  if(0 == denom) {
    it_warning("division by zero");
    return false;
  }
  coord_ms = (double*)ms_pos;
  coord_bs = (double*)bs_pos;/*BS0*/
  for(i = 0; i < 3; ++i) {
    dr[i] = (coord_ms[i] - coord_bs[i]) / denom;
  }
  for(n = 1; n < bs_nb; ++n) {
    denom = PosHelper::get_dist(ms_pos, bs_pos + n);
    if(0 == denom) {
      it_warning("division by zero");
      return false;
    }
    coord_bs = (double*)(bs_pos + n); /*BSn*/
    for(i = 0; i < 3; ++i) {
      grad[i + 3 * (n - 1)] = (coord_ms[i] - coord_bs[i]) / denom - dr[i];
    }
  }
  return true;
}

bool Multilateration::set_method(const itpp::bvec &method)
{
  int method_len = length(method);
  int n = 0;

  type_ = MULTI_FAILURE;
  method_.set_size(0);

  if((0 == nb_bs_) || (4 > method_len)) {
    it_warning("BSs positions are not set or too small method length");
    return false;
  }

  method_ = method;

  for(n = 0; n < method_len; ++n) {
    if(bin(1) == method[n]) {
      break;
    }
  }
  if(n == method_len) {
    type_ = MULTI_SPHERICAL;
    goto set_algo;
  }

  for(n = 0; n < method_len; ++n) {
    if(bin(0) == method[n]) {
      break;
    }
  }
  if(n == method_len) {
    type_ = MULTI_HYPERBOLIC;
    goto set_algo;
  }

  type_ = MULTI_HYBRID;

set_algo:
  //set algorithm
  if(NULL != algo_) {
    delete algo_;
  }
  switch(type_) {
  case MULTI_HYPERBOLIC:
    if((method_len + 1) != nb_bs_) {
      it_warning("For hyperbolic multilateration the number of BSs should exceed by one the number of measures");
      return false;
    }
    algo_ = new Hyperbolic(nb_bs_);
    break;
  case MULTI_HYBRID:
    if((method_len + 1) != nb_bs_) {
      it_warning("For hybrid multilateration the number of BSs should exceed by one the number of measures");
      return false;
    }
    algo_ = new Spherical(nb_bs_ - 1); /* after conversion the number of BSs is reduced by one */
    break;
  case MULTI_SPHERICAL:
    if(method_len != nb_bs_) {
      it_warning("For spherical multilateration the number of BSs should equal the number of measures");
      return false;
    }
    algo_ = new Spherical(nb_bs_);
    break;
  default:
    it_warning("unknown multilateration method");
    algo_ = NULL;
  }
  nb_fails_part = nb_fails_pos = 0;

  return (NULL != algo_) ? true : false;
}

bool Multilateration::set_bs_pos(const itpp::mat &bs_pos)
{
  int rows = bs_pos.rows();
  int cols = bs_pos.cols();

  if(((3 != cols) && (3 != rows)) || (cols == rows)) {
    it_warning("BS positions should be specified in 3D cartezian coordinates on either columns or rows");
    return false;
  }
  nb_bs_ = (3 == cols) ? rows : cols;
  if(NULL != bs_pos_) {
    delete[] bs_pos_;
  }
  bs_pos_ = new Point[nb_bs_];
  unsigned int n = 0;
  if(3 == cols) {
    for(n = 0; n < nb_bs_; ++n) {
      bs_pos_[n].x = bs_pos(n, 0);
      bs_pos_[n].y = bs_pos(n, 1);
      bs_pos_[n].z = bs_pos(n, 2);
    }
  }
  else if(3 == rows) {
    for(n = 0; n < nb_bs_; ++n) {
      bs_pos_[n].x = bs_pos(0, n);
      bs_pos_[n].y = bs_pos(1, n);
      bs_pos_[n].z = bs_pos(2, n);
    }
  }
  return true;
}

Multilateration::~Multilateration()
{
  delete algo_;
  delete bs_pos_;
}

bool Multilateration::get_pos(itpp::vec &ms_pos, const double *measures)
{
  unsigned int *subsets_idx = NULL;
  unsigned int subsets_nb = 0;
  unsigned int subset_len = 0;
  Point *tmp_bs_pos = NULL;
  Point *bs_pos = bs_pos_;
  unsigned int nb_bs = nb_bs_;
  bool out = false;
  const int method_len = length(method_);
  Point tmp_ms_pos;

  if((0 == method_len) || (NULL == bs_pos_) || (0 == nb_bs_) || (NULL == algo_)) {
    return false;
  }

  algo_->set_meas_mat(measures);

  switch(type_) {
  case MULTI_HYPERBOLIC:
    subset_len = 5;
    break;
  case MULTI_HYBRID:
    tmp_bs_pos = (Point*)malloc(nb_bs_ * sizeof(*tmp_bs_pos));
    memcpy(tmp_bs_pos, bs_pos_, nb_bs_ * sizeof(*tmp_bs_pos));
    if(false == hybrid2spherical(tmp_bs_pos, algo_->get_meas_mat())) {
      goto exit;
    }
    bs_pos = tmp_bs_pos;/* BS position might be changed after conversion */
    /* from this point on we deal with spherical multilateration (the number of BSs is reduced by one) */
    --nb_bs;
    /* fall through */
  case MULTI_SPHERICAL:
    subset_len = 4;
    break;
  default:
    (void)0;
  }

  /*ML algorithm*/
  if((false == partition(&subsets_idx, &subsets_nb, bs_pos, nb_bs, subset_len)) || (0 == subsets_nb)) {
    ++nb_fails_part;
  }
  if(0 != subsets_nb) {
    if(false == algo_->get_meas_mult_mat(subsets_idx, subset_len, subsets_nb)) {
      goto exit;
    }
    if(false == get_ml_pos(&tmp_ms_pos, bs_pos, nb_bs, subsets_idx, subsets_nb, subset_len)) {
      ++nb_fails_pos;
    }
    else {
      ms_pos.set_size(3);
      ms_pos(0) = tmp_ms_pos.x;
      ms_pos(1) = tmp_ms_pos.y;
      ms_pos(2) = tmp_ms_pos.z;
    }
  }

  out = (0 != subsets_nb) ? true : false;
exit:

  /* free subsets_idx */
  if(NULL != subsets_idx) {
    free(subsets_idx);
    subsets_idx = NULL;
  }
  /* free temporary array of BS positions */
  if(NULL != tmp_bs_pos) {
    free(tmp_bs_pos);
    tmp_bs_pos = NULL;
  }

  return out;
}

bool Multilateration::hybrid2spherical(Point *bs_pos, double *meas)
{
  int n;
  int k;
  Point *out_bs_pos = NULL;
  double *out_meas = NULL;
  int method_len = length(method_);
  static const bin zero = bin(0);

  if((4 > method_len) || (NULL == bs_pos) || (NULL == meas)) {
    return false;
  }

  if(MULTI_SPHERICAL == type_) {
    return true;/* no conversion is needed */
  }

  /* find first TOA measure */
  for(n = 0; n < method_len; ++n) {
    if(zero == method_[n]) {
      break;
    }
  }
  if(n == method_len) {
    return false;/* no TOA measure (must be hyperbolic multilateration) */
  }

  out_bs_pos = (Point*)malloc(method_len * sizeof(*out_bs_pos));
  if(NULL == out_bs_pos) {
    return false;
  }
  out_meas = (double*)malloc(method_len * sizeof(*out_meas));
  if(NULL == out_meas) {
    free(out_bs_pos);
    return false;
  }

  if(0 == n) {
    /* BS_0 is used for TOA */
    out_meas[0] = meas[0];
    out_bs_pos[0] = bs_pos[0];
    for(k = 1; k < method_len; ++k) {
      if(zero == method_[k]) {
        out_bs_pos[k] = bs_pos[k];
        out_meas[k] = meas[k];
      }
      else {
        out_meas[k] = meas[k] + out_meas[0];
        out_bs_pos[k] = bs_pos[k + 1];
      }
    }

  }
  else {
    /* BS_n is used for TOA */
    int idx = 1;
    out_meas[0] = meas[n] - meas[n - 1];
    out_bs_pos[0] = bs_pos[0];
    for(k = 0; k < method_len; ++k) {
      if((n - 1) == k) {
        continue; /* skip this measure */
      }
      if(zero == method_[k]) {
        out_bs_pos[idx] = bs_pos[k];
        out_meas[idx] = meas[k];
      }
      else {
        out_bs_pos[idx] = bs_pos[k + 1];
        out_meas[idx] = meas[k] + out_meas[0];
      }
      ++idx;
    }
  }

  /* copy to output */
  for(k = 0; k < method_len; ++k) {
    bs_pos[k] = out_bs_pos[k];
    meas[k] = out_meas[k];
  }

  free(out_bs_pos);
  free(out_meas);
  return true;
}

bool Multilateration::partition(unsigned int **subsets_idx, unsigned int *subsets_nb, const Point *bs_pos, unsigned int nb_bs, unsigned int subset_len)
{
  unsigned int nb_avail_bs = nb_bs;
  unsigned int *comb_mat = NULL;
  unsigned int nb_comb = 0;
  unsigned int i = 0;
  unsigned int n = 0;
  unsigned int k = 0;
  Point bs_pos_subset[5];/*enough for both hyperbolic and spherical trilateration*/
  unsigned int *bs_subset_idx = NULL;
  double meas[4];
  char subset_found = 0;
  int valid_res = Algorithm::GEO_POS_NO;
  bool res = true;

  if((nb_avail_bs < subset_len) || (NULL == subsets_idx) || (NULL == subsets_nb) || (NULL == bs_pos)) {
    it_warning("invalid input");
    return false;
  }
  *subsets_idx = NULL;/* make sure that the output is freed properly */
  *subsets_nb = 0;

  /* init subsets */
  IndexSet unused_subset_idx(nb_avail_bs);
  for(n = 0; n < nb_avail_bs; ++n) {
    unused_subset_idx[n] = n;
  }
  unused_subset_idx.set_size(nb_avail_bs);
  IndexSet used_subset_idx(nb_avail_bs);//size is zero
  nb_comb = PosHelper::comb_nb(nb_avail_bs - 1, subset_len - 1);
  comb_mat = (unsigned int*)malloc(subset_len * nb_comb * sizeof(*comb_mat));
  if(NULL == comb_mat) {
    it_warning("cannot allocate memory");
    return false;
  }
  bs_subset_idx = (unsigned int*)malloc(subset_len * sizeof(*bs_subset_idx));
  if(NULL == bs_subset_idx) {
    free(comb_mat);
    it_warning("cannot allocate memory");
    return false;
  }

  /*main loop*/
  subset_found = 0;
  while((0 != unused_subset_idx.get_size()) && (nb_avail_bs >= subset_len)) {
    /* check the unused set */
    if(false == PosHelper::combination(comb_mat, &nb_comb, subset_len, unused_subset_idx.get_size())) {
      //it_warning("error in PosHelper::combination");
      /* no error here */
      goto combine_both;
    }
    n = 0;
    while(n < nb_comb) {
      for(i = 0; i < subset_len; ++i) {
        bs_subset_idx[i] = unused_subset_idx[comb_mat[i + n * subset_len]];
        bs_pos_subset[i] = bs_pos[bs_subset_idx[i]];
      }
      algo_->get_meas(meas, bs_subset_idx, subset_len);
      valid_res = algo_->validate(bs_pos_subset, subset_len, meas);
      if((Algorithm::GEO_POS_EXACT == valid_res) || (Algorithm::GEO_POS_YES == valid_res)) {
        /*store the subset*/
        if(false == PosHelper::store_subset(subsets_idx, bs_subset_idx, subset_len, *subsets_nb)) {
          it_warning("error in geo_pos_store_subset");
          res = false;
          break;
        }
        /*update used and unused sets*/
        if(false == PosHelper::update_subsets(used_subset_idx, unused_subset_idx, bs_subset_idx, subset_len)) {
          it_warning("error in geo_pos_update_subsets");
          res = false;
          break;
        }
        /*increment index*/
        ++(*subsets_nb);
        /* next iteration */
        n = 0;
        if(false == PosHelper::combination(comb_mat, &nb_comb, subset_len, unused_subset_idx.get_size())) {
          /* no error code here, continue the algorithm */
          break;
        }
        continue;/*restart to iterate through new combinations set*/
      }
      ++n;
    }

    if((0 == unused_subset_idx.get_size()) || (false == res)) {
      break;
    }

    /* combine unused and used sets */
  combine_both:
    if(false == PosHelper::combination(comb_mat, &nb_comb, subset_len, nb_avail_bs)) {
      it_warning("error in geo_pos_combination");
      res = false;/* should not fail at this point */
      break;
    }
    for(n = 0; n < nb_comb; ++n) {
      if(comb_mat[subset_len - 1 + subset_len * n] < unused_subset_idx[0]) { /* TODO: check for correctness */
        continue;/*skip already searched sets*/
      }
      for(i = 0; i < subset_len; ++i) {
        k = comb_mat[i + n * subset_len];
        if(k < unused_subset_idx.get_size()) {
          bs_subset_idx[i] = unused_subset_idx[k];
        }
        else {
          /* ensure that recently added BSs are used first */
          bs_subset_idx[i] = used_subset_idx[used_subset_idx.get_size() - 1 - (k - unused_subset_idx.get_size())];
        }
        bs_pos_subset[i] = bs_pos[bs_subset_idx[i]];
      }
      algo_->get_meas(meas, bs_subset_idx, subset_len);
      valid_res = algo_->validate(bs_pos_subset, subset_len, meas);
      if((Algorithm::GEO_POS_EXACT == valid_res) || (Algorithm::GEO_POS_YES == valid_res)) {
        /*store the subset*/
        if(false == PosHelper::store_subset(subsets_idx, bs_subset_idx, subset_len, *subsets_nb)) {
          it_warning("error in geo_pos_store_subset");
          res = false;
          break;
        }
        /*update used and unused sets*/
        if(false == PosHelper::update_subsets(used_subset_idx, unused_subset_idx, bs_subset_idx, subset_len)) {
          it_warning("error in geo_pos_update_subsets");
          res = false;
          break;
        }
        /*increment index*/
        ++(*subsets_nb);
        /* stop after finding a subset */
        subset_found = 1;
        break;
      }
    }

    if((0 == unused_subset_idx.get_size()) || (false == res)) {
      break;
    }

    /* remove first BS if no subset found */
    if(0 == subset_found) {
      res = false;
      //it_warning("no subset found");
      if(false == unused_subset_idx.remove(unused_subset_idx[0])) {
        it_warning("error in utils_set_remove");
        break;
      }
      --nb_avail_bs;
    }
    else {
      subset_found = 0;
    }
  }
  free(comb_mat);
  free(bs_subset_idx);
  return res;
}

bool Multilateration::get_bs_pos_subset(Point *bs_pos_subset, const Point *bs_pos, unsigned int nb_bs, const unsigned int *subset_idx, unsigned int subset_len)
{
  unsigned int i;
  unsigned int k;
  for(i = 0; i < subset_len; ++i) {
    k = subset_idx[i];
    if(nb_bs <= k) {
      it_warning("index out of range");
      return false;
    }
    bs_pos_subset[i] = bs_pos[k];
  }
  return true;
}

bool Multilateration::get_ml_pos(Point *ms_pos, const Point *bs_pos, unsigned int nb_bs, const unsigned int *subsets_idx, unsigned int subsets_nb, unsigned int subset_len)
{
  Point bs_pos_subset[5];
  double meas[4];
  unsigned int meas_mult[4];
  double grad[3 * 4];
  double nom[3];
  double denom[9];
  double coeff[9];
  unsigned int n;
  int k;
  int i;
  double *coord = NULL;
  bool res = false;

  if((NULL == bs_pos) || (0 == nb_bs) || (NULL == subsets_idx) || (0 == subsets_nb) || (0 == subset_len) || (NULL == ms_pos)) {
    it_warning("invalid input");
    return false;
  }

  if(1 != subsets_nb) {
    memset(nom, 0, sizeof(nom));
    memset(denom, 0, sizeof(denom));
  }

  for(n = 0; n < subsets_nb; ++n) {
    /* compute MS position from current subset */
    if(false == get_bs_pos_subset(bs_pos_subset, bs_pos, nb_bs, subsets_idx + n * subset_len, subset_len)) {
      goto exit;
    }
    if(false == algo_->setup(bs_pos_subset, subset_len)) {
      it_warning("error in geo_hyper_setup");
      goto exit;
    }
    algo_->get_meas(meas, subsets_idx + n * subset_len, subset_len);
    if(false == algo_->get_pos(ms_pos, meas, 4)) {
      it_warning("error in geo_hyper_get_pos");
      goto exit;
    }
    if(1 == subsets_nb) {
      res = true;
      goto exit;
    }
    /* coefficient estimation */
    if(false == algo_->get_grad(grad, bs_pos_subset, subset_len, ms_pos)) {
      it_warning("error in geo_hyper_get_grad");
      goto exit;
    }
    algo_->get_meas_mult(meas_mult, subsets_idx + n * subset_len, 4);
    if(false == prod(coeff, grad, meas_mult, 3, 4)) {
      it_warning("error in geo_hyper_get_grad");
      goto exit;
    }
    /*nominator and denominator for the ML estimator*/
    coord = (double*)ms_pos;
    for(k = 0; k < 3; ++k) {
      for(i = 0; i < 3; ++i) {
        nom[k] += coeff[k + 3 * i] * coord[i];
        denom[k + 3 * i] += coeff[k + 3 * i];
      }
    }
  }
  /*ML estimator from several subsets*/
  if(true == PosHelper::inv_3by3(coeff, denom)) {
    coord = (double*)ms_pos;
    for(n = 0; n < 3; ++n) {
      coord[n] = 0;
      for(k = 0; k < 3; ++k) {
        coord[n] += coeff[n + 3 * k] * nom[k];
      }
    }
    res = true;
  }

exit:
  return res;
}

bool Multilateration::prod(double *out, const double *AT, const unsigned int *d, unsigned int cols, unsigned int rows)
{
  unsigned int i;
  unsigned int n;
  unsigned int k;

  if((NULL == out) || (NULL == AT) || (NULL == d) || (0 == cols) || (0 == rows)) {
    it_warning("invalid input");
    return false;
  }

  for(i = 0; i < cols; ++i) {  /* each line first input */
    for(n = 0; n < cols; ++n) {  /* each column of the second input */
      out[i + cols * n] = 0.0;
      for(k = 0; k < rows; ++k) {  /* each column of the first input or each row of the second input */
        if(0 == d[k]) {
          it_warning("division by zero");
          return false;
        }
        out[i + cols * n] += (AT[i + k * cols] / (double)d[k]) * AT[n + k * cols];
      }
    }
  }
  return true;
}

double Multilateration::get_crlb(const vec &ms_pos, double sigma2)
{
  if((3 != length(ms_pos)) || (0 == nb_bs_)) {
    it_error("invalid input");
  }
  if(0.0 == sigma2) {
    return 0.0;
  }

  vec pos(3);
  vec pos_ref(3);
  pos_ref(0) = bs_pos_[0].x;
  pos_ref(1) = bs_pos_[0].y;
  pos_ref(2) = bs_pos_[0].z;

  static const bin zero = bin(0);
  const unsigned int method_len = length(method_);
  mat dh(3, method_len);
  dh.zeros();
  unsigned int n;
  unsigned int i;
  for(n = 0; n < method_len; ++n) {
    if(zero == method_(n)) {
      pos(0) = bs_pos_[n].x;
      pos(1) = bs_pos_[n].y;
      pos(2) = bs_pos_[n].z;
      for(i = 0; i < 3; ++i) {
        dh(i, n) = (ms_pos(i) - pos(i)) / std::sqrt(sum_sqr(ms_pos - pos));
      }
    }
    else {
      pos(0) = bs_pos_[n + 1].x;
      pos(1) = bs_pos_[n + 1].y;
      pos(2) = bs_pos_[n + 1].z;
      for(i = 0; i < 3; ++i) {
        dh(i, n) = (ms_pos(i) - pos(i)) / std::sqrt(sum_sqr(ms_pos - pos)) -
                   (ms_pos(i) - pos_ref(i)) / std::sqrt(sum_sqr(ms_pos - pos_ref));
      }
    }
  }

  mat info_mat(3, 3);
  for(n = 0; n < 3; ++n) {
    for(i = 0; i < 3; ++i) {
      info_mat(n, i) = dh.get_row(n) * dh.get_row(i);
    }
  }
  info_mat = info_mat / sigma2;
  return std::sqrt(sum(diag(inv(info_mat))));
}

}
