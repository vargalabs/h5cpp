/*!
 * \file
 * \brief Implementation of special vectors and matrices
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger, Adam Piatyszek
 *         and Erik G. Larsson
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

#include <itpp/base/specmat.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/itcompat.h>
#include <itpp/base/matfunc.h>


namespace itpp
{

ivec find(const bvec &invector)
{
    it_assert(invector.size() > 0, "find(): vector cannot be empty");
    ivec temp(invector.size());
    int pos = 0;
    for (int i = 0;i < invector.size();i++) {
        if (invector(i) == bin(1)) {
            temp(pos) = i;
            pos++;
        }
    }
    temp.set_size(pos, true);
    return temp;
}

//! \cond

#define CREATE_SET_FUNS(typef,typem,name,value) \
  typef name(int size)    \
  {      \
    typef t(size);    \
    t = value;     \
    return t;     \
  }      \
      \
    typem name(int rows, int cols)  \
    {      \
      typem t(rows, cols);   \
      t = value;    \
      return t;     \
    }

#define CREATE_EYE_FUN(type,name,zero,one) \
  type name(int size) {    \
    type t(size,size);    \
    t = zero;     \
    for (int i=0; i<size; i++)   \
      t(i,i) = one;    \
    return t;     \
  }

CREATE_SET_FUNS(vec, mat, ones, 1.0)
CREATE_SET_FUNS(bvec, bmat, ones_b, bin(1))
CREATE_SET_FUNS(ivec, imat, ones_i, 1)
CREATE_SET_FUNS(cvec, cmat, ones_c, std::complex<double>(1.0))

CREATE_SET_FUNS(vec, mat, zeros, 0.0)
CREATE_SET_FUNS(bvec, bmat, zeros_b, bin(0))
CREATE_SET_FUNS(ivec, imat, zeros_i, 0)
CREATE_SET_FUNS(cvec, cmat, zeros_c, std::complex<double>(0.0))

CREATE_EYE_FUN(mat, eye, 0.0, 1.0)
CREATE_EYE_FUN(bmat, eye_b, bin(0), bin(1))
CREATE_EYE_FUN(imat, eye_i, 0, 1)
CREATE_EYE_FUN(cmat, eye_c, std::complex<double>(0.0), std::complex<double>(1.0))

//! \endcond

vec impulse(int size)
{
    vec t(size);
    t.clear();
    t[0] = 1.0;
    return t;
}

vec linspace(double from, double to, int points)
{
    if (points < 2) {
        // This is the "Matlab definition" of linspace
        vec output(1);
        output(0) = to;
        return output;
    }
    else {
        vec output(points);
        double step = (to - from) / double(points - 1);
        int i;
        for (i = 0; i < (points-1); i++)
            output(i) = from + i * step;
        output(i) = to;
        return output;
    }
}

vec zigzag_space(double t0, double t1, int K)
{
    it_assert(K > 0, "zigzag_space:() K must be positive");
    ivec N = "0 1";

    int n = 2;
    for (int k = 0; k < K; k++) {
        ivec Nn = 2 * N;
        for (int i = 1; i < length(Nn); i += 2)  {
            Nn = concat(Nn, i);
            n++;
        }
        N = Nn;
    }

    vec T0 = linspace(t0, t1, n);
    vec Tt = zeros(n);
    for (int i = 0; i < n; i++) {
        Tt(i) = T0(N(i));
    }
    return Tt;
}

// Construct a Hadamard-imat of size "size"
imat hadamard(int size)
{
    it_assert(size > 0, "hadamard(): size is not a power of 2");
    int logsize = ceil_i(::log2(static_cast<double>(size)));
    it_assert(pow2i(logsize) == size, "hadamard(): size is not a power of 2");

    imat H(size, size);
    H(0, 0) = 1;

    for (int i = 0; i < logsize; ++i) {
        int pow2 = 1 << i;
        for (int k = 0; k < pow2; ++k) {
            for (int l = 0; l < pow2; ++l) {
                H(k, l) = H(k, l);
                H(k + pow2, l) = H(k, l);
                H(k, l + pow2) = H(k, l);
                H(k + pow2, l + pow2) = (-1) * H(k, l);
            }
        }
    }
    return H;
}

imat jacobsthal(int p)
{
    int quadratic_residue;
    imat out(p, p);
    int i, j;

    out = -1; // start with all elements equal to "-1"

    // Generate a complete list of quadratic residues
    for (i = 0; i < (p - 1) / 2; i++) {
        quadratic_residue = ((i + 1) * (i + 1)) % p;
        // set this element in all rows (col-row) = quadratic_residue
        for (j = 0; j < p; j++) {
            out(j, (j + quadratic_residue) % p) = 1;
        }
    }

    // set diagonal elements to zero
    for (i = 0; i < p; i++) {
        out(i, i) = 0;
    }
    return out;
}

imat conference(int n)
{
    it_assert_debug(n % 4 == 2, "conference(int n); wrong size");
    int pm = n - 1; // p must be odd prime, not checked
    imat out(n, n);

    out.set_submatrix(1, 1, jacobsthal(pm));
    out.set_submatrix(0, 0, 1, n - 1, 1);
    out.set_submatrix(1, n - 1, 0, 0, 1);
    out(0, 0) = 0;

    return out;
}

const cmat toeplitz(const cvec &c)
{
    int s = c.size();
    cmat output(s, s);
    cvec c_conj = conj(c);
    for (int i = 1; i < s; ++i) {
        for (int j = 0; j < s - i; ++j) {
            output(i + j, j) = c_conj(i);
        }
    }
    // start from j = 0 here, because the main diagonal is not conjugated
    for (int j = 0; j < s; ++j) {
        for (int i = 0; i < s - j; ++i) {
            output(i, i + j) = c(j);
        }
    }
    return output;
}

mat rotation_matrix(int dim, int plane1, int plane2, double angle)
{
    mat m;
    double c = std::cos(angle), s = std::sin(angle);

    it_assert(plane1 >= 0 && plane2 >= 0 &&
              plane1 < dim && plane2 < dim && plane1 != plane2,
              "Invalid arguments to rotation_matrix()");

    m.set_size(dim, dim, false);
    m = 0.0;
    for (int i = 0; i < dim; i++)
        m(i, i) = 1.0;

    m(plane1, plane1) = c;
    m(plane1, plane2) = -s;
    m(plane2, plane1) = s;
    m(plane2, plane2) = c;

    return m;
}

void house(const vec &x, vec &v, double &beta)
{
    double sigma, mu;
    int n = x.size();

    v = x;
    if (n == 1) {
        v(0) = 1.0;
        beta = 0.0;
        return;
    }
    sigma = sum(sqr(x(1, n - 1)));
    v(0) = 1.0;
    if (sigma == 0.0)
        beta = 0.0;
    else {
        mu = std::sqrt(sqr(x(0)) + sigma);
        if (x(0) <= 0.0)
            v(0) = x(0) - mu;
        else
            v(0) = -sigma / (x(0) + mu);
        beta = 2 * sqr(v(0)) / (sigma + sqr(v(0)));
        v /= v(0);
    }
}

void givens(double a, double b, double &c, double &s)
{
    double t;

    if (b == 0) {
        c = 1.0;
        s = 0.0;
    }
    else {
        if (fabs(b) > fabs(a)) {
            t = -a / b;
            s = -1.0 / std::sqrt(1 + t * t);
            c = s * t;
        }
        else {
            t = -b / a;
            c = 1.0 / std::sqrt(1 + t * t);
            s = c * t;
        }
    }
}

void givens(double a, double b, mat &m)
{
    double t, c, s;

    m.set_size(2, 2);

    if (b == 0) {
        m(0, 0) = 1.0;
        m(1, 1) = 1.0;
        m(1, 0) = 0.0;
        m(0, 1) = 0.0;
    }
    else {
        if (fabs(b) > fabs(a)) {
            t = -a / b;
            s = -1.0 / std::sqrt(1 + t * t);
            c = s * t;
        }
        else {
            t = -b / a;
            c = 1.0 / std::sqrt(1 + t * t);
            s = c * t;
        }
        m(0, 0) = c;
        m(1, 1) = c;
        m(0, 1) = s;
        m(1, 0) = -s;
    }
}

mat givens(double a, double b)
{
    mat m(2, 2);
    givens(a, b, m);
    return m;
}

void givens_t(double a, double b, mat &m)
{
    double t, c, s;

    m.set_size(2, 2);

    if (b == 0) {
        m(0, 0) = 1.0;
        m(1, 1) = 1.0;
        m(1, 0) = 0.0;
        m(0, 1) = 0.0;
    }
    else {
        if (fabs(b) > fabs(a)) {
            t = -a / b;
            s = -1.0 / std::sqrt(1 + t * t);
            c = s * t;
        }
        else {
            t = -b / a;
            c = 1.0 / std::sqrt(1 + t * t);
            s = c * t;
        }
        m(0, 0) = c;
        m(1, 1) = c;
        m(0, 1) = -s;
        m(1, 0) = s;
    }
}

mat givens_t(double a, double b)
{
    mat m(2, 2);
    givens_t(a, b, m);
    return m;
}

//! Template instantiation of eye
template void eye(int, mat &);
//! Template instantiation of eye
template void eye(int, bmat &);
//! Template instantiation of eye
template void eye(int, imat &);
//! Template instantiation of eye
template void eye(int, cmat &);

//! Template instantiation of linspace_fixed_step
template vec linspace_fixed_step(double, double, double);
//! Template instantiation of linspace_fixed_step
template ivec linspace_fixed_step(int, int, int);
//! Template instantiation of linspace_fixed_step
template svec linspace_fixed_step(short int, short int, short int);

} // namespace itpp
