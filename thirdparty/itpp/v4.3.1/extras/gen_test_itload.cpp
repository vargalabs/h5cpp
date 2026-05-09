/*!
 * \file
 * \brief Generates test_pyitpp.it itfile used for unit tests by test_pyitpp.py script.
 * \author Bogdan Cristea
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;

int main()
{
  it_file ff("test_pyitpp.it");

  //generate scalars
  {
    bin a = 1;
    ff << Name("a") << a;
  }

  {
    char b = '2';
    ff << Name("b") << b;
  }

  {
    short c = -3;
    ff << Name("c") << c;
  }

  {
    int d = -4;
    ff << Name("d") << d;
  }

  {
    float e = 5.67;
    ff << Name("e") << e;
  }

  {
    double f = 8.91234567;
    ff << Name("f") << f;
  }

  {
    complex<float> g(-1.234, 4.321);
    ff << Name("g") << g;
  }

  {
    complex<double> h(-1.234567, 4.321654);
    ff << Name("h") << h;
  }

  //generate vectors
  {
    bvec i = "0 1 1 0 0 1";
    ff << Name("i") << i;
  }

  {
    string j = "abc";
    ff << Name("j") << j;
  }

  {
    svec k = "10:19";
    ff << Name("k") << k;
  }

  {
    ivec l = "20:29";
    ff << Name("l") << l;
  }

  //vec of floats not supported

  {
    vec m = "30:1:39";
    ff << Name("m") << m;
  }

  //cvec of floats not supported

  {
    cvec n(10);
    complex<double> step(0.5, -0.5);
    n[0] = complex<double>(0.0, 0.0);
    for(unsigned int i = 1; i < 10; ++i) {
      n[i] = n[i-1] + step;
    }
    ff << Name("n") << n;
  }

  //generate matrices
  {
    bmat o = "0 1 1; 0 0 1";
    ff << Name("o") << o;
  }

  //no string

  {
    smat p = "1 2 3; 4 5 6";
    ff << Name("p") << p;
  }

  {
    imat q = "11 12 13; 14 15 16";
    ff << Name("q") << q;
  }

  //mat of floats not supported

  {
    mat r = "1.5 1.6 1.7; 2.3 2.4 2.5";
    ff << Name("r") << r;
  }

  //cmat of floats not supported

  {
    cmat s(3, 2);
    complex<double> step(0.5, -0.5);
    s(0) = complex<double>(0.0, 0.0);
    for(unsigned int i = 1; i < 6; ++i) {
      s(i) = s(i - 1) + step;
    }
    ff << Name("s") << s;
  }

  //simple arrays
  {
    Array<bin> t("{0 1 0 1 1 0 0 0 1}");
    ff << Name("t") << t;
  }

  {
    Array<short> u("{0 1 2 3 4 5 6 7 8}");
    ff << Name("u") << u;
  }

  {
    Array<int> v("{10 11 12 13 14 15 16 17 18}");
    ff << Name("v") << v;
  }

  {
    Array<float> w("{1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8}");
    ff << Name("w") << w;
  }

  {
    Array<double> x("{1.2222 2.3333 3.44444 4.55555}");
    ff << Name("x") << x;
  }

  {
    Array<complex<float> > y(10);
    complex<float> step(0.5, -0.5);
    y(0) = complex<float>(0, 0);
    for(int i = 1; i < 10; ++i) {
      y(i) = y(i - 1) + step;
    }
    ff << Name("y") << y;
  }

  {
    Array<complex<double> > z(10);
    complex<double> step(-0.5, 0.5);
    z(0) = complex<double>(0, 0);
    for(int i = 1; i < 10; ++i) {
      z(i) = z(i - 1) + step;
    }
    ff << Name("z") << z;
  }

  //arrays of vectors
  {
    Array<bvec> aa("{[1 1 1 0 0 0 1], [0 0 1 1], [1 0 1]}");
    ff << Name("aa") << aa;
  }

  {
    Array<svec> bb("{[1 1 1 2 3 0 1], [0 6 1 1], [1 7 1]}");
    ff << Name("bb") << bb;
  }

  {
    Array<ivec> cc("{[10 10 10 20 30 0 10], [0 60 10 10], [10 70 10]}");
    ff << Name("cc") << cc;
  }

  {
    Array<vec> dd("{[0.1 0.1 0.1 0.2 0.3 0.0 0.1], [0.0 0.6 0.1 0.1], [0.1 0.7 0.1]}");
    ff << Name("dd") << dd;
  }


  {
    Array<cvec> ee(3);
    cvec tmp(5);
    tmp(0) = complex<double>(0.0, 0.0);
    complex<double> step(0.5, -0.5);
    for(int i = 1; i < 5; ++i) {
      tmp(i) = tmp(i - 1) + step;
    }
    ee(0) = tmp;
    ee(1) = tmp.get(0, 2);
    ee(2) = tmp.get(1, 4);
    ff << Name("ee") << ee;
  }

  {
    Array<string> gg(3);
    gg(0) = "abcd";
    gg(1) = "abc";
    gg(2) = "defghijk";
    ff << Name("gg") << gg;
  }

  //arrays of matrices
  {
    Array<bmat> hh("{[1 1; 1 0; 0 1], [0 0 1; 1 0 1], [1; 0; 1]}");
    ff << Name("hh") << hh;
  }

  {
    Array<smat> ii("{[1 2; 1 3; 4 7], [4 7 1; 8 6 1], [3; 8; 5]}");
    ff << Name("ii") << ii;
  }

  {
    Array<imat> jj("{[11 21; 11 31; 41 71], [41 71 11; 81 61 11], [31; 81; 51]}");
    ff << Name("jj") << jj;
  }

  {
    Array<mat> kk("{[1.1 2.1; 1.1 3.1; 4.1 7.1], [4.1 7.1 1.1; 8.1 6.1 1.1], [3.1; 8.1; 5.1]}");
    ff << Name("kk") << kk;
  }

  {
    Array<cmat> ll(3);
    cmat s(3, 2);
    complex<double> step(0.5, -0.5);
    s(0) = complex<double>(0.0, 0.0);
    for(unsigned int i = 1; i < 6; ++i) {
      s(i) = s(i - 1) + step;
    }
    ll(0) = s;
    ll(1) = s(0, 1, 0, 1);
    ll(2) = s(0, 0, 0, 1);
    ff << Name("ll") << ll;
  }

  ff.close();

  cout << "done" << endl;
}
