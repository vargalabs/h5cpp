/*!
 * \file
 * \brief Parser test program
 * \author Pal Frenger and Adam Piatyszek
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

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

static
void assert_vec_p(const vec &ref, const vec &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol) << line;
  }
}
#define assert_vec(ref, act) assert_vec_p(ref, act, __LINE__)
static
void assert_mat_p(const mat &ref, const mat &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.rows(), act.rows()) << line;
  ASSERT_EQ(ref.cols(), act.cols()) << line;
  for (int n = 0; n < ref.rows(); ++n) {
    for (int k = 0; k < ref.cols(); ++k) {
      ASSERT_NEAR(ref(n,k), act(n,k), tol) << line;
    }
  }
}
#define assert_mat(ref, act) assert_mat_p(ref, act, __LINE__)
static
void assert_cvec_p(const cvec &ref, const cvec &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n).real(), act(n).real(), tol) << line;
    ASSERT_NEAR(ref(n).imag(), act(n).imag(), tol) << line;
  }
}
#define assert_cvec(ref, act) assert_cvec_p(ref, act, __LINE__)

TEST(Parser, all)
{
#ifndef PARSER_TEST_FILE
  ASSERT_TRUE(false) << "PARSER_TEST_FILE not defined. Test skipped.";
#else
  const double tol = 1e-12;
  Parser p;
  int a;
  double b;
  string c;
  ivec d;
  vec e;
  svec f;
  bvec g;
  imat h;
  mat i;
  smat j;
  bmat k;
  bool l, m;
  string n;
  Array<Array<cvec> > o, q;

  p.set_silentmode(true);

  // Use the Parser class on a parameter file
  p.init(string(PARSER_TEST_FILE));
  a  = p.get_int("a");
  ASSERT_EQ(4, a);
  b = p.get_double("b");   //The default value of b
  ASSERT_NEAR(2.35, b, tol);
  b = p.get_double("b", 1); //The first alternative value of b
  ASSERT_NEAR(2.36, b, tol);
  b = p.get_double("b", 2); //The second alternative value of b
  ASSERT_NEAR(2.37, b, tol);
  c  = p.get_string("c");
  ASSERT_TRUE(string("Hello World") == c);
  d  = p.get_ivec("d");
  ivec ref_i = "1 2 3 4 3 2 -1";
  ASSERT_TRUE(ref_i == d);
  e  = p.get_vec("e");
  vec ref = "1.2 2 3.33 4.01 3.2 2 -1.2";
  assert_vec(ref, e);
  f  = p.get_svec("f");
  svec ref_s = "1 2 3 4 3 2 -1";
  ASSERT_TRUE(ref_s == f);
  g  = p.get_bvec("g");
  bvec ref_b = "0 1 0 0 1";
  ASSERT_TRUE(ref_b == g);
  h  = p.get_imat("h");
  imat ref_imat = "1 2 3; 4 5 6";
  ASSERT_TRUE(ref_imat == h);
  i  = p.get_mat("i");
  mat ref_mat = "1 2; -3 4.2";
  assert_mat(ref_mat, i);
  j  = p.get_smat("j");
  smat ref_smat = "1 -2 -3; 4 -5 6";
  ASSERT_TRUE(ref_smat == j);
  k  = p.get_bmat("k");
  bmat ref_bmat = "0 1 0 1; 1 0 1 0; 0 0 0 0; 1 1 1 1";
  ASSERT_TRUE(ref_bmat == k);
  l  = p.get_bool("l");
  ASSERT_FALSE(l);
  m  = p.get_bool("m");
  ASSERT_TRUE(m);
  n  = p.get_string("n");
  ASSERT_TRUE(string("Hello World") == n);

  // Use the Parser class on a char pointer (usually the command line input)
  char argv1_0[] = "a=345";
  char argv1_1[] = "b=1.001";
  char argv1_2[] = "c=\"Hello Bird\"";
  char argv1_3[] = "d=[1,2,3,-1,-2,-3]";
  int  argc1 = 4;
  char *argv1[4] = { argv1_0, argv1_1, argv1_2, argv1_3 };
  p.init(argc1, argv1);
  a = p.get_int("a");
  ASSERT_EQ(345, a);
  b = p.get_double("b");
  ASSERT_DOUBLE_EQ(1.001, b);
  c = p.get_string("c");
  ASSERT_TRUE(string("Hello Bird") == c);
  d = p.get_ivec("d");
  ref_i = "1 2 3 -1 -2 -3";
  ASSERT_TRUE(ref_i == d);

  // Use the Parser class on a parameter file and a char pointer
  // The data in the char pointer are selected first. The data in the
  // file are used as default values
  char argv2_0[] = "a=345";
  char argv2_1[] = "d=[1,-2,3,-4,5,-6,7,-8]";
  int argc2 = 2;
  char *argv2[2] = { argv2_0, argv2_1 };
  p.init(string(PARSER_TEST_FILE), argc2, argv2);
  a = p.get_int("a");
  ASSERT_EQ(345, a);
  b = p.get_double("b");
  ASSERT_DOUBLE_EQ(2.35, b);
  c = p.get_string("c");
  ASSERT_TRUE(string("Hello World") == c);
  d = p.get_ivec("d");
  ref_i = "1 -2 3 -4 5 -6 7 -8";
  ASSERT_TRUE(ref_i == d);

  // Use the Parser class on an Array of strings
  Array<string> parser_data(4);
  parser_data(0) = "a=-11";
  parser_data(1) = "b=3.14";
  parser_data(2) = "c=\"Hello Nerd\"";
  parser_data(3) = "d=[0,1,2,0,3,4,7]";
  p.init(parser_data);
  a = p.get_int("a");
  ASSERT_EQ(-11, a);
  b = p.get_double("b");
  ASSERT_DOUBLE_EQ(3.14, b);
  c = p.get_string("c");
  ASSERT_TRUE(string("Hello Nerd") == c);
  d = p.get_ivec("d");
  ref_i = "0 1 2 0 3 4 7";
  ASSERT_TRUE(ref_i == d);

  // Use the Parser::get() method on a parameter file
  p.init(string(PARSER_TEST_FILE));
  ASSERT_TRUE(p.get(a, "a"));
  ASSERT_TRUE(p.get(b, "b"));
  ASSERT_TRUE(p.get(b, "b", 1));
  ASSERT_TRUE(p.get(b, "b", 2));
  ASSERT_TRUE(p.get(c, "c"));
  ASSERT_TRUE(p.get(d, "d"));
  ASSERT_TRUE(p.get(e, "e"));
  ASSERT_TRUE(p.get(f, "f"));
  ASSERT_TRUE(p.get(g, "g"));
  ASSERT_TRUE(p.get(h, "h"));
  ASSERT_TRUE(p.get(i, "i"));
  ASSERT_TRUE(p.get(j, "j"));
  ASSERT_TRUE(p.get(k, "k"));
  ASSERT_TRUE(p.get(l, "l"));
  ASSERT_TRUE(p.get(m, "m"));
  ASSERT_TRUE(p.get(n, "n"));
  ASSERT_TRUE(p.get(o, "o"));
  // The following tries can not be found in the data file, so their initial
  // values should not be changed:
  set_array(q, "{{[1][2 3][4 5 6]}}");
  ASSERT_FALSE(p.get(q, "q"));
  cvec ref_c = "1";
  assert_cvec(ref_c, q(0)(0));
  ref_c = "2 3";
  assert_cvec(ref_c, q(0)(1));
  ref_c = "4 5 6";
  assert_cvec(ref_c, q(0)(2));
  bool r = true;
  ASSERT_FALSE(p.get(r, "r"));
  ASSERT_TRUE(r);
  int s = 7;
  ASSERT_FALSE(p.get(s, "s"));
  ASSERT_EQ(7, s);
  double t = 3.14;
  ASSERT_FALSE(p.get(t, "t"));
  ASSERT_DOUBLE_EQ(3.14, t);
  std::string u = "test string";
  p.get(u, "u");
  ASSERT_TRUE(string("test string") == u);

  // Check if variables exist in the Parser
  ASSERT_TRUE(p.exist("a"));
  ASSERT_FALSE(p.exist("aa"));
#endif   /* ----- #ifndef PARSER_TEST_FILE  ----- */
}
