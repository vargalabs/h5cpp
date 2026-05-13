/*!
 * \file
 * \brief IT file endianness test program
 * \author Andy Panov and Adam Piatyszek
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
#include <iomanip>

using namespace itpp;
using namespace std;

// To rewrite the ITFILE_TEST_FILE uncomment the following definition
// #define SAVE_DATA

const string test_file_name = ITFILE_TEST_FILE;

//reference data;
const char c_ref = 'c';
const bool bo_ref = true;
const bin b_ref = 0;
const short s_ref = 1234;
const int i_ref = -1234567890;
const float f_ref = -12345.6f;
const double d_ref = 2.1e-8;
const complex<double> cd_ref = std::complex<double>(1.0, -1.5);
const string st_ref = "abcdefghij 0123456789";

const bvec bv_ref = "0 1 0 1 1";
const ivec iv_ref = "2 5 -400 2 -10";
const vec v_ref = "1e-9 0.2 0.7 -1.0 0.0";
const cvec cv_ref = "(0,2) (1.5,7.2)";

const bmat bm_ref = "0 1 0; 1 1 1";
const imat im_ref = "2 5; 4 -10; 0 3";
const mat m_ref = "1e-9 0.2 0.7; -1.0 0.0 3e10";
const cmat cm_ref = "(0,2) (-1.5,7.2); (1.1,2) (7,-4e-5)";

const Array<bvec> abv_ref = "{[0 1] [0 0 0] [1 0] [1 1 0]}";
const Array<ivec> aiv_ref = "{[1 2 3] [4 5 6] [7 8 9]}";
const Array<vec> av_ref = "{[1 2e4 -3.0] [-0.5 6e-9] [7e-3]}";
const Array<cvec> acv_ref = "{[(0,2) (0.5,-0.5)] [(2,1)] [(0,0) (-4.6,2)]}";

const Array<bmat> abm_ref = "{[0 1 0; 0 0 1] [1 1; 0 0; 1 0] [1; 1; 0]}";
const Array<imat> aim_ref = "{[0 2 3; 0 -1 9; 2 3 -1] [1 10 100] [0 4; -2 3]}";
const Array<mat> am_ref = "{[0.5 2e7; 0.5 -0.5] [1e-4 3 4; 0.1 0.2 .3]}";
const Array<cmat> acm_ref = "{[(0,2) (0.5,-0.5); (2,1) (0,0)] "
                            "[(1.1,2) (7,-4e-5); (0,2) (1.5,7.2)]}";

const string contents_list_ref[] = {
  "   c        int8      1   char variable",
  "  bo        bool      1   bool variable",
  "   b         bin      1   bin variable",
  "   s       int16      2   short int variable",
  "   i       int32      4   int variable",
  "   f     float32      4   ",
  "   d     float64      8   ",
  "  cd    cfloat64     16   ",
  "  st      string     29   ",
  "  bv        bvec     13   ",
  "  iv        ivec     28   ",
  "   v        dvec     48   ",
  "  cv       dcvec     40   ",
  "  bm        bmat     22   ",
  "  im        imat     40   ",
  "   m        dmat     64   ",
  "  cm       dcmat     80   ",
  " abv   bvecArray     50   ",
  " aiv   ivecArray     68   ",
  "  av    vecArray     80   ",
  " acv   cvecArray    112   ",
  " abm   bmatArray     71   ",
  " aim   imatArray    120   ",
  "  am    matArray    120   ",
  " acm   cmatArray    168   "
};
namespace itpp{
//comparison of two arrays
  template<typename T>
  bool operator==(const itpp::Array<T> l, const itpp::Array<T> r)
  {
    if (l.length() != r.length()) return false;
    for(int i = 0; i < l.length(); ++i)
      if(l(i) != r(i)) return false;
    return true;
  }
}

//checker and item list extraction helpers
template <typename T> inline
void check_equality(const T& val, const T& ref)
{
  ASSERT_EQ(val, ref) << "equality failure for " << typeid(T).name() << endl;
}

int list_file_contents(it_ifile& f, ostream& str)
{
  std::string name, type, desc;
  uint64_t size;
  int n = 0;
  while (f.seek(n++)) {
    f.info(name, type, desc, size);
    str << setw(4) << name << setw(12) << type << setw(7) << size
         << "   " << desc << endl;
  }
  return n - 1;
}

//base class for ItFile tests
class ItFileTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
#ifdef SAVE_DATA
    it_file fw;
    fw.open(test_file_name, true);
    fw << Name("c", "char variable") << c_ref;
    fw << Name("bo", "bool variable") << bo_ref;
    fw << Name("b", "bin variable") << b_ref;
    fw << Name("s", "short int variable") << s_ref;
    fw << Name("i", "int variable") << i_ref;
    fw << Name("f") << f_ref;
    fw << Name("d") << d_ref;
    fw << Name("cd") << cd_ref;
    fw << Name("st") << st_ref;
    fw << Name("bv") << bv_ref;
    fw << Name("iv") << concat(iv_ref, iv_ref, iv_ref);
    fw << Name("v") << v_ref;
    fw << Name("cv") << cv_ref;
    fw << Name("bm") << bm_ref;
    fw << Name("im") << im_ref;
    fw << Name("m") << m_ref;
    fw << Name("cm") << cm_ref;
    fw << Name("abv") << abv_ref;
    fw << Name("aiv") << aiv_ref;
    fw << Name("av") << av_ref;
    fw << Name("acv") << acv_ref;
    fw << Name("abm") << abm_ref;
    fw << Name("aim") << aim_ref;
    fw << Name("am") << am_ref;
    fw << Name("acm") << acm_ref;
    fw.remove("iv");
    fw << Name("iv") << iv_ref;
    fw.close();
#endif
  }

};


//it_file contents
TEST_F(ItFileTest, Contents)
{
  stringstream str;
  it_ifile ff(test_file_name);
  int items_read = list_file_contents(ff,str);
  string cur; int i = 0;
  while(getline(str,cur))
  {
    ASSERT_EQ(cur, contents_list_ref[i])<< "item info mismatch for item #" << i << endl;
    ++i;
  }
}

//it_file extraction test
TEST_F(ItFileTest, Extraction)
{

  it_ifile ff(test_file_name);

  char c;
  bool bo;
  bin b;
  short s;
  int i;
  float f;
  double d;
  complex<double> cd;
  string st;

  bvec bv;
  ivec iv;
  vec v;
  cvec cv;

  bmat bm;
  imat im;
  mat m;
  cmat cm;

  Array<bvec> abv;
  Array<ivec> aiv;
  Array<vec> av;
  Array<cvec> acv;

  Array<bmat> abm;
  Array<imat> aim;
  Array<mat> am;
  Array<cmat> acm;

  ff >> Name("abm") >> abm;
  check_equality(abm,abm_ref);
  ff >> Name("abv") >> abv;
  check_equality(abv,abv_ref);
  ff >> Name("acm") >> acm;
  check_equality(acm,acm_ref);
  ff >> Name("acv") >> acv;
  check_equality(acv,acv_ref);
  ff >> Name("aim") >> aim;
  check_equality(aim,aim_ref);
  ff >> Name("aiv") >> aiv;
  check_equality(aiv,aiv_ref);
  ff >> Name("am") >> am;
  check_equality(am,am_ref);
  ff >> Name("av") >> av;
  check_equality(av,av_ref);
  ff >> Name("b") >> b;
  check_equality(b,b_ref);
  ff >> Name("bm") >> bm;
  check_equality(bm,bm_ref);
  ff >> Name("bo") >> bo;
  check_equality(bo,bo_ref);
  ff >> Name("bv") >> bv;
  check_equality(bv,bv_ref);
  ff >> Name("c") >> c;
  check_equality(c,c_ref);
  ff >> Name("cd") >> cd;
  check_equality(cd,cd_ref);
  ff >> Name("cm") >> cm;
  check_equality(cm,cm_ref);
  ff >> Name("cv") >> cv;
  check_equality(cv,cv_ref);
  ff >> Name("d") >> d;
  check_equality(d,d_ref);
  ff >> Name("f") >> f;
  check_equality(f,f_ref);
  ff >> Name("i") >> i;
  check_equality(i,i_ref);
  ff >> Name("im") >> im;
  check_equality(im,im_ref);
  ff >> Name("iv") >> iv;
  check_equality(iv,iv_ref);
  ff >> Name("m") >> m;
  check_equality(m,m_ref);
  ff >> Name("s") >> s;
  check_equality(s,s_ref);
  ff >> Name("st") >> st;
  check_equality(st,st_ref);
  ff >> Name("v") >> v;
  check_equality(v,v_ref);
  ff.close();
}

//this test verifies deletion and packing of it_file items
TEST_F(ItFileTest, Packing)
{
  ivec iv0 = "0 0";
  ivec iv1 = ones_i(100);
  ivec iv2 = "2 2 2 2";
  ivec iv3 = "3";

  it_file ff1("itfile_test_pack.it", true);
  ff1 << Name("iv0") << iv0 << flush;
  ff1 << Name("iv1") << iv1 << flush;
  ff1 << Name("iv2") << iv2 << flush;
  ff1.remove("iv1");
  ff1 << Name("iv1") << ivec("1") << flush;
  ff1 << Name("iv3") << iv3 << flush;
  ff1 << Name("iv4") << iv3 << flush;
  ff1.remove("iv3");
  ff1.low_level().seekg(0, std::ios::end);
  streamsize before_packing = ff1.low_level().tellg();
  ff1.pack();
  //ensure that file was squeezed by the previous pack operation
  ff1.low_level().seekg(0, std::ios::end);
  ASSERT_LT(ff1.low_level().tellg(), before_packing) << "size squeeze failure after packing" << std::endl;
  ff1.close();
}

