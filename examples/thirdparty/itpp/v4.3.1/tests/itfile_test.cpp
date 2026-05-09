/*!
 * \file
 * \brief IT file endianness test program
 * \author Adam Piatyszek
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
#include <iomanip>

using namespace itpp;
using namespace std;

// To run extensive tests uncomment the following definition
// #define EXTENSIVE_TESTS
// To rewrite the ITFILE_TEST_FILE uncomment the following definition
// #define SAVE_DATA

#if (!defined(ITFILE_TEST_FILE) || defined(ITPP_ITFILE_EXCLUDED))

int main()
{
  cerr << "ITFILE_TEST_FILE not defined or ITPP_ITFILE_EXCLUDED defined. Test skipped." << endl;
  return 1;
}

#else

int main()
{
  char c, c_ref = 'c';
  bool bo, bo_ref = true;
  bin b, b_ref = 0;
  short s, s_ref = 1234;
  int i, i_ref = -1234567890;
  float f, f_ref = -12345.6f;
  double d, d_ref = 2.1e-8;
  complex<double> cd, cd_ref = std::complex<double>(1.0, -1.5);
  string st, st_ref = "abcdefghij 0123456789";

  bvec bv, bv_ref = "0 1 0 1 1";
  ivec iv, iv_ref = "2 5 -400 2 -10";
  vec v, v_ref = "1e-9 0.2 0.7 -1.0 0.0";
  cvec cv, cv_ref = "(0,2) (1.5,7.2)";

  bmat bm, bm_ref = "0 1 0; 1 1 1";
  imat im, im_ref = "2 5; 4 -10; 0 3";
  mat m, m_ref = "1e-9 0.2 0.7; -1.0 0.0 3e10";
  cmat cm, cm_ref = "(0,2) (-1.5,7.2); (1.1,2) (7,-4e-5)";

  Array<bvec> abv, abv_ref = "{[0 1] [0 0 0] [1 0] [1 1 0]}";
  Array<ivec> aiv, aiv_ref = "{[1 2 3] [4 5 6] [7 8 9]}";
  Array<vec> av, av_ref = "{[1 2e4 -3.0] [-0.5 6e-9] [7e-3]}";
  Array<cvec> acv, acv_ref = "{[(0,2) (0.5,-0.5)] [(2,1)] [(0,0) (-4.6,2)]}";

  Array<bmat> abm, abm_ref = "{[0 1 0; 0 0 1] [1 1; 0 0; 1 0] [1; 1; 0]}";
  Array<imat> aim, aim_ref = "{[0 2 3; 0 -1 9; 2 3 -1] [1 10 100] [0 4; -2 3]}";
  Array<mat> am, am_ref = "{[0.5 2e7; 0.5 -0.5] [1e-4 3 4; 0.1 0.2 .3]}";
  Array<cmat> acm, acm_ref = "{[(0,2) (0.5,-0.5); (2,1) (0,0)] "
                             "[(1.1,2) (7,-4e-5); (0,2) (1.5,7.2)]}";

#ifdef SAVE_DATA
  it_file fw;
  fw.open(string(ITFILE_TEST_FILE), true);
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
  std::string name, type, desc;
  uint64_t size;
  int n = 0;
  cout << "Name |      Type | Size | Description\n"
       << "------------------------------------------------\n";
  it_ifile ff(string(ITFILE_TEST_FILE));
  while (ff.seek(n++)) {
    ff.info(name, type, desc, size);
    cout << setw(4) << name << setw(12) << type << setw(7) << size
         << "   " << desc << endl;
  }
  cout << "------------------------------------------------\n\n";

  ff >> Name("abm") >> abm;
  ff >> Name("abv") >> abv;
  ff >> Name("acm") >> acm;
  ff >> Name("acv") >> acv;
  ff >> Name("aim") >> aim;
  ff >> Name("aiv") >> aiv;
  ff >> Name("am") >> am;
  ff >> Name("av") >> av;
  ff >> Name("b") >> b;
  ff >> Name("bm") >> bm;
  ff >> Name("bo") >> bo;
  ff >> Name("bv") >> bv;
  ff >> Name("c") >> c;
  ff >> Name("cd") >> cd;
  ff >> Name("cm") >> cm;
  ff >> Name("cv") >> cv;
  ff >> Name("d") >> d;
  ff >> Name("f") >> f;
  ff >> Name("i") >> i;
  ff >> Name("im") >> im;
  ff >> Name("iv") >> iv;
  ff >> Name("m") >> m;
  ff >> Name("s") >> s;
  ff >> Name("st") >> st;
  ff >> Name("v") >> v;
  ff.close();

  cout << "char    : '" << fixed << c << "'" << endl
       << "          '" << fixed << c_ref << "'" << endl
       << "bool    : " << fixed << bo << endl
       << "          " << fixed << bo_ref << endl
       << "bin     : " << fixed << b << endl
       << "          " << fixed << b_ref << endl
       << "short   : " << fixed << s << endl
       << "          " << fixed << s_ref << endl
       << "int     : " << fixed << i << endl
       << "          " << fixed << i_ref << endl
       << "float   : " << fixed << f << endl
       << "          " << fixed << f_ref << endl
       << "double  : " << fixed << d << endl
       << "          " << fixed << d_ref << endl
       << "complex : " << fixed << cd << endl
       << "          " << fixed << cd_ref << endl
       << "string  : \"" << fixed << st << "\"" << endl
       << "          \"" << fixed << st_ref << "\"" << endl << endl;

  cout << "bvec    : " << fixed << bv << endl
       << "          " << fixed << bv_ref << endl
       << "ivec    : " << fixed << iv << endl
       << "          " << fixed << iv_ref << endl
       << "vec     : " << fixed << v << endl
       << "          " << fixed << v_ref << endl
       << "cvec    : " << fixed << cv << endl
       << "          " << fixed << cv_ref << endl << endl;

  cout << "bmat    :\n" << fixed << bm << endl << fixed << bm_ref << endl
       << "imat    :\n" << fixed << im << endl << fixed << im_ref << endl
       << "mat     :\n" << fixed << m << endl << fixed << m_ref << endl
       << "cmat    :\n" << fixed << cm << endl << fixed << cm_ref << endl << endl;

  cout << "Array<bvec> :\n" << fixed << abv << endl << fixed << abv_ref << endl
       << "Array<ivec> :\n" << fixed << aiv << endl << fixed << aiv_ref << endl
       << "Array<vec>  :\n" << fixed << av << endl << fixed << av_ref << endl
       << "Array<cvec> :\n" << fixed << acv << endl << fixed << acv_ref << endl << endl;

  cout << "Array<bmat> :\n" << fixed << abm << endl << fixed << abm_ref << endl
       << "Array<imat> :\n" << fixed << aim << endl << fixed << aim_ref << endl
       << "Array<mat> :\n" << fixed << am << endl << fixed << am_ref << endl
       << "Array<cmat> :\n" << fixed << acm << endl << fixed << acm_ref << endl << endl;

#ifdef EXTENSIVE_TESTS
  ivec iv0 = "0 0";
  ivec iv1 = ones_i(100);
  ivec iv2 = "2 2 2 2";
  ivec iv3 = "3";

  it_file ff1("itfile_test_extensive.it", true);
  ff1 << Name("iv0") << iv0 << flush;
  ff1 << Name("iv1") << iv1 << flush;
  ff1 << Name("iv2") << iv2 << flush;
  ff1.remove("iv1");
  ff1 << Name("iv1") << ivec("1") << flush;
  ff1 << Name("iv3") << iv3 << flush;
  ff1 << Name("iv4") << iv3 << flush;
  ff1.remove("iv3");
  ff1.low_level().seekg(0, std::ios::end);
  it_info("Size before packing: " << ff1.low_level().tellg());
  ff1.pack();
  ff1.low_level().seekg(0, std::ios::end);
  it_info("Size after packing:  " << ff1.low_level().tellg());
  ff1.close();

  it_ifile ff2("itfile_test_extensive.it");
  n = 0;
  while (ff2.seek(n++)) {
    ff2.info(name, type, desc, size);
    ff2 >> iv1;
    cout << "Name = " << name << "  Type = " << type << "  Size = " << size
         << "  Desc = \"" << desc << "\"  Data = " << iv1 << endl;
  }
  ff2.close();
#endif

  return 0;
}

#endif
