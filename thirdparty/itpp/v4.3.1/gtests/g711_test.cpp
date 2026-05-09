/*!
 * \file
 * \brief Test program for G.711 logarithmic codecs
 * \author Andy Panov
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

/*
 * This code uses repacked G.191 test vectors distributed under the
 * ITU-T SOFTWARE TOOLS' GENERAL PUBLIC LICENSE.

 * Following is the information about original authors of G.711 test suite:
 * Simao Ferraz de Campos Neto         Rudolf Hofmann
 * CPqD/Telebras                       PHILIPS KOMMUNIKATIONS INDUSTRIE AG
 * DDS/Pr.11                           Kommunikationssysteme
 * Rd. Mogi Mirim-Campinas Km.118      Thurn-und-Taxis-Strasse 14
 * 13.080-061 - Campinas - SP (Brazil) D-8500 Nuernberg 10 (Germany)

 * Phone : +55-192-39-6637             Phone : +49-911-526-2603
 * FAX   : +55-192-39-2179             FAX   : +49-911-526-3385
 * EMail : tdsimao@cpqd.ansp.br        EMail : hf@pkinbg.uucp

 * In order to get license and original test vectors you should download G.191 software tools
 * archive. License can be found in \Software\gen-lic.txt file. Original test vectors are stored in
 * /Software/stl2009/g711/g711-tst.zip
*/


#include "gtest/gtest.h"
#include <iomanip>
#include <itpp/itbase.h>
#include <itpp/itsrccode.h>

using namespace itpp;

//This test program uses original G.191 test vectors repacked to itpp file format.
//Please uncomment following define in order to update/regenerate the test vectors. Also you should
//unpack the G.191 original files (sweep.src, sweep-r.u, sweep-r.u-u, sweep-r.a, sweep-r.a-a)
//to the test executable folder.
//Original test vectors cover the whole range of 16-bit integers (2^16 samples). Please adjust the
//num_samples constant accordingly if you want to run the test on the reduced range and speed things up.

//#define UPDATE_TEST_VECTORS

const std::string test_file_name = G711_TEST_FILE;
const int num_samples = 256*256; //use full test set

//equality checker
template <typename T> inline
void check_equality(const Array<T>& left, const Array<T>& right)
{
  int size = left.size();
  ASSERT_EQ (size , right.size()) << "array size mismatch";
  for(int i = 0; i < size; ++i) {
    ASSERT_EQ(left(i) , right(i)) << "equality failure for element # " << i << std::endl;
  }
}

//base class for ItFile tests
class G711 : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
#ifdef UPDATE_TEST_VECTORS
  it_file g711_refs(test_file_name);
  Array<int16_t> samples(num_samples);
  int cnt;
  //load 16-bit sweep reference
  {
    bifstream g191_file("sweep.src"); int16_t s;
    for(cnt = 0; cnt < num_samples && g191_file; ++cnt){
      g191_file >> s ; samples(cnt) = s;
    }
    //it seems to be a bug in gtest (see issue 247, http://code.google.com/p/googletest/issues/detail?id=247)
    //gtest does not seem to mark fixture failed on failures in SetUpTestCase(), so we have to terminate here.
    it_error_if(cnt != num_samples,"Failed to read sweep.src file contents.");
    g711_refs << Name("g191_Sweep") << samples;
  }

  //ulaw-encoded g191_Sweep
  {
    bifstream g191_file("sweep-r.u"); int16_t s;
    for(cnt = 0; cnt < num_samples && g191_file; ++cnt) {
      g191_file >> s; samples(cnt) = s;
    }
    it_error_if(cnt != num_samples,"Failed to read sweep-r.u file contents.");
    g711_refs << Name("g191_ulaw") << samples;
  }

  //load ulaw-decoded samples.  Test samples are left-justified, so scale them down to fit into 13 LSBs
  {
    bifstream g191_file("sweep-r.u-u"); int16_t s;
    for(cnt = 0; cnt < num_samples && g191_file; ++cnt) {
      g191_file >> s; samples(cnt) = (s>>2);
    }

    it_error_if(cnt != num_samples, "Failed to read sweep-r.u-u file contents.");
    g711_refs << Name("g191_ulaw_dec") << samples;
  }

  //alaw-encoded g191_Sweep
  {
    bifstream g191_file("sweep-r.a"); int16_t s;
    for(cnt = 0; cnt < num_samples && g191_file; ++cnt) {
      g191_file >> s; samples(cnt) = s;
    }
    it_error_if(cnt != num_samples,"Failed to read sweep-r.a file contents.");
    g711_refs << Name("g191_alaw") << samples;
  }

  //load alaw-decoded samples. Test samples are left-justified, so scale them down to fit into 13 LSBs
  {
    bifstream g191_file("sweep-r.a-a"); int16_t s;
    for(cnt = 0; cnt < num_samples && g191_file; ++cnt) {
      g191_file >> s; samples(cnt) = (s>>3);
    }
    it_error_if(cnt != num_samples,"Failed to read sweep-r.a-a file contents.");
    g711_refs << Name("g191_alaw_dec") << samples;
  }
#endif
  }
};

//G.711 u-codec test
TEST_F(G711, uCodec)
{
  it_file ref_data(test_file_name);

  Array<int16_t> ref_in_samples(num_samples);
  Array<int16_t> ref_enc_samples(num_samples);
  Array<int16_t> ref_dec_samples(num_samples);

  ref_data >> Name("g191_Sweep")>>ref_in_samples;
  ref_data >> Name("g191_ulaw") >> ref_enc_samples;
  ref_data >> Name("g191_ulaw_dec") >> ref_dec_samples;

  Array<int16_t> enc_samples(num_samples);
  Array<int16_t> dec_samples(num_samples);

  //encode
  for(int cnt = 0; cnt < num_samples; ++cnt){
    enc_samples(cnt) = ulaw_compress(ref_in_samples(cnt)>>2);
  }
  check_equality(enc_samples, ref_enc_samples);
  //decode
  for(int cnt = 0; cnt < num_samples; ++cnt){
    dec_samples(cnt) = ulaw_expand(static_cast<uint8_t>(ref_enc_samples(cnt)));
  }
  check_equality(dec_samples, ref_dec_samples);
}

//G.711 a-codec test
TEST_F(G711, aCodec)
{
  it_file ref_data(test_file_name);

  Array<int16_t> ref_in_samples(num_samples);
  Array<int16_t> ref_enc_samples(num_samples);
  Array<int16_t> ref_dec_samples(num_samples);

  ref_data >> Name("g191_Sweep")>>ref_in_samples;
  ref_data >> Name("g191_alaw") >> ref_enc_samples;
  ref_data >> Name("g191_alaw_dec") >> ref_dec_samples;

  Array<int16_t> enc_samples(num_samples);
  Array<int16_t> dec_samples(num_samples);

  //encode
  for(int cnt = 0; cnt < num_samples; ++cnt){
    enc_samples(cnt) = alaw_compress(ref_in_samples(cnt)>>3);
  }
  check_equality(enc_samples, ref_enc_samples);
  //decode
  for(int cnt = 0; cnt < num_samples; ++cnt){
    dec_samples(cnt) = alaw_expand(static_cast<uint8_t>(ref_enc_samples(cnt)));
  }
  check_equality(dec_samples, ref_dec_samples);
}

