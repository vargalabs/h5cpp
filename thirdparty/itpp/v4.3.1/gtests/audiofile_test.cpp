/*!
 * \file
 * \brief Audio streaming test program
 * \author Andy Panov
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2013  (see AUTHORS file for a list of contributors)
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
#include <itpp/itsrccode.h>
#include "gtest/gtest.h"
#include <iomanip>

using namespace itpp;
using namespace std;

// To rewrite the AUDIOFILE_REF_FILE uncomment the following definition
// #define SAVE_DATA

//reference data
const string ref_file_name = AUDIOFILE_REF_FILE;
const int ref_num_samples = 5000;
const int ref_num_channels = 2;
mat ref_channels(ref_num_samples,ref_num_channels);
std::string ref_annotation("ITPP reference audio file contains 1004 Hz tone at 8 kHz sampling rate.");
const Audio_Encoding ref_encoding = enc_double;
const int ref_sampling_rate = 8000;
const int ref_tone_frequency = 1004;
Audio_Stream_Description ref_description(ref_encoding,ref_sampling_rate,ref_num_channels);

static
const double tol = 1e-9;

//checkers
void check_sample(double s, double ref)
{
  ASSERT_NEAR(ref,s, tol) << "sample test failed." << endl;
}

void check_samples(const vec& s, const vec& ref)
{
  ASSERT_EQ(s.length(),ref.length()) <<"check_samples: vector length test failed" << endl;
  for(int j = 0; j < s.length(); ++j){
    ASSERT_NEAR(ref(j), s(j), tol) <<"check_samples: failed for sample number: "<< j << endl;
  }
}

void check_samples(const mat& s, const mat& ref)
{
  ASSERT_EQ(s.rows(),ref.rows()) <<"check_samples: matrix rows test failed" << endl;
  ASSERT_EQ(s.cols(),ref.cols()) <<"check_samples: matrix columns test failed" << endl;
  for(int i = 0; i < s.rows(); ++i)
    for(int j = 0; j < s.cols(); ++j){
      ASSERT_NEAR(ref(i,j), s(i,j), tol) <<"check_samples: failed for row: "<< i << " column: " << j << endl;
    }
}

void check_description(const Audio_Stream_Description d, const Audio_Stream_Description ref)
{
  ASSERT_TRUE(is_valid(d)) << "check_description: invalid description" << endl;
  ASSERT_TRUE(is_valid(ref)) << "check_description: invalid reference" << endl;
  ASSERT_EQ(d.get_encoding(), ref.get_encoding()) << "check_description: encodings mismatch" << endl;
  ASSERT_EQ(d.get_sampling_rate(), ref.get_sampling_rate()) << "check_description: sampling rates mismatch" << endl;
  ASSERT_EQ(d.get_num_channels(), ref.get_num_channels()) << "check_description: number of channels mismatch" << endl;
  ASSERT_EQ(d.get_description(), ref.get_description()) << "check_description: annotations mismatch" << endl;
}

//base class for AudioFile tests
class AudioFileTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
  //generate reference tone
  for(int i = 0; i < ref_num_samples; ++i){
    ref_channels(i,0) = Audio_Sample<ref_encoding>::decode(Audio_Sample<ref_encoding>::encode(sin(2*pi*ref_tone_frequency/ref_sampling_rate*i)));
    ref_channels(i,1) = Audio_Sample<ref_encoding>::decode(Audio_Sample<ref_encoding>::encode(cos(2*pi*ref_tone_frequency/ref_sampling_rate*i)));
  }
  ref_description.set_description(ref_annotation.c_str());
#ifdef SAVE_DATA
  Audio_Stream_Description d(ref_encoding,ref_sampling_rate,ref_num_channels);
  d.set_description(ref_annotation.c_str());
  SND_Out_File ref(ref_file_name.c_str(),d);
  ref.write(ref_channels);
#endif
  }
};


//read test with reference file (refence file should be audible with .au files player (e.g. audacity))
TEST_F(AudioFileTest, In)
{
  SND_In_File ref(ref_file_name.c_str());
  check_description(ref.get_description(),ref_description);
  ASSERT_EQ(ref.num_samples(),ref_num_samples) << "Unexpected number of samples" <<endl;
  //read all data from stream
  mat m = ref.read((int)ref.num_samples());
  check_samples(m,ref_channels);
  //test positioning and sample-wise reads - read sample numbwer 500 from channel 0;
  ASSERT_TRUE(ref.seek_read(500)) << "Seek failed" << endl;
  double s;
  ASSERT_TRUE(ref.read_sample(s,0)) << "Read failed" << endl;
  check_sample(s,ref_channels(500,0));
  //test channel-wise reads - read 200 samples from channel 1 starting from sample 501;
  vec v = ref.read_channel(200,1);
  check_samples(v,ref_channels.get_col(1)(501,700));
}

//write test
TEST_F(AudioFileTest, Out)
{
  SND_Out_File f_out("itpp_test.au",ref_description);
  //write single sample
  ASSERT_TRUE(f_out.write_sample(1.0,0)) << "Write failed" << endl;
  ASSERT_TRUE(f_out.seek_write(0)) << "Seek failed" << endl;
  ASSERT_TRUE(f_out.write_sample(0.5,1)) << "Write failed" << endl;
  ASSERT_EQ(f_out.num_samples(),1)<< "Unexpected number of samples" <<endl;
  //write vector
  ASSERT_TRUE(f_out.seek_write(0)) << "Seek failed" << endl;
  vec ch0 = zeros(100);
  vec ch1 = 0.5*ones(100);
  ASSERT_TRUE(f_out.write_channel(ch0,0)) << "Write failed" << endl;
  ASSERT_TRUE(f_out.seek_write(0)) << "Seek failed" << endl;
  ASSERT_TRUE(f_out.write_channel(ch1,1)) << "Write failed" << endl;
  ASSERT_EQ(f_out.num_samples(),100)<< "Unexpected number of samples" <<endl;
  //write matrix
  ASSERT_TRUE(f_out.seek_write(0)) << "Seek failed" << endl;
  ASSERT_TRUE(f_out.write(ref_channels)) << "Write failed" << endl;
  ASSERT_EQ(f_out.num_samples(),ref_num_samples)<< "Unexpected number of samples" <<endl;
}

//test interleaved i/o with the audio stream
TEST_F(AudioFileTest, InOut)
{
  SND_IO_File f("itpp_test.au");
  check_description(f.get_description(),ref_description);
  ASSERT_EQ(f.num_samples(),ref_num_samples) << "Unexpected number of samples" <<endl;
  //read all data from stream
  mat m = f.read((int)f.num_samples());
  check_samples(m,ref_channels);
  //write and read single sample from position 10 of channel 1
  double s_ref = 0.5, s;
  ASSERT_TRUE(f.seek_write(10)) << "Seek failed" << endl;
  ASSERT_TRUE(f.write_sample(s_ref,1)) << "Write failed" << endl;
  ASSERT_TRUE(f.seek_read(10)) << "Seek failed" << endl;
  ASSERT_TRUE(f.read_sample(s,1)) << "Read failed" << endl;
  check_sample(s,s_ref);

  //write and read 100 samples from position 50 of channel 0
  vec v_ref = zeros(100);
  ASSERT_TRUE(f.seek_write(50)) << "Seek failed" << endl;
  ASSERT_TRUE(f.write_channel(v_ref,1)) << "Write failed" << endl;
  ASSERT_TRUE(f.seek_read(50)) << "Seek failed" << endl;
  vec v = f.read_channel(100,1);
  check_samples(v,v_ref);

  //write from and read to the matrix (position 150, length 10)
  mat m_ref = ref_channels.get(150,159,0,1);
  ASSERT_TRUE(f.seek_write(150)) << "Seek failed" << endl;
  ASSERT_TRUE(f.write(m_ref)) << "Write failed" << endl;
  ASSERT_TRUE(f.seek_read(150)) << "Seek failed" << endl;
  m = f.read(10);
  check_samples(m,m_ref);
}
