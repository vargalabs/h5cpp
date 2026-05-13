/*!
 * \file
 * \brief Implementation of audio Audio classes and functions
 * \author Tobias Ringstrom, Adam Piatyszek and Andy Panov
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

#include <itpp/srccode/audiofile.h>
#include <itpp/base/converters.h>
#include <itpp/base/ittypes.h>
#include <itpp/base/itassert.h>
#include <iostream>

//! \cond
namespace itpp
{

//magic id of snd file header
static const uint32_t snd_magic = 0x2e736e64;
//maximum length of annotation to extract from snd file
static const std::size_t max_annotation_length = 1024;

//////////////////////////////////////////////////
//
// Audio_Samples_Reader - templated implementation of Audio_Samples_Reader_If
//
//////////////////////////////////////////////////
template<typename Binary_In_Stream, Audio_Encoding Encoding>
class Audio_Samples_Reader : public audiofile_details::Audio_Samples_Reader_If
{
public:
  Audio_Samples_Reader(Binary_In_Stream& str, std::streamoff start, int nc):
    _str(str), _start_pos(start), _num_channels(nc), _cur_pos(0){}
  bool read_sample(double& s, int ch);
  //Read n samples from audio channel ch
  vec read_channel(int n, int ch);
  mat read(int n);
  virtual std::streamoff tell() const;
  virtual bool seek(std::streamoff n);
  virtual std::streamoff num_samples();
private:
  static const std::size_t sample_size = Audio_Sample<Encoding>::enc_sample_size;
  typedef typename Audio_Sample<Encoding>::enc_sample_type sample_type;
  //Number of audio channels
  int _num_channels;
  //! First sample offset from the start of the file
  std::streamoff _start_pos;
  //! Current position in samples
  std::streamoff _cur_pos;
  //! Binary stream
  Binary_In_Stream& _str;
};

template<typename Binary_In_Stream, Audio_Encoding Encoding>
std::streamoff Audio_Samples_Reader<Binary_In_Stream,Encoding>::tell() const
{
  return _cur_pos;
}

template<typename Binary_In_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Reader<Binary_In_Stream,Encoding>::seek(std::streamoff n)
{
  _str.seekg(_start_pos + (_cur_pos * _num_channels *sample_size), std::ios_base::beg);
  if(_str){
    _cur_pos = n;
    return true;
  }
  else{
    return false;
  }
}

template<typename Binary_In_Stream, Audio_Encoding Encoding>
std::streamoff Audio_Samples_Reader<Binary_In_Stream,Encoding>::num_samples()
{
  _str.seekg(0, std::ios_base::end);
  if(!_str) return -1;
  std::streamoff end_pos = _str.tellg();
  return (end_pos - _start_pos)/(_num_channels * sample_size);
}

//read single channel samples starting at current position
template<typename Binary_In_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Reader<Binary_In_Stream,Encoding>::read_sample(double& s, int ch)
{
  if(ch >= _num_channels) return false;
  std::streamoff read_pos = _start_pos + (_cur_pos * _num_channels + ch )*sample_size;
  _str.seekg(read_pos, std::ios_base::beg);
  if(!_str) return false;
  sample_type raw_sample;
  _str >> raw_sample;
  if(_str){
    s = Audio_Sample<Encoding>::decode(raw_sample);
    _cur_pos++;
    return true;
  }
  return false;
}

//read n samples from channel ch starting at current position
template<typename Binary_In_Stream, Audio_Encoding Encoding>
vec Audio_Samples_Reader<Binary_In_Stream,Encoding>::read_channel(int n, int ch)
{
  //ignore threshold - ignore() is used instead of seekg()
  //if stride between samples is smaller then this threshold
  static const std::streamsize ignore_threshold = 64;

  if((n <= 0) || (ch >= _num_channels)) return vec();
  vec ret(n);

  //read first n-1 samples
  const std::streamsize stride = sample_size*(_num_channels - 1);
  _str.seekg(_start_pos + (_cur_pos *_num_channels + ch) * sample_size, std::ios_base::beg);
  for(int i = 0; (i < (n-1)) && _str; ++i) {
    sample_type raw_sample; _str >> raw_sample;
    ret(i) = Audio_Sample<Encoding>::decode(raw_sample);
    if(stride > ignore_threshold)
      _str.seekg(stride, std::ios_base::cur);
    else
      _str.ignore(stride);
  }

  //read last sample
  if(_str){
    sample_type raw_sample; _str >> raw_sample;
    ret(n-1) = Audio_Sample<Encoding>::decode(raw_sample);
  }

  if(_str){
    _cur_pos += n;
  }
  else{
    ret.set_size(0);
  }

  return ret;
}

//read n samples from all channels starting at current position
template<typename Binary_In_Stream, Audio_Encoding Encoding>
mat Audio_Samples_Reader<Binary_In_Stream,Encoding>::read(int n)
{
  if(n <= 0) return mat();
  mat ret(n,_num_channels);

  //read samples
  const std::streamsize stride = sample_size*(_num_channels - 1);
  _str.seekg(_start_pos + _cur_pos * sample_size *_num_channels, std::ios_base::beg);
  for(int i = 0; (i < n) && _str; ++i) {
    for(int j = 0; j < _num_channels && _str; ++j) {
      sample_type raw_sample; _str >> raw_sample;
      ret(i,j) = Audio_Sample<Encoding>::decode(raw_sample);
    }
  }

  if(_str){
    _cur_pos += n;
  }
  else{
    ret.set_size(0,0);
  }
  return ret;
}

//make audio samples reader to input data from stream
template<typename Binary_In_Stream>
audiofile_details::Audio_Samples_Reader_If* make_reader(Binary_In_Stream& str,
    std::streamoff start_pos, Audio_Stream_Description* d)
{
  Audio_Encoding encoding = d->get_encoding();
  int num_channels = d->get_num_channels();
  switch(encoding){
    case enc_mulaw8:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_mulaw8>(str,start_pos,num_channels);
    case enc_alaw8:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_alaw8>(str,start_pos,num_channels);
    case enc_linear8:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_linear8>(str,start_pos,num_channels);
    case enc_linear16:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_linear16>(str,start_pos,num_channels);
    case enc_linear24:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_linear24>(str,start_pos,num_channels);
    case enc_linear32:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_linear32>(str,start_pos,num_channels);
    case enc_float:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_float>(str,start_pos,num_channels);
    case enc_double:
      return new Audio_Samples_Reader<Binary_In_Stream, enc_double>(str,start_pos,num_channels);
    case enc_unknown:
    default:
      return 0;
  }
}
//////////////////////////////////////////////////
//
// Audio_Samples_Writer - templated implementation of Audio_Samples_Writer_If
//
//////////////////////////////////////////////////
template<typename Binary_Out_Stream, Audio_Encoding Encoding>
class Audio_Samples_Writer : public audiofile_details::Audio_Samples_Writer_If
{
public:
  Audio_Samples_Writer(Binary_Out_Stream& str, std::streamoff start, int nc):
    _str(str), _start_pos(start), _num_channels(nc), _cur_pos(0),
    _zero(Audio_Sample<Encoding>::encode(0.0)){}
  virtual bool write_sample(const double& s, int ch);
  virtual bool write_channel(const vec& s, int ch);
  //Write n samples to audio channel ch
  virtual bool write(const mat& s);
  virtual std::streamoff tell() const;
  virtual bool seek(std::streamoff n);
  virtual std::streamoff num_samples();
private:
  static const std::size_t sample_size = Audio_Sample<Encoding>::enc_sample_size;
  typedef typename Audio_Sample<Encoding>::enc_sample_type sample_type;
  //Number of audio channels
  int _num_channels;
  //! First sample offset from the start of the file
  std::streamoff _start_pos;
  //! Current position in samples
  std::streamoff _cur_pos;
  //! Binary stream
  Binary_Out_Stream& _str;
  //Zero sample
  sample_type _zero;
};

template<typename Binary_Out_Stream, Audio_Encoding Encoding>
std::streamoff Audio_Samples_Writer<Binary_Out_Stream,Encoding>::tell() const
{
  return _cur_pos;
}

template<typename Binary_Out_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Writer<Binary_Out_Stream,Encoding>::seek(std::streamoff n)
{
  _str.seekp(_start_pos + (n * _num_channels *sample_size), std::ios_base::beg);
  if(_str){
    _cur_pos = n;
    return true;
  }
  else{
    return false;
  }
}

template<typename Binary_Out_Stream, Audio_Encoding Encoding>
std::streamoff Audio_Samples_Writer<Binary_Out_Stream,Encoding>::num_samples()
{
  _str.seekp(0, std::ios_base::end);
  if(!_str) return -1;
  std::streamoff end_pos = _str.tellp();
  return (end_pos - _start_pos)/(_num_channels * sample_size);
}

//write single sample starting at current position
template<typename Binary_Out_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Writer<Binary_Out_Stream,Encoding>::write_sample(const double& s, int ch)
{
  if(ch >= _num_channels) return false;
  std::streamoff write_pos = _start_pos + (_cur_pos * _num_channels + ch )*sample_size;
  _str.seekp(write_pos, std::ios_base::beg);
  if(!_str) return false;
  _str << Audio_Sample<Encoding>::encode(s);
  if(_str){
    _cur_pos++;
    return true;
  }
  return false;
}

//write single channel samples starting at current position
template<typename Binary_Out_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Writer<Binary_Out_Stream,Encoding>::write_channel(const vec& s, int ch)
{
  if(ch >= _num_channels) return false;

  int len = s.length();
  //overall number of already written samples
  std::streamoff ns = num_samples();
  if(ns < 0) return false;
  //compute number of samples to overwrite and number of samples to add to the end of file
  int to_overwrite = (int)std::min(ns - _cur_pos, (std::streamoff)len);

  const std::streamoff stride = sample_size*(_num_channels - 1);

  int i = 0;
  //overwrite samples
  if(to_overwrite)
  {
    //process first to_overwrite-1 samples
    _str.seekp(_start_pos + (_cur_pos * _num_channels + ch )*sample_size, std::ios_base::beg);
    for(i = 0; (i < (to_overwrite-1)) && _str; ++i) {
      _str << Audio_Sample<Encoding>::encode(s(i));
      if(stride) _str.seekp(stride,std::ios_base::cur);
    }
    if(_str){
      _str << Audio_Sample<Encoding>::encode(s(i)); ++i;
    }
  }

  //add samples to the end of file
  if(i < len)
  {
    _str.seekp(_start_pos + ns * _num_channels * sample_size, std::ios_base::beg);
    for(; (i < len) && _str; ++i) {
      for(int j = 0; (j < _num_channels) && _str; ++j){
        if(j == ch)
          _str << Audio_Sample<Encoding>::encode(s(i));
        else
          _str << _zero;
      }
    }
  }

  if(_str){
    _cur_pos += len;
    return true;
  }
  return false;
}

//write samples starting at current position
template<typename Binary_Out_Stream, Audio_Encoding Encoding>
bool Audio_Samples_Writer<Binary_Out_Stream,Encoding>::write(const mat& s)
{
  if(s.cols() < _num_channels) return false;
  int len = s.rows();
  for(int i = 0; (i < len) && _str; ++i){
    for(int j = 0; (j < _num_channels) && _str; ++j){
      sample_type raw_sample = Audio_Sample<Encoding>::encode(s(i,j));
      _str << raw_sample;
    }
  }
  if(_str){
    _cur_pos += len;
    return true;
  }
  return false;
}

//make audio samples writer for stream output
template<typename Binary_Out_Stream>
audiofile_details::Audio_Samples_Writer_If* make_writer(Binary_Out_Stream& str,
    std::streamoff start_pos, Audio_Stream_Description* d)
{
  Audio_Encoding encoding = d->get_encoding();
  int num_channels = d->get_num_channels();
  switch(encoding){
    case enc_mulaw8:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_mulaw8>(str,start_pos,num_channels);
    case enc_alaw8:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_alaw8>(str,start_pos,num_channels);
    case enc_linear8:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_linear8>(str,start_pos,num_channels);
    case enc_linear16:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_linear16>(str,start_pos,num_channels);
    case enc_linear24:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_linear24>(str,start_pos,num_channels);
    case enc_linear32:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_linear32>(str,start_pos,num_channels);
    case enc_float:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_float>(str,start_pos,num_channels);
    case enc_double:
      return new Audio_Samples_Writer<Binary_Out_Stream, enc_double>(str,start_pos,num_channels);
    case enc_unknown:
    default:
      return 0;
  }
}

//////////////////////////////////////////////////
//
// SND Header helpers
//
//////////////////////////////////////////////////
//snd header consists of 6 uint32_t fixed fields possibly followed
//by the variable length annotation
static const std::size_t snd_fixed_header_size = 24;

//read_header() reads audio stream information from snd file header.
//returns true if successfull.
//d - pointer to stream description to collect info from snd header
//audio_offset - byte offset of first audio sample inside snd file
//number of audio samples stored in snd file, num_samples
template<typename Binary_In_Stream>
bool read_header(Binary_In_Stream& _str, Audio_Stream_Description* d,
  std::streamoff& audio_offset, std::streamoff& num_samples)
{
  //SND header fields
  uint32_t magic, hdr_size, data_size, encoding, sampling_rate, num_channels;
  //encoded sample size
  std::size_t sample_size;
  //annotation
  std::string annotation;

  //read fixed fields of snd file header
  _str.seekg(0, std::ios_base::beg);
  _str >> magic >> hdr_size>> data_size>> encoding >> sampling_rate>> num_channels;
  if(!_str) return false;

  //check magic
  if(magic != snd_magic) return false; //invalid magic

  //check header size (we do not verify divisibility of header size by 8 since we still able to
  //read unaligned data samples)
  if(hdr_size < snd_fixed_header_size) return false; //header is too short
  audio_offset = hdr_size;

  //check encoding
  sample_size = encoded_sample_size((Audio_Encoding)encoding);
  if(!sample_size) return false; //unknown or invalid encoding

  //read annotation
  if(hdr_size > snd_fixed_header_size)
  {//annotation is present
    //get annotation length
    std::streamsize ann_length = (std::streamsize)std::min(hdr_size - snd_fixed_header_size, max_annotation_length);
  for(int i = 0; i < ann_length; ++i){
    char s; _str>>s;
    if(_str && s)
      annotation += s;
    else
      break;
  }
  if(!_str) return false; //failed to read annotation
  }

  //compute number of audio samples based on the file length
  _str.seekg(0, std::ios_base::end);
  if(!_str) return false; //failed to seek to the end
  std::streamoff ns = ((std::streamoff)_str.tellg() - hdr_size)/(num_channels * sample_size);

  //update number of samples just read from header
  if(data_size = 0xffffffff){
    //data size was set to unknown in file header, use number of smaples obtained from file length
    num_samples = ns;
  }
  else{
    num_samples = std::min<long unsigned int>((long unsigned int)ns,(std::streamoff)data_size/(num_channels * sample_size));
  }

  //update start position of audio samples
  audio_offset = hdr_size;

  //update stream description
  d->set_encoding((Audio_Encoding) encoding);
  d->set_num_channels(num_channels);
  d->set_sampling_rate(sampling_rate);
  d->set_description(annotation);
  return true;
}

//write_header() writes audio stream information to snd file header.
//it is assumed that stream description pointed by d is initialized with correct values
template<typename Binary_Out_Stream>
bool write_header(Binary_Out_Stream& _str, const Audio_Stream_Description* const d, std::streamoff& audio_offset)
{
  uint32_t hdr_size, data_size, encoding, sampling_rate, num_channels;
  data_size = 0xffffffff;
  encoding = (uint32_t)d->get_encoding();
  sampling_rate = (uint32_t)d->get_sampling_rate();
  num_channels = (uint32_t)d->get_num_channels();

  //compute header size based on fixed fields length and length of the annotation
  uint32_t ann_length = (uint32_t)std::min(d->get_description().length(), max_annotation_length);
  uint32_t padding_length = (8 - ((ann_length+1) % 8)) % 8; //compute padding length
  hdr_size = (uint32_t)snd_fixed_header_size + ann_length + padding_length + 1;

  //position stream pointer at the beginning of the file
  _str.seekp(0,std::ios_base::beg);
  if(!_str) return false;

  //write fixed-sized part of the header
  _str << snd_magic << hdr_size << data_size << encoding << sampling_rate << num_channels;
  if(!_str) return false;

  //write annotationn and padding
  _str.write(d->get_description().c_str(), ann_length);
  for(uint32_t i = 0; (i < (padding_length + 1)) && _str; ++i) _str << '\0';

  if(!_str) return false;
  audio_offset = hdr_size;

  return true;
}

//update_num_samples_in_header() updates numder of samples in snd file header
template<typename Binary_Out_Stream>
bool update_num_samples_in_header(Binary_Out_Stream& _str, const Audio_Stream_Description* const d, std::streamoff num_samples)
{
  uint32_t data_size = (uint32_t) std::min<long unsigned int>((long unsigned int)(num_samples * 
    d->get_num_channels() * encoded_sample_size(d->get_encoding())), 
    (long unsigned int)0xffffffff);
  _str.seekp(2*sizeof(uint32_t),std::ios_base::beg);
  if(!_str) return false;
  _str << data_size;
  if(!_str) return false;
  return true;
}


//////////////////////////////////////////////////
//
// SND_In_File
//
//////////////////////////////////////////////////

SND_In_File::SND_In_File():_samples_reader(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{
}

SND_In_File::SND_In_File(const char *fname):_samples_reader(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{
  open(fname);
}

bool SND_In_File::open(const char *fname)
{
  //try to reopen the stream
  if (_str.is_open()) close();
  _str.clear();
  _str.open(fname, bfstream_base::b_endian);
  if (!_str) return false;

  //read header and update description
  std::streamoff audio_offset;
  if (!read_header(_str,_description, audio_offset,_num_samples)) {
    _str.close();
    return false;
  }

  //create samples reader
  it_assert(_samples_reader == 0, "SND_In_File::open: samples reader was not deallocated properly.");
  _samples_reader = make_reader(_str,audio_offset,_description);
  return true;
}

void SND_In_File::close()
{
  //close stream
  if(_str.is_open()) _str.close();
  //dispose reader
  if(_samples_reader){
    delete _samples_reader;
    _samples_reader = 0;
  }
  //reset description
  _num_samples = 0;
  *_description = Audio_Stream_Description();
}

SND_In_File::~SND_In_File()
{
  //close file and dispose description
  close();
  delete _description;
}


//////////////////////////////////////////////////
//
// SND_Out_File
//
//////////////////////////////////////////////////
SND_Out_File::SND_Out_File():_samples_writer(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{
}

SND_Out_File::SND_Out_File(const char *fname, const Audio_Stream_Description& d):
  _samples_writer(0), _description(new Audio_Stream_Description), _num_samples(0)
{
  open(fname, d);
}

bool SND_Out_File::open(const char *fname, const Audio_Stream_Description& d)
{
  //check if we have a valid description
  if(!is_valid(d)) return false;
  //try to reopen the stream
  if (_str.is_open()) close();
  _str.clear();
  _str.open(fname, true, bfstream_base::b_endian);
  if (!_str) return false;
  //init description and write header
  *_description = d;
  std::streamoff audio_offset;
  if (!write_header(_str, &d, audio_offset)){
    _str.close();
    return false;
  }
  //create samples writer
  it_assert(_samples_writer == 0, "SND_Out_File::open: samples writer was not deallocated properly.");
  _samples_writer = make_writer(_str,audio_offset,_description);
  _num_samples = 0;
  return true;
}

void SND_Out_File::close()
{
  //update file header with written samples counter and close the stream
  if(_str.is_open()){
    update_num_samples_in_header(_str,_description,_num_samples);
    _str.close();
  }
  //dispose samples writer
  if(_samples_writer){
    delete _samples_writer;
    _samples_writer = 0;
  }
  //reset description
  _num_samples = 0;
  *_description = Audio_Stream_Description();
}

SND_Out_File::~SND_Out_File()
{
  //close file and dispose description
  close();
  delete _description;
}
//////////////////////////////////////////////////
//
// SND_IO_File
//
//////////////////////////////////////////////////

SND_IO_File::SND_IO_File():
  _samples_reader(0), _samples_writer(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{

}

SND_IO_File::SND_IO_File(const char *fname):
  _samples_reader(0), _samples_writer(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{
  open(fname);
}

SND_IO_File::SND_IO_File(const char *fname, const Audio_Stream_Description& d):
  _samples_reader(0), _samples_writer(0),
  _description(new Audio_Stream_Description), _num_samples(0)
{
  open(fname, d);
}

//open IO file. Pick up information about audio stream from file header
bool SND_IO_File::open(const char *fname)
{
  //try to reopen the stream
  if (_str.is_open()) close();
  _str.clear();
  _str.open(fname, false, bfstream_base::b_endian);
  if (!_str) return false;
  //reader header and update description
  std::streamoff audio_offset;
  if (!read_header(_str,_description, audio_offset,_num_samples)) {
    _str.close();
    return false;
  }
  //create reader and writer for audio samples
  it_assert(_samples_reader == 0, "SND_IO_File::open: samples reader was not deallocated properly.");
  _samples_reader = make_reader(_str,audio_offset,_description);
  it_assert(_samples_writer == 0, "SND_IO_File::open: samples writer was not deallocated properly.");
  _samples_writer = make_writer(_str,audio_offset,_description);
  return true;

}

//open file for IO with new descrition - truncates file contents
bool SND_IO_File::open(const char *fname, const Audio_Stream_Description& d)
{
  //check if provided decription is valid
  if(!is_valid(d)) return false;

  //try reopen the stream
  if (_str.is_open()) close();
  _str.clear();
  _str.open(fname, true, bfstream_base::b_endian);
  if (!_str) return false;

  //init description and write file header
  *_description = d;
  std::streamoff audio_offset;
  if (!write_header(_str, &d, audio_offset)){
    _str.close();
    return false;
  }
  //create reader and writer for audio samples
  it_assert(_samples_reader == 0, "SND_IO_File::open: samples reader was not deallocated properly.");
  _samples_reader = make_reader(_str,audio_offset,_description);
  it_assert(_samples_writer == 0, "SND_IO_File::open: samples writer was not deallocated properly.");
  _samples_writer = make_writer(_str,audio_offset,_description);

  _num_samples = 0;
  return true;

}


void SND_IO_File::close()
{
  //close the stream and update number of written samples
  if(_str.is_open()){
    update_num_samples_in_header(_str,_description,_num_samples);
    _str.close();
  }
  //dispose reader and writer of audio samples
  if(_samples_writer){
    delete _samples_writer;
    _samples_writer = 0;
  }

  if(_samples_reader){
    delete _samples_reader;
    _samples_reader = 0;
  }
  //reset description
  _num_samples = 0;
  *_description = Audio_Stream_Description();
}

SND_IO_File::~SND_IO_File()
{
  //close audio file and dispose description
  close();
  delete _description;
}

} // namespace itpp

//! \endcond
