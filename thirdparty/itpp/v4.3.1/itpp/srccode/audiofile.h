/*!
 * \file
 * \brief Definitions of audio Audio classes and functions
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

#ifndef AUDIOFILE_H
#define AUDIOFILE_H

#include <string>
#include <algorithm>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/binfile.h>
#include <itpp/srccode/audiosample.h>
#include <itpp/itexports.h>


namespace itpp
{


/*!
\addtogroup audio
\section audiostreams Audio Streams.

Audio streams are used to handle audio files in Sun/NeXT format (the .au format). SND file
consists of header followed by the audio samples. Number of channels, samples encoding
and sampling rate are stored inside the file header. Header can also contain annotation of
variable length describing the stream contents.
 Library provides Audio_Stream_Description class to describe the audio stream and
 SND_In_File/SND_Out_File/SND_InOut_File classes to support handling of SND files.
*/


/*!
\ingroup audio
\brief Description of audio stream.

This class holds information about the stream of audio samples. Information includes
samples encoding, number of channels and sampling rate. Stream can be annotated
via set_description() method.
*/
class Audio_Stream_Description
{
public:
  //! Default ctor - creates uninitialized description
  Audio_Stream_Description():_encoding(enc_unknown), _sampling_rate(0), _num_channels(0){}
  //! Construct with stream parameters: encoding \a e, sampling rate \a sr and number of audio channels \a nc
  Audio_Stream_Description(Audio_Encoding e, int sr, int nc = 1):
  _encoding(e), _sampling_rate(sr), _num_channels(nc){}
  //! Set encoding of audio samples
  Audio_Stream_Description& set_encoding(Audio_Encoding e) {_encoding = e; return *this;}
  //! Set sampling rate (samples per second)
  Audio_Stream_Description& set_sampling_rate(int sr) {_sampling_rate = sr; return *this;}
  //! Set number of audio channels
  Audio_Stream_Description& set_num_channels(int nc) {_num_channels = nc; return *this;}
  //! Set stream annotation
  Audio_Stream_Description& set_description(const std::string& d) {_description = d; return *this;}
  //! Get encoding of audio samples
  Audio_Encoding get_encoding() const {return _encoding;}
  //! Get sampling rate (samples per second)
  int get_sampling_rate() const {return _sampling_rate;}
  //! Get number of audio channels
  int get_num_channels() const {return _num_channels;}
  //! Get stream annotation
  const std::string& get_description() const {return _description;}
private:
  //! Encoding of audio samples
  Audio_Encoding _encoding;
  //! Sampling rate
  int _sampling_rate;
  //! Number of audio channels
  int _num_channels;
  //! Stream annotation text
  std::string _description;
};

//! validity check for stream description \a d
inline bool is_valid(const Audio_Stream_Description& d)
{
  if(!encoded_sample_size(d.get_encoding())) return false;
  if(!d.get_num_channels()) return false;
  return true;
}

//! \cond
namespace audiofile_details{

  //abstract interfaces to read and write audio samples to streams

  class Audio_Samples_Reader_If
  {
  public:
    virtual bool read_sample(double& s, int ch) = 0;
    virtual vec read_channel(int n, int ch) = 0;
    virtual mat read(int n) = 0;
    virtual std::streamoff tell() const = 0;
    virtual bool seek(std::streamoff n) = 0;
    virtual std::streamoff num_samples() = 0;
    virtual ~Audio_Samples_Reader_If() {}
  };


  class Audio_Samples_Writer_If
  {
  public:
    virtual bool write_sample(const double& s, int ch) = 0;
    virtual bool write_channel(const vec& s, int ch) = 0;
    //Write n samples to audio channel ch
    virtual bool write(const mat& s) = 0;
    virtual std::streamoff tell() const = 0;
    virtual bool seek(std::streamoff n) = 0;
    virtual std::streamoff num_samples() = 0;
    virtual ~Audio_Samples_Writer_If() {}
  };
}
//! \endcond

/*!
\ingroup audio
\brief Class to read audio data from au file

Input stream of audio samples uses binary stream to get encoded audio data from
snd audio file. Audio can be read as single sample from current read position in audio stream,
as a vector of samples containing data from single audio channel or as matrix with audio channels
stored columnwise.

Following example illustratates read operations with SND_In_File
\code
#include <itpp/srccode/audiofile.h>
using namespace itpp;

int main() {
  //create audio stream
  SND_In_File f_in("inptut.au");
  //get description
  Audio_Stream_Description d = f_in.get_description();
  //read 100 audio samples if file contains stereo data on 8 kHz sampling rate
  mat in(100,2); vec first_channel(100); vec second_channel(100);
  if((d.get_num_channels == 2) && (d.get_sampling_rate() == 8000) && (f_in.num_samples()>=100))
  {
     in = f_in.read(100);
     f_in.seek_read(0); //reposition to the first sample
     first_channel = f_in.read_channel(100); //read first channel
     f_in.seek_read(0); //reposition to the first sample
     second_channel = f_in.read_channel(100,1); //read second channel
  }
  return 0;
}
\endcode

*/
class ITPP_EXPORT SND_In_File
{
public:
  //! Default constructor - creates uninitialized stream
  SND_In_File();
  //! Constructor from file name \a fname
  SND_In_File(const char* fname);
  //! Stream destructor
  ~SND_In_File();
  //! Open the file \a fname
  bool open(const char* fname);
  //! Close the file.
  void close();
  //! Get stream description
  Audio_Stream_Description get_description() const {return *_description;}
  //! Go to sample number \a pos
  bool seek_read(std::streamoff pos)
  {
    if((pos > _num_samples) || (pos < 0))
      return false;

    if(_samples_reader)
      return _samples_reader->seek(pos);
    else
      return false;
  }
  //! Get current position in samples.
  std::streamoff tell_read()
  {
    if(_samples_reader)
      return _samples_reader->tell();
    else
      return -1;
  }
  //! Get number of samples in stream
  std::streamoff num_samples() const {return _num_samples;}
  //! Read single sample \a s at current position to channel \a ch.
  bool read_sample(double& s, int ch = 0)
  {
    if(_samples_reader)
      return _samples_reader->read_sample(s,ch);
    else
      return false;
  }
  //! Read \a n samples from channel \a ch starting at current position
  vec read_channel(int n, int ch = 0)
  {
    if(_samples_reader)
      return _samples_reader->read_channel(n,ch);
    else
      return vec();
  }
  //! Read \a n samples from all channels starting at current position into matrix
  mat read(int n)
  {
    if(_samples_reader)
      return _samples_reader->read(n);
    else
      return mat();
  }
private:
  //! Binary stream
  bifstream _str;
  //! Number of samples
  std::streamoff _num_samples;
  //! Samples Reader
  audiofile_details::Audio_Samples_Reader_If* _samples_reader;
  //! Stream Description
  Audio_Stream_Description* _description;
};

/*!
\ingroup audio
\brief A class to write SND-files (the .au format)

Output stream uses underlying binary stream to write audio data to the audio file in au format.
Audio can be written sample-by sample, channel-wise or taken column-wise from the input matrix

Following example illustrates write operations with SND_Out_File:
\code
#include <itpp/srccode/audiofile.h>
using namespace itpp;

int main() {
  //create stream description to store 8kHz 16-bit PCM stereo audio data
  d = Audio_Stream_Description(enc_linear16, 8000,2);
  //create audio stream
  SND_Out_File f_out("inptut.au",d);
  //fill input with 10004 Hz test tone
  mat in(100,2);
  for(int i = 0; i < 100; ++i){
    in(i,0) = sin(2*pi*i*1004/8000);
    in(i,1) = cos(2*pi*i*1004/8000);
  }
  //write audio samples
  f_out.write(in);
  //zero sample 5 in channel 0
  f_out.seek_write(5); f_out.write_sample(0,0);
  //zero samples 5..7 in channel 1
  vec z(3); z(0) = 0.0; z(1) = 0.0; z(2) = 0.0;
  f_out.seek_write(5); f_out.write_channel(z,1);

  return 0;
}
\endcode
*/
class ITPP_EXPORT SND_Out_File
{
public:
  //! Default constructor - creates uninitialized stream
  SND_Out_File();
  //! Constructor from file name \a fname and stream description \a d
  SND_Out_File(const char *fname, const Audio_Stream_Description& d);
  //! Stream destructor
  ~SND_Out_File();
  //! Open the file \a fname with stream description \a d
  bool open(const char *fname, const Audio_Stream_Description& d);
  //! Close the file.
  void close();
  //! Get stream description
  Audio_Stream_Description get_description() const {return *_description;}
  //! Set current position to write to \a pos (samples).
  bool seek_write(std::streamoff pos)
  {
    if((pos > _num_samples) || (pos < 0))
      return false;

    if(_samples_writer)
      return _samples_writer->seek(pos);
    else
      return false;
  }
  //! Get current position in samples.
  std::streamoff tell_write()
  {
    if(_samples_writer)
      return _samples_writer->tell();
    else
      return -1;
  }
  //! Get number of samples in stream
  std::streamoff num_samples() const {return _num_samples;}
  //! Write single sample \a s at current position to channel \a ch
  bool write_sample(const double &s, int ch = 0)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write_sample(s,ch);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
  //! Write the vector \a v to channel \a ch starting at current position
  bool write_channel(const vec &v, int ch = 0)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write_channel(v,ch);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
  //! Write audio channels from columns of the matrix \a m starting at current position
  bool write(const mat &m)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write(m);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
private:
  //! Binary stream
  bofstream _str;
  //! Number of samples
  std::streamoff _num_samples;
  //! Samples Writer
  audiofile_details::Audio_Samples_Writer_If* _samples_writer;
  //! Stream Description
  Audio_Stream_Description* _description;
};

/*!
\ingroup audio
\brief A class for doing both input and output of audio samples.

SND_IO_File provides facilities for doing both input and output of audio samples.
*/
class ITPP_EXPORT SND_IO_File
{
public:
  //! Constructor - creates uninitialized stream
  SND_IO_File();
  //! Open the file \a fname, check file header.
  SND_IO_File(const char *fname);
  //! Open the file \a fname, truncate and overwrite header with description \a d.
  SND_IO_File(const char *fname, const Audio_Stream_Description& d);
  //! Stream destructor
  ~SND_IO_File();
  //! Open the file \a fname, check file header.
  bool open(const char *fname);
  //! Open the file \a fname, truncate and overwrite header with description \a d.
  bool open(const char *fname, const Audio_Stream_Description& d);
  //! Close the file.
  void close();
  //! Get stream description
  Audio_Stream_Description get_description() const {return *_description;}
  //! Set current position to read from \a pos (samples).
  bool seek_read(std::streamoff pos)
  {
    if((pos > _num_samples) || (pos < 0))
      return false;

    if(_samples_reader)
      return _samples_reader->seek(pos);
    else
      return false;
  }
  //! Get current position to read from  in samples.
  std::streamoff tell_read()
  {
    if(_samples_reader)
      return _samples_reader->tell();
    else
      return -1;
  }
  //! Set current position to write to \a pos (samples).
  bool seek_write(std::streamoff pos)
  {
    if((pos > _num_samples) || (pos < 0))
      return false;

    if(_samples_writer)
      return _samples_writer->seek(pos);
    else
      return false;
  }
  //! Get current position to write in samples.
  std::streamoff tell_write()
  {
    if(_samples_writer)
      return _samples_writer->tell();
    else
      return -1;
  }
  //! Get number of samples in stream
  std::streamoff num_samples() const {return _num_samples;}
  //! Read single sample \a s at current position to channel \a ch.
  bool read_sample(double& s, int ch = 0)
  {
    if(_samples_reader)
      return _samples_reader->read_sample(s,ch);
    else
      return false;
  }
  //! Read \a n samples from channel \a ch starting at current position
  vec read_channel(int n, int ch = 0)
  {
    if(_samples_reader)
      return _samples_reader->read_channel(n,ch);
    else
      return vec();
  }
  //! Read \a n samples from all channels starting at current position
  mat read(int n)
  {
    if(_samples_reader)
      return _samples_reader->read(n);
    else
      return mat();
  }

  //! Write single sample \a s at current position to channel \a ch
  bool write_sample(const double &s, int ch = 0)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write_sample(s,ch);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
  //! Write the vector \a v to channel \a ch starting at current position
  bool write_channel(const vec &v, int ch = 0)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write_channel(v,ch);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
  //! Write audio channels from columns of the matrix \a m  starting at current position
  bool write(const mat &m)
  {
    if(_samples_writer){
      bool ret = _samples_writer->write(m);
      if(ret){
        _num_samples = std::max(_num_samples, _samples_writer->tell());
      }
      return ret;
    }
    else
      return false;
  }
private:
  //! Binary stream
  bfstream _str;
  //! Number of samples
  std::streamoff _num_samples;
  //! Samples Reader
  audiofile_details::Audio_Samples_Reader_If* _samples_reader;
  //! Samples Writer
  audiofile_details::Audio_Samples_Writer_If* _samples_writer;
  //! Stream Description
  Audio_Stream_Description* _description;
};

//! \cond
/*
   \brief SAP audio file input class
   \ingroup audio

   ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_In_File : virtual public Audio_File {
  public:
  // Constructor
  SAP_In_File();
  // Open the file {\em fname}.
  SAP_In_File(const char *fname);
  // Destructor
  virtual ~SAP_In_File() { close(); }

  // Open the file {\em fname}.
  virtual bool open(const char *fname);
  // Close the file.
  virtual void close();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual bool seek_read(int pos);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual int tell_read();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool read(vec &v);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool read(vec &v, int n);

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  const char *get_header() { return header; }

  protected:
  char header[SAP_HEADER_SIZE];
  };
*/

/*
  \brief SAP audio file output class
  \ingroup audio

  ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_Out_File : virtual public Audio_File {
  public:
  // Constructor
  SAP_Out_File();
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  SAP_Out_File(const char *fname, const char *hdr);
  // Destructor
  virtual ~SAP_Out_File() { close(); }

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool open(const char *fname, const char *hdr);

  // Old def. Removed since Sun CC gave warning.
  //virtual bool open(const char *fname, const char *hdr);

  // Close the file
  virtual void close();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool seek_write(int pos);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  int tell_write();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual bool write(const vec &v);
  };
*/

/*
  \brief SAP audio file input and output class
  \ingroup audio

  ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_IO_File : public SAP_In_File, public SAP_Out_File {
  public:
  // Constructor
  SAP_IO_File() { }
  // Open the file {\em fname}.
  SAP_IO_File(const char *fname) { open(fname); }
  // Destructor
  virtual ~SAP_IO_File() { close(); }

  // Open the file {\em fname}.
  virtual bool open(const char *fname);
  // Close the file
  virtual void close();
  };
*/
//! \endcond


/*! \addtogroup audio */
//!@{


//! Read audio channel
inline vec snd_read_channel(const char *fname, int ch = 0)
{
  SND_In_File f(fname);
  int ns = (int)std::min(f.num_samples(), (std::streamoff)std::numeric_limits<int>::max());
  return f.read_channel(ns,ch);
}
//! Read \a len audio channel samples starting at position \a beg
inline vec snd_read_channel(const char *fname, int ch, int len, std::streamoff beg = 0)
{
  vec ret; SND_In_File f(fname);
  if(f.seek_read(beg)) ret = f.read_channel(len,ch);
  return ret;
}

//! Read audio data
inline mat snd_read(const char *fname)
{
  SND_In_File f(fname);
  int ns = (int)std::min(f.num_samples(), (std::streamoff)std::numeric_limits<int>::max());
  return f.read(ns);
}
//! Read \a len audio samples starting at position \a beg
inline mat snd_read(const char *fname, int len, std::streamoff beg = 0)
{
  mat ret; SND_In_File f(fname);
  if(f.seek_read(beg)) ret = f.read(len);
  return ret;
}


//! Write audio channel from vector \a s using stream description \a descr
inline bool snd_write_channel(const char *fname, const Audio_Stream_Description& descr, const vec& s, int ch = 0)
{
  SND_Out_File f(fname,descr);
  return f.write_channel(s,ch);
}

//! Write audio data
inline bool snd_write(const char *fname, const Audio_Stream_Description& descr, const mat& s)
{
  SND_Out_File f(fname,descr);
  return f.write(s);
}

//!@}

//! \cond
/*
// Read SAP audio data
bool sap_read(const char *fname, vec &v);
// Read SAP audio data
bool sap_read(const char *fname, vec &v, int beg, int len);
// Write SAP audio data
bool sap_write(const char *fname, const vec &v, const char *hdr);
*/
//! \endcond

} // namespace itpp

#endif // #ifndef AUDIOFILE_H
