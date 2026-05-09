/*!
 * \file
 * \brief Encoding and decoding of audio samples
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

#ifndef AUDIOSAMPLE_H
#define AUDIOSAMPLE_H

#include <climits>
#include <itpp/base/ittypes.h>
#include <itpp/srccode/g711.h>

/*!
\addtogroup audio

\section audiorepresentation Representaion of Audio samples.

ITPP supports various types of representation for audio samples:
- 8/16/32-bit PCM encoding formats
-  IEEE 754 Floating point formats with single and double precision
-  G.711 u-law and A-law encoded samples

*/

namespace itpp
{

/*!
\ingroup audio
\brief Supported encoding types for audio samples.

Numerical values of these constants correspond to the encoding
type identifier in snd file format, introduced by Sun Microsystems.
*/
enum Audio_Encoding {enc_unknown = 0, enc_mulaw8 = 1, enc_alaw8 = 27,
        enc_linear8 = 2,enc_linear16 = 3,enc_linear24 = 4,
        enc_linear32 = 5,enc_float = 6,enc_double = 7};


/*!
\ingroup audio
\brief Helper function for scaling and limiting of audio samples.

This function maps [-1.0, 1.0] doubles to [-max_abs, max_abs] values of type T.
Input values are limited before mapping.
*/

template<typename T, T max_abs> T limit_audio_sample(double s)
{
   //ternary operators allow optimizer to deploy SIMD floating-point instructions
   s < -1.0 ? -1.0 : s > 1.0 ? 1.0 : s;
   return (T)(s*max_abs);
}

/*!
\ingroup audio
\brief Helper function for scaling and limiting of audio samples.

This function maps [-max_abs, max_abs] values of type T to doubles in [-1.0,1.0] interval
*/

template<typename T, T down_scaling> double audio_sample_to_double(T s)
{
  return (1.0/down_scaling) * s;
}

/*!
\ingroup audio
\brief Generic template class for Audio samples.

Specializations of this class provide encoding and decoding facilities for various
representations of audio samples. Encoding inputs are limited to [-1.0,1.0] range.
Decoding outputs are scaled to [-1.0,1.0] range.
*/

template<Audio_Encoding> class Audio_Sample;

/*!
\ingroup audio
\brief uLaw-encoded Audio samples.
*/
template<> class Audio_Sample<enc_mulaw8>
{
public:
  typedef uint8_t enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {
    int16_t l = limit_audio_sample<int16_t, SHRT_MAX>(s);
    return ulaw_compress(l);
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<int16_t, SHRT_MAX>((ulaw_expand(s)));
  }
};

/*!
\ingroup audio
\brief 8-bit PCM encoded audio samples.
*/
template<> class Audio_Sample<enc_linear8>
{
public:
  typedef int8_t enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {
    return limit_audio_sample<enc_sample_type, SCHAR_MAX>(s);
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<enc_sample_type, SCHAR_MAX>(s);
  }
};

/*!
\ingroup audio
\brief 16-bit PCM encoded audio samples
*/
template<> class Audio_Sample<enc_linear16>
{
public:
  typedef int16_t enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {
    return limit_audio_sample<enc_sample_type, SHRT_MAX>(s);
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<enc_sample_type, SHRT_MAX>(s);
  }
};

//! Small class to represent 24-bit PCM samples.
class Sample_24
{
public:
  static const int32_t max_abs_value = (1<<23) - 1;
  explicit Sample_24(uint32_t v = 0):_value(v){}
  uint32_t value() const {return _value;}
  void value(uint32_t v){_value = v;}
private:
  uint32_t _value;
};


//! insertion operator for 24-bit PCM sample
template<typename Binary_Out_Stream>
Binary_Out_Stream& operator<<(Binary_Out_Stream& s, Sample_24 v)
{
  uint32_t sample = v.value();
  char *c = reinterpret_cast<char *>(&sample);
  if(s.get_endianity() == s.get_native_endianity()){
    //stream endian matches machine endian
    s.write(c,3);
  }
  else{
    //stream endian differs from machine endian - reverse order of bytes
    s.put(c[2]); s.put(c[1]); s.put(c[0]);
  }
  return s;
}

//! extraction operator for 24-bit PCM sample
template<typename Binary_In_Stream>
Binary_In_Stream& operator>>(Binary_In_Stream& s, Sample_24& v)
{
  uint32_t sample;
  char *c = reinterpret_cast<char *>(&sample);
  if(s.get_endianity() == s.get_native_endianity()){
    //stream endian matches machine endian
    s.read(c,3);
  }
  else{
    //stream endian differs from machine endian - reverse order of bytes
    s.get(c[2]); s.get(c[1]); s.get(c[0]);
  }
  if(s) v.value(sample);
  return s;
}

/*!
\ingroup audio
\brief 24-bit PCM encoded audio samples.
*/
template<> class Audio_Sample<enc_linear24>
{
public:
  typedef Sample_24 enc_sample_type;
  static const std::size_t enc_sample_size = 3; //3 bytes per sample
  static enc_sample_type encode(const double& s)
  {
    return Sample_24(limit_audio_sample<int32_t, Sample_24::max_abs_value>(s));
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<int32_t, Sample_24::max_abs_value>(s.value());
  }
};


/*!
\ingroup audio
\brief 32-bit PCM encoded audio samples.
*/
template<> class Audio_Sample<enc_linear32>
{
public:
  typedef int32_t enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {
    return limit_audio_sample<enc_sample_type, INT_MAX>(s);
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<enc_sample_type, INT_MAX>(s);
  }
};


/*!
\ingroup audio
\brief Audio samples encoded as floats.

Samples are NOT saturated to +/- 1.0 during conversion to this format.
Encoded values are limited to [INT_MIN,INT_MAX] to avoid overflow on conversion.
*/
template<> class Audio_Sample<enc_float>
{
public:
  typedef float enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {//saturate here to avoid Infinity values
    return (enc_sample_type)(s < -INT_MAX ? -INT_MAX : s > INT_MAX ? INT_MAX : s);
  }
  static double decode(const enc_sample_type& s){return s;}
};

/*!
\ingroup audio
\brief Audio samples encoded as doubles.

Samples are NOT saturated to +/- 1.0 during conversion to this format.
*/
template<> class Audio_Sample<enc_double>
{
public:
  typedef double enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s) {return s;}
  static double decode(const enc_sample_type& s){return s;}
};

/*!
\brief aLaw-encoded Audio samples.
\ingroup audio
*/
template<> class Audio_Sample<enc_alaw8>
{
public:
  typedef uint8_t enc_sample_type;
  static const std::size_t enc_sample_size = sizeof(enc_sample_type);
  static enc_sample_type encode(const double& s)
  {
    int16_t l = limit_audio_sample<int16_t, SHRT_MAX>(s);
    return alaw_compress(l);
  }
  static double decode(const enc_sample_type& s)
  {
    return audio_sample_to_double<int16_t, SHRT_MAX>((alaw_expand(s)));
  }
};

//! Size of encoded sample based on the encoding type \a e.
inline std::size_t encoded_sample_size(Audio_Encoding e)
{
  switch(e) {
  case enc_mulaw8:
    return Audio_Sample<enc_mulaw8>::enc_sample_size;
  case enc_linear8:
    return Audio_Sample<enc_linear8>::enc_sample_size;
  case enc_linear16:
    return Audio_Sample<enc_linear16>::enc_sample_size;
  case enc_linear24:
    return Audio_Sample<enc_linear24>::enc_sample_size;
  case enc_linear32:
    return Audio_Sample<enc_linear32>::enc_sample_size;
  case enc_float:
    return Audio_Sample<enc_float>::enc_sample_size;
  case enc_double:
    return Audio_Sample<enc_double>::enc_sample_size;
  case enc_alaw8:
    return Audio_Sample<enc_alaw8>::enc_sample_size;
  case enc_unknown:
  default:
    return 0;
  }
}

} // namespace itpp

#endif // #ifndef AUDIOFILE_H
