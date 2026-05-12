/*!
 * \file
 * \brief Implementation of G.711 logarithmic codecs
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
#ifndef G711_H
#define G711_H

#include <itpp/itexports.h>
#include <utility>
#include <itpp/base/ittypes.h>

/*!
 * \addtogroup audio
 *
 * \section g711 G.711 audio codecs
 * ITU-T G.711 defines two logarithmic codecs used for PCM audio data. Both codecs are widely used
 * for 64 kbps serial transmission in telephone networks.
 *
 * Codec algorithms are based on the perceptual properties of the human ear - weak signals are masked
 * by strong ones, so only the MSBs of the original sample can be kept without significant loss of quality.
 *
 * u-law compression encodes 14-bit input sample as floating point value with 1 bit of sign,
 * 3 bits of exponent and 4 bits of mantissa, representing the MSBs of the magnitude of the
 * original sample.
 *
 * a-law compression generates output in the same format but operates on 13-bit samples, providing
 * somewhat higher output dynamic range at the expense of larger quantization error [SPRA634].
 *
 * Codecs differ in particular details of mantissa and exponent encoding (see [G.711] for details).
 * G.191 specification provides reference implementation and test vectors for G.711 codecs.
 *
 * \section g711_refs G.711 references
 *
 * - [G711] G.711 : Pulse code modulation (PCM) of voice frequencies; ITU-T Recommendation (11/1988)
 * - [G191] G.191 : Software tools for speech and audio coding; ITU-T Recommendation (03/2010)
 * - [SPRA634] Mark A. Castellano, Todd Hiers, Rebecca Ma TMS320C6000 u-Law and A-Law Companding with
 *     Software or the McBSP. Application report, Texas Instruments, April 2000.
 */


namespace itpp {

//! \cond

//forward declarations to make friends visible
std::pair<int16_t,int16_t> ulaw_range();
uint8_t ulaw_compress(int16_t s);
int16_t ulaw_expand(uint8_t s);
std::pair<int16_t,int16_t> alaw_range();
uint8_t alaw_compress(int16_t s);
int16_t alaw_expand(uint8_t s);


namespace g711_details {
  //This makes compression-expansion tables inaccessible in user's code
  //while compression-expansion functions are still defined inline. Also, it
  //hides property classes from itpp namespace.

  //Base properties of G.711 codecs
  class G711_Base_Properties
  {
  protected:
    //! look-up table used in compression functions
    static ITPP_EXPORT uint8_t compression_table[128];
  };

  //u-law algorithm properties
  class MuLaw_Properties : public G711_Base_Properties
  {
    //u-law codec input bitwidth
    static const int input_bitwidth = 14;
    //offset applied to magnitude of the input sample
    static const int16_t magnitude_offset = 33;
    //maximum input value according to G.711
    static const int16_t input_max = 8158;
    //minimum input value according to G.711
    static const int16_t input_min = -8159;
    //table used in u-law expansion
    static ITPP_EXPORT int16_t expansion_table[256];

    friend std::pair<int16_t,int16_t> itpp::ulaw_range();
    friend uint8_t itpp::ulaw_compress(int16_t s);
    friend int16_t itpp::ulaw_expand(uint8_t s);
  };

  //a-law algorithm properties
  class ALaw_Properties : public G711_Base_Properties
  {
    //a-law codec input bitwidth
    static const int input_bitwidth = 13;
    //maximum input value according to G.711
    static const int16_t input_max = 4095;
    //minimum input value according to G.711
    static const int16_t input_min = -4096;
    //table used in u-law expansion
    static ITPP_EXPORT int16_t expansion_table[256];

    friend std::pair<int16_t,int16_t> itpp::alaw_range();
    friend uint8_t itpp::alaw_compress(int16_t s);
    friend int16_t itpp::alaw_expand(uint8_t s);
  };
}
//! \endcond

/*!
  \brief G.711 u-Law compressor input range. Returns (min,max) input values in std::pair.
  \ingroup audio
*/
inline std::pair<int16_t,int16_t> ulaw_range()
{
  using namespace g711_details;
  return std::make_pair(MuLaw_Properties::input_min,MuLaw_Properties::input_max);
}


/*!
  \brief G.711 u-Law compression function. Returns encoded value for sample \a s.
  \ingroup audio
*/
inline uint8_t ulaw_compress(int16_t s)
{
  using namespace g711_details;
  //Limiting  and shifting. Negative samples are 1's complemented to align dynamic
  //ranges of positive and negative numbers. Compressed negative
  //and positive values of equal magnitude are distinguished by the sign bit.
  //As per G.711 spec, resulting magnitude is shifted before compression.
  uint8_t sign; uint16_t shifted_magnitude;
  if(s >= 0){
    if(s > MuLaw_Properties::input_max) s = MuLaw_Properties::input_max;
    shifted_magnitude = s + MuLaw_Properties::magnitude_offset;
    sign = 0xff;
  }
  else{
    if(s < MuLaw_Properties::input_min) s = MuLaw_Properties::input_min;
    shifted_magnitude = (MuLaw_Properties::magnitude_offset - 1) - s;
    sign = 0x7f;
  }
  //use compression table to get the segment number. Segment number corresponds
  //to the exponent value stored in compressed sample.
  uint8_t seg_no = MuLaw_Properties::compression_table[shifted_magnitude>>6];
  //extract 4 MSBs of magnitude, except leading 1, compose it with segment number
  //and sign, store 1's complement value as compressed sample
  uint8_t ret = (seg_no << 4) | ((shifted_magnitude >> (seg_no + 1)) & 0x0f);
  ret ^= sign; //1's complement and add sign
  return ret;
}

/*!
  \brief G.711 u-Law expansion function. Returns decoded value for previously compressed sample \a s.
  \ingroup audio
  Expansion is performed by the table look up.
*/
inline int16_t ulaw_expand(uint8_t s){return g711_details::MuLaw_Properties::expansion_table[s];}

/*!
  \brief G.711 a-Law compressor input range. Returns (min,max) input values in std::pair.
  \ingroup audio
*/
inline std::pair<int16_t,int16_t> alaw_range()
{
  using namespace g711_details;
  return std::make_pair(ALaw_Properties::input_min,ALaw_Properties::input_max);
}

/*!
  \brief G.711 a-Law compression function. Returns encoded value for sample \a s.
  \ingroup audio
*/
inline uint8_t alaw_compress(int16_t s)
{
  using namespace g711_details;
  //Limiting. Negative samples are 1's complemented to align dynamic
  //ranges of positive and negative numbers. Compressed negative
  //and positive values of equal magnitude are distinguished by the sign bit.
  uint8_t sign; uint16_t magnitude;
  if(s >= 0){
    if(s > ALaw_Properties::input_max) s = ALaw_Properties::input_max;
    magnitude = s; sign = 0xd5;
  }
  else{//1's complement to get magnitude
    if(s < ALaw_Properties::input_min) s = ALaw_Properties::input_min;
    magnitude = -1 - s; sign = 0x55;
  }
  //use compression table to get the exponent.
  uint8_t exp_val = ALaw_Properties::compression_table[magnitude>>5];
  uint8_t ret;
  if(exp_val > 0) {
    //extract 4 MSBs of magnitude, except leading 1 and compose it with exponent
    ret = (exp_val << 4) | ((magnitude >> (exp_val)) & 0x0f);
  }
  else{//exp is 0, store 4 MSBs of magnitude
    ret = (uint8_t)magnitude >> 1;
  }
  ret ^= sign; //toggle even bits and add sign
  return ret;
}

/*!
  \brief G.711 u-Law expansion function. Returns decoded value for previously compressed sample \a s.
  \ingroup audio
  Expansion is performed by the table look up.
*/
inline int16_t alaw_expand(uint8_t s){return g711_details::ALaw_Properties::expansion_table[s];}

}

#endif
