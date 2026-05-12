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
#include <itpp/srccode/g711.h>

namespace itpp
{

namespace g711_details {
  //define tables used in compression  and expansion algorithms
  uint8_t G711_Base_Properties::compression_table[128] =
  {
    0,  1,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
    5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
    6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
    6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7
  };
  int16_t MuLaw_Properties::expansion_table[256] =
  {
   -8031, -7775, -7519, -7263, -7007, -6751, -6495, -6239,
   -5983, -5727, -5471, -5215, -4959, -4703, -4447, -4191,
   -3999, -3871, -3743, -3615, -3487, -3359, -3231, -3103,
   -2975, -2847, -2719, -2591, -2463, -2335, -2207, -2079,
   -1983, -1919, -1855, -1791, -1727, -1663, -1599, -1535,
   -1471, -1407, -1343, -1279, -1215, -1151, -1087, -1023,
    -975,  -943,  -911,  -879,  -847,  -815,  -783,  -751,
    -719,  -687,  -655,  -623,  -591,  -559,  -527,  -495,
    -471,  -455,  -439,  -423,  -407,  -391,  -375,  -359,
    -343,  -327,  -311,  -295,  -279,  -263,  -247,  -231,
    -219,  -211,  -203,  -195,  -187,  -179,  -171,  -163,
    -155,  -147,  -139,  -131,  -123,  -115,  -107,   -99,
     -93,   -89,   -85,   -81,   -77,   -73,   -69,   -65,
     -61,   -57,   -53,   -49,   -45,   -41,   -37,   -33,
     -30,   -28,   -26,   -24,   -22,   -20,   -18,   -16,
     -14,   -12,   -10,    -8,    -6,    -4,    -2,     0,
    8031,  7775,  7519,  7263,  7007,  6751,  6495,  6239,
    5983,  5727,  5471,  5215,  4959,  4703,  4447,  4191,
    3999,  3871,  3743,  3615,  3487,  3359,  3231,  3103,
    2975,  2847,  2719,  2591,  2463,  2335,  2207,  2079,
    1983,  1919,  1855,  1791,  1727,  1663,  1599,  1535,
    1471,  1407,  1343,  1279,  1215,  1151,  1087,  1023,
     975,   943,   911,   879,   847,   815,   783,   751,
     719,   687,   655,   623,   591,   559,   527,   495,
     471,   455,   439,   423,   407,   391,   375,   359,
     343,   327,   311,   295,   279,   263,   247,   231,
     219,   211,   203,   195,   187,   179,   171,   163,
     155,   147,   139,   131,   123,   115,   107,    99,
    93,    89,    85,    81,    77,    73,    69,    65,
    61,    57,    53,    49,    45,    41,    37,    33,
    30,    28,    26,    24,    22,    20,    18,    16,
    14,    12,    10,     8,     6,     4,     2,     0
  };
  int16_t ALaw_Properties::expansion_table[256] =
  {
    -688,  -656,  -752,  -720,  -560,  -528,  -624,  -592,
    -944,  -912, -1008,  -976,  -816,  -784,  -880,  -848,
    -344,  -328,  -376,  -360,  -280,  -264,  -312,  -296,
    -472,  -456,  -504,  -488,  -408,  -392,  -440,  -424,
   -2752, -2624, -3008, -2880, -2240, -2112, -2496, -2368,
   -3776, -3648, -4032, -3904, -3264, -3136, -3520, -3392,
   -1376, -1312, -1504, -1440, -1120, -1056, -1248, -1184,
   -1888, -1824, -2016, -1952, -1632, -1568, -1760, -1696,
     -43,   -41,   -47,   -45,   -35,   -33,   -39,   -37,
     -59,   -57,   -63,   -61,   -51,   -49,   -55,   -53,
     -11,    -9,   -15,   -13,    -3,    -1,    -7,    -5,
     -27,   -25,   -31,   -29,   -19,   -17,   -23,   -21,
    -172,  -164,  -188,  -180,  -140,  -132,  -156,  -148,
    -236,  -228,  -252,  -244,  -204,  -196,  -220,  -212,
     -86,   -82,   -94,   -90,   -70,   -66,   -78,   -74,
    -118,  -114,  -126,  -122,  -102,   -98,  -110,  -106,
     688,   656,   752,   720,   560,   528,   624,   592,
     944,   912,  1008,   976,   816,   784,   880,   848,
     344,   328,   376,   360,   280,   264,   312,   296,
     472,   456,   504,   488,   408,   392,   440,   424,
    2752,  2624,  3008,  2880,  2240,  2112,  2496,  2368,
    3776,  3648,  4032,  3904,  3264,  3136,  3520,  3392,
    1376,  1312,  1504,  1440,  1120,  1056,  1248,  1184,
    1888,  1824,  2016,  1952,  1632,  1568,  1760,  1696,
    43,    41,    47,    45,    35,    33,    39,    37,
    59,    57,    63,    61,    51,    49,    55,    53,
    11,     9,    15,    13,     3,     1,     7,     5,
    27,    25,    31,    29,    19,    17,    23,    21,
     172,   164,   188,   180,   140,   132,   156,   148,
     236,   228,   252,   244,   204,   196,   220,   212,
    86,    82,    94,    90,    70,    66,    78,    74,
     118,   114,   126,   122,   102,    98,   110,   106
  };

  /*
  //These two functions generate and print compression-expansion tables defined above.
  //They can be used to create these tables dynamically.
  void generate_tables()
  {
    //generation of table used in compression - compute segment
    //numbers for all possible values of MSBs
    uint8_t cmp_val = 1;
    uint8_t segno = 0;
    for(uint8_t i = 0; i < 128; ++i) {
      if(i == cmp_val) { //detect leading 1
        cmp_val <<= 1; segno++;
      }
      G711_Base_Properties::compression_table[i] = segno;
    }

    //generation of expansion table (inverse encoding)
    for(uint16_t i = 0; i < 256; ++i){
      uint16_t mantissa = ~i; // 1's complement of the input value
      uint8_t exponent = (mantissa >> 4) & (0x0007); // extract exponent
      uint8_t segment = exponent + 1; //compute segment number
      mantissa = (mantissa & 0x0f) | 0x10; // extract mantissa and add leading 1

      int16_t rounding_value = 1 << exponent; //rounding value is equal to the half of the LSB
      int16_t restored_magnitude = (mantissa << segment) + rounding_value
        - MuLaw_Properties::magnitude_offset; //correct magnitude offset introduced during compression
      MuLaw_Properties::expansion_table[i] = i & 0x80 ? restored_magnitude : -restored_magnitude;
    }

    //generation of expansion table (inverse encoding)
    for(uint16_t i = 0; i < 256; ++i){
      uint16_t mantissa = i ^ 0x55; //remove even bits toggle during compression
      uint8_t exponent = (mantissa >> 4) & (0x0007); //extract 3 bits of exponent from bits 6..4
      mantissa = (mantissa & 0x0f); //get mantissa from bits 3..0
      int16_t rounding_value;
      if(exponent > 0) {
        mantissa |= 0x10; ///add leading 1
        rounding_value = 1 << (exponent - 1) ; //rounding value is equal to the half of the LSB
      }
      else {
        mantissa <<= 1; //scale mantissa a to align with initial dynamic range
        rounding_value = 1; //rounding value is equal to the LSB
      }
      //restore magnitude and sign
      int16_t restored_magnitude = (mantissa << (exponent)) + rounding_value;
      ALaw_Properties::expansion_table[i] = i & 0x80 ? restored_magnitude : -restored_magnitude;
    }
  }

  void print_tables()
  {
    std::cout << "compression:" << std::endl;
    int j = 0;
    for(int i = 0; i < 128; ++i, ++j) {
      if(j == 16){
        std::cout << std::endl; j = 0;
      }
      std::cout <<"  "<< (int)G711_Base_Properties::compression_table[i];
      if(i != 127) std::cout<<',';
    }
    std::cout << std::endl;
    std::cout << "u expnasion:" << std::endl;
    j = 0;
    for(int i = 0; i < 256; ++i, ++j) {
      if(j == 8) {
        std::cout << std::endl; j = 0;
      }
      std::cout <<" "<< std::setw(5) << (int)MuLaw_Properties::expansion_table[i];
      if(i != 255) std::cout<<',';
    }
    std::cout << std::endl;
    std::cout << "a expansion:" << std::endl;
    j = 0;
    for(int i = 0; i < 256; ++i, ++j) {
      if(j == 8) {
        std::cout << std::endl; j = 0;
      }
      std::cout <<" "<< std::setw(5) << (int)ALaw_Properties::expansion_table[i];
      if(i != 255) std::cout<<',';
    }
    std::cout << std::endl;
  }
  */
}

}

