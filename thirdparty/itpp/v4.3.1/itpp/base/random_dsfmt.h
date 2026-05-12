/*!
 * \file
 * \brief C++ implementation of dSFMT random number generator
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

#ifndef RANDOM_DSFMT_H
#define RANDOM_DSFMT_H

#include <itpp/base/ittypes.h>
#include <itpp/base/vec.h>
#include <cstring> // required for memset()
#include <ctime>
#include <limits>
#include <itpp/itexports.h>

#if defined(__SSE2__)
#  include <emmintrin.h>
#endif

namespace itpp
{

namespace random_details
{

/*!
 * \ingroup randgen
 * \brief C++ implementation of dSFMT random number generator.
 *
 * The DSFMT class implements parts of the Double precision SIMD-oriented
 * Fast Mersenne Twister (dSFM) random number generator. DSFMT directly
 * generates double precision floating point random numbers, which have the
 * IEEE Standard for Binary Floating-Point Arithmetic (ANSI/IEEE Std
 * 754-1985) format. DSFMT does not support integer outputs.
 *
 * Visit http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html for
 * more details on the original dSFMT implementation.
 *
 * Here is a copy of the LICENSE.txt file from the dSFMT-src-2.0.tar.gz
 * package, on which C++ DSFMT implementation is based:
 * \verbatim
 * Copyright (c) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of the Hiroshima University nor the names of
 *       its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * \endverbatim
 */
template < int MEXP, int POS1, int SL1, uint64_t MSK1, uint64_t MSK2,
         uint64_t FIX1_V, uint64_t FIX2_V, uint64_t PCV1_V, uint64_t PCV2_V >
class ITPP_EXPORT DSFMT
{

public:
  //make usefull constants available in typedefed definitions
  static const int N = (MEXP - 128) / 104 + 1;
  static const uint64_t FIX1 = FIX1_V;
  static const uint64_t FIX2 = FIX2_V;
  static const uint64_t PCV1 = PCV1_V;
  static const uint64_t PCV2 = PCV2_V;

#if defined(__SSE2__)
  static const uint32_t MSK32_1 =  static_cast<uint32_t>((MSK1 >> 32) & (0xffffffffULL));
  static const uint32_t MSK32_2 =  static_cast<uint32_t>(MSK1 & (0xffffffffULL));
  static const uint32_t MSK32_3 =  static_cast<uint32_t>((MSK2 >> 32) & (0xffffffffULL));
  static const uint32_t MSK32_4 =  static_cast<uint32_t>(MSK2 & (0xffffffffULL));
#endif

  /*!
  \brief DSFMT context structure.

  Shall be a POD type since we want to allocate it in thread-local storage
  gcc and msvc may have problems with non-POD types and threadprivate pragma
  */
  struct Context {
    //! \brief Data structure to hold 128-bit values
    union W128_T {
#if defined(__SSE2__)
      __m128i si;
      __m128d sd;
#endif // __SSE2__
      uint64_t u[2];
      uint32_t u32[4];
      double d[2];
    };
    //! 128-bit data type
    typedef union W128_T w128_t;
    //! 128-bit internal state array
    w128_t status[N + 1];
    //! State array indexing
    int idx;
    //! Last known seed used to initialize context
    unsigned int last_seed;
  };

public:
  //! Constructor using a certain context
  DSFMT(Context& c): _context(c) {}

  /*!
   * \brief Initialise the generator with a new seed.
   *
   * This function initializes the internal state array with a 32-bit
   * integer seed.
   * \param seed a 32-bit integer used as the seed.
   */
  void init_gen_rand(unsigned int seed) {
    uint32_t *psfmt = &_context.status[0].u32[0];
    psfmt[idxof(0)] = seed;
    for(int i = 1; i < (N + 1) * 4; i++) {
      psfmt[idxof(i)] = 1812433253UL
                        * (psfmt[idxof(i - 1)] ^ (psfmt[idxof(i - 1)] >> 30)) + i;
    }
    initial_mask();
    period_certification();
    _context.idx = Nx2;
    _context.last_seed = seed;
  }

  //! Generate uniform [0, UINT_MAX) integer pseudorandom number.
  uint32_t genrand_uint32() {
    uint64_t *psfmt64 = &_context.status[0].u[0];
    if(_context.idx >= Nx2) {
      dsfmt_gen_rand_all();
      _context.idx = 0;
    }
    return (uint32_t)(psfmt64[_context.idx++] & 0xffffffffU);
  }

  /*!
   * \brief Generate uniform [1, 2) double pseudorandom number.
   *
   * This function generates and returns double precision pseudorandom
   * number which distributes uniformly in the range [1, 2).  This is
   * the primitive and faster than generating numbers in other ranges.
   * \c init_gen_rand() must be called before this function.
   * \return double precision floating point pseudorandom number
   */
  double genrand_close1_open2() {
    double *psfmt64 = &_context.status[0].d[0];
    if(_context.idx >= Nx2) {
      dsfmt_gen_rand_all();
      _context.idx = 0;
    }
    return psfmt64[_context.idx++];
  }

  /*!
   * \brief Generate uniform (0, 1) double pseudorandom number.
   *
   * This function generates and returns double precision pseudorandom
   * number which distributes uniformly in the range (0, 1).
   * \c init_gen_rand() must be called before this function.
   * \return double precision floating point pseudorandom number
   */
  double genrand_open_open() {
    double *dsfmt64 = &_context.status[0].d[0];
    union {
      double d;
      uint64_t u;
    } r;

    if(_context.idx >= Nx2) {
      dsfmt_gen_rand_all();
      _context.idx = 0;
    }
    r.d = dsfmt64[_context.idx++];
    r.u |= 1;
    return r.d - 1.0;
  }

private:
  static const int Nx2 = N * 2;
  static const unsigned int SR = 12U;
  //! Endianness flag
  static const bool bigendian;

#if defined(__SSE2__)
  //! Mask data for sse2
  static const __m128i sse2_param_mask;
#endif // __SSE2__

  //! Computations context
  Context& _context;


  /*!
   * This function simulate a 32-bit array index overlapped to 64-bit
   * array of LITTLE ENDIAN in BIG ENDIAN machine.
   */
  static int idxof(int i) { return (bigendian ? (i ^ 1) : i); }

  /*!
   * This function initializes the internal state array to fit the IEEE
   * 754 format.
   */
  void initial_mask() {

    const uint64_t LOW_MASK = 0x000fffffffffffffULL;
    const uint64_t HIGH_CONST = 0x3ff0000000000000ULL;

    uint64_t *psfmt = &_context.status[0].u[0];
    for(int i = 0; i < Nx2; i++) {
      psfmt[i] = (psfmt[i] & LOW_MASK) | HIGH_CONST;
    }
  }

  //! This function certificate the period of 2^{MEXP}-1.
  void period_certification() {
    uint64_t pcv[2] = {PCV1, PCV2};
    uint64_t tmp[2];
    uint64_t inner;

    tmp[0] = (_context.status[N].u[0] ^ FIX1);
    tmp[1] = (_context.status[N].u[1] ^ FIX2);

    inner = tmp[0] & pcv[0];
    inner ^= tmp[1] & pcv[1];
    for(int i = 32; i > 0; i >>= 1) {
      inner ^= inner >> i;
    }
    inner &= 1;
    /* check OK */
    if(inner == 1) {
      return;
    }
    /* check NG, and modification */
#if (PCV2 & 1) == 1
    _context.status[N].u[1] ^= 1;
#else
    uint64_t work;
    for(int i = 1; i >= 0; i--) {
      work = 1;
      for(int j = 0; j < 64; j++) {
        if((work & pcv[i]) != 0) {
          _context.status[N].u[i] ^= work;
          return;
        }
        work = work << 1;
      }
    }
#endif // (PCV2 & 1) == 1
    return;
  }

  /*!
   * This function represents the recursion formula.
   * \param r output 128-bit
   * \param a a 128-bit part of the internal state array
   * \param b a 128-bit part of the internal state array
   * \param lung a 128-bit part of the internal state array (I/O)
   */
  static void do_recursion(typename Context::w128_t *r, typename Context::w128_t *a, typename Context::w128_t *b, typename Context::w128_t *lung) {
#if defined(__SSE2__)
    const unsigned int SSE2_SHUFF = 0x1bU;

    __m128i x = a->si;
    __m128i z = _mm_slli_epi64(x, SL1);
    __m128i y = _mm_shuffle_epi32(lung->si, SSE2_SHUFF);
    z = _mm_xor_si128(z, b->si);
    y = _mm_xor_si128(y, z);

    __m128i v = _mm_srli_epi64(y, SR);
    __m128i w = _mm_and_si128(y, sse2_param_mask);
    v = _mm_xor_si128(v, x);
    v = _mm_xor_si128(v, w);
    r->si = v;
    lung->si = y;
#else // standard C++
    uint64_t t0 = a->u[0];
    uint64_t t1 = a->u[1];
    uint64_t L0 = lung->u[0];
    uint64_t L1 = lung->u[1];
    lung->u[0] = (t0 << SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
    lung->u[1] = (t1 << SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
    r->u[0] = (lung->u[0] >> SR) ^ (lung->u[0] & MSK1) ^ t0;
    r->u[1] = (lung->u[1] >> SR) ^ (lung->u[1] & MSK2) ^ t1;
#endif // __SSE2__
  }

  /*!
   * This function fills the internal state array with double precision
   * floating point pseudorandom numbers of the IEEE 754 format.
   */
  void dsfmt_gen_rand_all() {
    int i;
    typename Context::w128_t *status = _context.status;
    typename Context::w128_t lung = status[N];
    do_recursion(&status[0], &status[0], &status[POS1], &lung);
    for(i = 1; i < N - POS1; i++) {
      do_recursion(&status[i], &status[i], &status[i + POS1], &lung);
    }
    for(; i < N; i++) {
      do_recursion(&status[i], &status[i], &status[i + POS1 - N], &lung);
    }
    status[N] = lung;
  }

};


// ----------------------------------------------------------------------
// typedefs of different RNG
// ----------------------------------------------------------------------

typedef DSFMT < 521, 3, 25,
        0x000fbfefff77efffULL, 0x000ffeebfbdfbfdfULL,
        0xcfb393d661638469ULL, 0xc166867883ae2adbULL,
        0xccaa588000000000ULL, 0x0000000000000001ULL > DSFMT_521_RNG;

typedef DSFMT < 1279, 9, 19,
        0x000efff7ffddffeeULL, 0x000fbffffff77fffULL,
        0xb66627623d1a31beULL, 0x04b6c51147b6109bULL,
        0x7049f2da382a6aebULL, 0xde4ca84a40000001ULL > DSFMT_1279_RNG;

typedef DSFMT < 2203, 7, 19,
        0x000fdffff5edbfffULL, 0x000f77fffffffbfeULL,
        0xb14e907a39338485ULL, 0xf98f0735c637ef90ULL,
        0x8000000000000000ULL, 0x0000000000000001ULL > DSFMT_2203_RNG;

typedef DSFMT < 4253, 19, 19,
        0x0007b7fffef5feffULL, 0x000ffdffeffefbfcULL,
        0x80901b5fd7a11c65ULL, 0x5a63ff0e7cb0ba74ULL,
        0x1ad277be12000000ULL, 0x0000000000000001ULL > DSFMT_4253_RNG;

typedef DSFMT < 11213, 37, 19,
        0x000ffffffdf7fffdULL, 0x000dfffffff6bfffULL,
        0xd0ef7b7c75b06793ULL, 0x9c50ff4caae0a641ULL,
        0x8234c51207c80000ULL, 0x0000000000000001ULL > DSFMT_11213_RNG;

typedef DSFMT < 19937, 117, 19,
        0x000ffafffffffb3fULL, 0x000ffdfffc90fffdULL,
        0x90014964b32f4329ULL, 0x3b8d12ac548a7c7aULL,
        0x3d84e1ac0dc82880ULL, 0x0000000000000001ULL > DSFMT_19937_RNG;



/*!
 * \ingroup randgen
 * \brief Active Generator for random (stochastic) sources.
 *
 * ActiveDSFMT is a typedef of DSFMT class specialization using 19937
 * generation period. Library shall be recompiled if switched to other
 * available algorithm.
 *
 * \sa DSFMT
 */
typedef DSFMT_19937_RNG ActiveDSFMT;



/*! \addtogroup randgen
Some functions to deal with thread-local RNG generation context:
\code
ActiveDSFMT::Context& lc_get();
bool lc_is_initialized();
void lc_mark_initialized();
\endcode
  @{
*/

//! Function to access thread-local context for random numbers generation
ITPP_EXPORT ActiveDSFMT::Context& lc_get();
//! Function to check if thread-local context is initialized
ITPP_EXPORT bool lc_is_initialized();
//! Function to mark thread-local context as initialized
ITPP_EXPORT void lc_mark_initialized();

//!@}

}

} // namespace itpp

#endif // #ifndef RANDOM_DSFMT_H
