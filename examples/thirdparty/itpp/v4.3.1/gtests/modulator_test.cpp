/*!
 * \file
 * \brief 1D and 2D modulators test program
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;

static
void assert_vec_p(const vec &ref, const vec &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol) << line;
  }
}
#define assert_vec(ref, act) assert_vec_p(ref, act, __LINE__)
static
void assert_cvec_p(const cvec &ref, const cvec &act, int line)
{
  static const double tol = 1e-4;
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n).real(), act(n).real(), tol) << line;
    ASSERT_NEAR(ref(n).imag(), act(n).imag(), tol) << line;
  }
}
#define assert_cvec(ref, act) assert_cvec_p(ref, act, __LINE__)

TEST(Modulator, All)
{
  RNG_reset(12345);

  // Test of Modulators

  const int no_symbols = 5;
  const double N0 = 0.1;

  {
    // Modulator_1D (configured as BPSK)
    Modulator_1D  mod("1.0 -1.0", "0 1");
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    vec noise = sqrt(N0) * randn(no_symbols);

    vec tx_symbols = mod.modulate_bits(tx_bits);
    vec rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);

    // modulating bits
    vec ref = "-1 -1 1 -1 -1";
    assert_vec(ref, tx_symbols);
    ref = "-0.903346 -0.634069 0.890343 -1.15583 -0.830113";
    assert_vec(ref, rx_symbols);
    ASSERT_TRUE(tx_bits == decbits);

    tx_symbols = mod.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    // modulating symbol numbers
    ref = "-1 -1 1 -1 -1";
    assert_vec(ref, tx_symbols);
    ref = "-0.903346 -0.634069 0.890343 -1.15583 -0.830113";
    assert_vec(ref, rx_symbols);
    ASSERT_TRUE(tx_sym_numbers == dec_sym_numbers);

    // BPSK (real signal)
    BPSK bpsk;

    bpsk.modulate_bits(tx_bits, tx_symbols);
    rx_symbols = tx_symbols + noise;
    bpsk.demodulate_bits(rx_symbols, decbits);

    // modulating bits
    ref = "-1 -1 1 -1 -1";
    assert_vec(ref, tx_symbols);
    ref = "-0.903346 -0.634069 0.890343 -1.15583 -0.830113";
    assert_vec(ref, rx_symbols);
    ASSERT_TRUE(tx_bits == decbits);

    // BPSK (complex signal)
    BPSK_c bpsk_c;

    cvec tx_csymbols = bpsk_c.modulate_bits(tx_bits);
    cvec rx_csymbols = tx_csymbols + to_cvec(noise, -noise);
    decbits = bpsk_c.demodulate_bits(rx_csymbols);
    vec softbits_approx = bpsk_c.demodulate_soft_bits(rx_csymbols, N0, APPROX);
    vec softbits = bpsk_c.demodulate_soft_bits(rx_csymbols, N0, LOGMAP);

    // modulating bits
    cvec ref_c = "-1+0i -1+0i 1+0i -1+0i -1+0i";
    assert_cvec(ref_c, tx_csymbols);
    ref_c = "-0.903346-0.0966544i -0.634069-0.365931i 0.890343+0.109657i -1.15583+0.155832i -0.830113-0.169887i";
    assert_cvec(ref_c, rx_csymbols);
    ASSERT_TRUE(tx_bits == decbits);
    ref = "-36.1338 -25.3628 35.6137 -46.2333 -33.2045";    
    assert_vec(ref, softbits);
    assert_vec(ref, softbits_approx);
  }

  {
    // Modulator_1D (configured as 4-PAM)
    Modulator_1D mod("-3.0 -1.0 1.0 3.0", "0 1 3 2");
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    vec noise = sqrt(N0) * randn(no_symbols);

    vec tx_symbols = mod.modulate_bits(tx_bits);
    vec rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);

    // modulating bits
    vec ref = "3 -1 3 1 -1";
    assert_vec(ref, tx_symbols);
    ref = "2.81898 -1.09523 2.8743 1.56882 -0.880341";
    assert_vec(ref, rx_symbols);
    ASSERT_TRUE(tx_bits == decbits);

    tx_symbols = mod.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    // modulating symbol numbers
    ref = "1 -3 3 3 1";
    assert_vec(ref, tx_symbols);
    ref = "0.818976 -3.09523 2.8743 3.56882 1.11966";
    assert_vec(ref, rx_symbols);
    ASSERT_TRUE(tx_sym_numbers == dec_sym_numbers);

    // 4-PAM (real signal)
    PAM pam(4);

    pam.modulate_bits(tx_bits, tx_symbols);
    rx_symbols = tx_symbols + noise;
    pam.demodulate_bits(rx_symbols, decbits);

    ref = "-1.34164 0.447214 -1.34164 -0.447214 0.447214";
    assert_vec(ref, tx_symbols);
    ref = "-1.52266 0.351985 -1.46734 0.121605 0.566872";
    assert_vec(ref, rx_symbols);
    bvec ref_b = "1 0 0 1 1 0 0 1 0 1";
    ASSERT_TRUE(ref_b == decbits);

    // 4-PAM (complex signal)
    PAM_c pam_c(4);

    cvec tx_csymbols = pam_c.modulate_bits(tx_bits);
    cvec rx_csymbols = tx_csymbols + to_cvec(noise, -noise);
    decbits = pam_c.demodulate_bits(rx_csymbols);
    vec softbits_approx = pam_c.demodulate_soft_bits(rx_csymbols, N0, APPROX);
    vec softbits = pam_c.demodulate_soft_bits(rx_csymbols, N0, LOGMAP);

    cvec ref_c = "-1.34164+0i 0.447214+0i -1.34164+0i -0.447214+0i 0.447214+0i";
    assert_cvec(ref_c, tx_csymbols);
    ref_c = "-1.52266+0.181024i 0.351985+0.0952288i -1.46734+0.125701i 0.121605-0.568819i 0.566872-0.119659i";
    assert_cvec(ref_c, rx_csymbols);
    ref_b = "1 0 0 1 1 0 0 1 0 1";
    ASSERT_TRUE(ref_b == decbits);
    ref = "-38.4765 11.2383 6.29656 -9.70535 -36.4973 10.2486 2.17534 -13.9308 10.1434 -5.85952";
    assert_vec(ref, softbits);
    ref = "-38.4765 11.2383 6.2965 -9.7035 -36.4972 10.2486 2.17534 -13.8247 10.1405 -5.85948";
    assert_vec(ref, softbits_approx);
  }

  {
    // Modulator_2D (configured as 256-QAM)
    QAM qam(256);
    Modulator_2D mod(qam.get_symbols(), qam.get_bits2symbols());
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = mod.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    cvec ref_c = "-0.536875+0.690268i -0.536875-0.0766965i -0.843661+0.997054i 0.843661-0.536875i -0.690268+0.0766965i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "-0.254627+0.564075i -0.440218+0.144641i -1.01567+1.09629i 0.848682-0.37757i -0.759703-0.322026i";
    assert_cvec(ref_c, rx_symbols);
    ivec ref_i = "73 122 14 162 172";
    ASSERT_TRUE(ref_i == dec_sym_numbers);

    tx_symbols = mod.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);
    vec softbits_approx = mod.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = mod.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    ref_c = "-0.843661+0.230089i -0.997054+1.15045i 0.997054+0.0766965i -0.383482-0.0766965i -0.230089+0.536875i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "-0.561413+0.103896i -0.900397+1.37179i 0.825046+0.175927i -0.378462+0.0826089i -0.299524+0.138153i";
    assert_cvec(ref_c, rx_symbols);
    bvec ref_b = "0 1 0 0 1 1 1 0 0 0 0 0 1 0 1 1 0 1 0 1 0 0 1 1 0 1 0 0 1 1 1 1 0 1 0 0 1 1 0 1";
    ASSERT_TRUE(ref_b == decbits);
    vec ref = "0.764361 -4.53879 1.33914 0.0998816 -5.24815 -0.383786 -1.56908 0.175588 20.8995 6.83599 2.47707 0.707791 "
    "-10.7329 2.17238 -0.289121 -0.303773 1.31869 -3.74479 0.901794 -0.0470311 9.34327 1.56432 -0.762366 -0.165231 "
    "0.605491 -4.76275 1.44412 0.13652 -3.12443 -1.80187 -0.499773 -0.153056 1.02428 -4.16087 1.14289 0.0319552 "
    "-2.36526 -2.49883 0.0507592 -0.205765";
    assert_vec(ref, softbits);
    ref = "0.318739 -3.43093 0.774286 0.151849 -4.06582 -0.160015 -1.09173 0.310573 20.4911 6.48084 2.29924 0.914328 "
    "-9.51492 1.28929 -0.061239 -0.409349 0.608852 -2.61613 0.401456 -0.0691321 8.12793 0.826956 -0.292404 "
    "-0.178184 0.253433 -3.69215 0.9049 0.217156 -2.07144 -0.971982 -0.219892 -0.250697 0.423836 -3.01054 0.564094 "
    "0.0467527 -1.36721 -1.4786 0.0222797 -0.448308";
    assert_vec(ref, softbits_approx);

    // 256-QAM

    tx_symbols = qam.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    dec_sym_numbers = qam.demodulate(rx_symbols);

    ref_c = "-0.536875+0.690268i -0.536875-0.0766965i -0.843661+0.997054i 0.843661-0.536875i -0.690268+0.0766965i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "-0.254627+0.564075i -0.440218+0.144641i -1.01567+1.09629i 0.848682-0.37757i -0.759703-0.322026i";
    assert_cvec(ref_c, rx_symbols);
    ref_i = "73 122 14 162 172";
    ASSERT_TRUE(ref_i == dec_sym_numbers);

    tx_symbols = qam.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    decbits = qam.demodulate_bits(rx_symbols);
    softbits_approx = qam.demodulate_soft_bits(rx_symbols, N0, APPROX);
    softbits = qam.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    ref_c = "-0.843661+0.230089i -0.997054+1.15045i 0.997054+0.0766965i -0.383482-0.0766965i -0.230089+0.536875i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "-0.561413+0.103896i -0.900397+1.37179i 0.825046+0.175927i -0.378462+0.0826089i -0.299524+0.138153i";
    assert_cvec(ref_c, rx_symbols);
    ref_b = "0 1 0 0 1 1 1 0 0 0 0 0 1 0 1 1 0 1 0 1 0 0 1 1 0 1 0 0 1 1 1 1 0 1 0 0 1 1 0 1";
    ASSERT_TRUE(ref_b == decbits);
    ref = "0.764361 -4.53879 1.33914 0.0998816 -5.24815 -0.383786 -1.56908 0.175588 20.8995 6.83599 2.47707 0.707791 "
    "-10.7329 2.17238 -0.289121 -0.303773 1.31869 -3.74479 0.901794 -0.0470311 9.34327 1.56432 -0.762366 -0.165231 "
    "0.605491 -4.76275 1.44412 0.13652 -3.12443 -1.80187 -0.499773 -0.153056 1.02428 -4.16087 1.14289 0.0319552 "
    "-2.36526 -2.49883 0.0507592 -0.205765";
    assert_vec(ref, softbits);
    ref = "0.318739 -3.43093 0.774286 0.151849 -4.06582 -0.160015 -1.09173 0.310573 20.4911 6.48084 2.29924 0.914328 "
    "-9.51492 1.28929 -0.061239 -0.409349 0.608852 -2.61613 0.401456 -0.0691321 8.12793 0.826956 -0.292404 -0.178184 "
    "0.253433 -3.69215 0.9049 0.217156 -2.07144 -0.971982 -0.219892 -0.250697 0.423836 -3.01054 0.564094 0.0467527 "
    "-1.36721 -1.4786 0.0222797 -0.448308";
    assert_vec(ref, softbits_approx);
  }

  {
    // 8-PSK
    PSK psk(8);
    int bps = round_i(psk.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = psk.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = psk.demodulate(rx_symbols);

    cvec ref_c = "0.707107-0.707107i 0.707107-0.707107i 0.707107-0.707107i 1+0i 0.707107+0.707107i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "0.8977-0.790151i 0.613921-1.00098i 0.597143-0.839004i 1.14134+0.25144i 0.978761+0.560701i";
    assert_cvec(ref_c, rx_symbols);
    ASSERT_TRUE(tx_sym_numbers == dec_sym_numbers);

    tx_symbols = psk.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = psk.demodulate_bits(rx_symbols);
    vec softbits_approx = psk.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = psk.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    ref_c = "0-1i 0.707107+0.707107i 0.707107-0.707107i 0-1i 0+1i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "0.190593-1.08304i 0.613921+0.413229i 0.597143-0.839004i 0.141345-0.74856i 0.271654+0.853594i";
    assert_cvec(ref_c, rx_symbols);
    ASSERT_TRUE(tx_bits == decbits);
    vec ref = "-17.8748 9.0654 -3.64905 11.79 6.36194 -2.2495 -8.39621 16.9189 3.53033 -12.2338 6.47187 -2.38768 22.7189 "
    "-1.15859 -9.05583";
    assert_vec(ref, softbits);
    ref = "-17.849 9.03972 -3.64894 11.6879 6.26152 -2.24768 -8.36733 16.8897 3.5301 -12.1443 6.38388 -2.38605 22.505 "
    "-1.15847 -8.842";
    assert_vec(ref, softbits_approx);
  }

  {
    // 16-QAM
    QAM qam(16);
    int bps = round_i(qam.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = qam.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = qam.demodulate(rx_symbols);

    cvec ref_c = "0.948683-0.316228i 0.316228-0.316228i -0.948683-0.948683i -0.316228-0.316228i -0.316228+0.948683i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "0.761834-0.196896i 0.818726-0.145779i -1.08588-1.63427i -0.239432-0.0518489i -0.318584+0.783665i";
    assert_cvec(ref_c, rx_symbols);
    ivec ref_i = "8 8 15 10 2";
    ASSERT_TRUE(ref_i == dec_sym_numbers);

    tx_symbols = qam.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = qam.demodulate_bits(rx_symbols);
    vec softbits_approx = qam.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = qam.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    ref_c = "0.316228-0.316228i 0.948683-0.948683i -0.316228-0.948683i 0.316228-0.948683i -0.948683+0.316228i";
    assert_cvec(ref_c, tx_symbols);
    ref_c = "0.129379-0.196896i 1.45118-0.778235i -0.453421-1.63427i 0.393024-0.684304i -0.95104+0.151209i";
    assert_cvec(ref_c, rx_symbols);
    ASSERT_TRUE(tx_bits == decbits);
    vec ref = "-2.49457 -5.58849 1.63818 -6.53399 -11.8348 1.84392 28.7123 10.3562 -33.3441 12.6721 -5.83419 "
    "-2.26785 -9.72974 0.655668 5.01865 -3.03551 1.91488 -6.22187 -16.0772 4.0298";
    assert_vec(ref, softbits);
    ref = "-2.49055 -5.50945 1.63652 -6.36348 -11.688 1.84398 28.7123 10.3562 -33.3441 12.6721 -5.73537 -2.26463 "
    "-9.31168 0.655842 4.9714 -3.0286 1.91266 -6.08734 -16.0596 4.02981";
    assert_vec(ref, softbits_approx);
  }
}
