/*!
 * \file
 * \brief Vector ("MIMO") modulator classes test program
 * \author Erik G. Larsson, Adam Piatyszek, and Mirsad Cirkic
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
#include <itpp/itstat.h>
#include "gtest/gtest.h"

using namespace itpp;

template<class T>
void assert_vec_p(const Vec<T> &ref, const Vec<T> &act, int line)
{
  static const double tol = 1e-3;
  ASSERT_EQ(ref.length(), ref.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol) << line;
  }
}
#define assert_vec(ref,act) assert_vec_p(ref, act, __LINE__)

TEST(ModulatorNd, All)
{
	// Test of ND (MIMO) Modulators
	RNG_reset(12345);
	double sigma2 = 0.05;

	{ // --- Test for M-PAM, M^2-QAM, and M^2-PSK constellations
		ND_UPAM pammod;
		ND_UQAM qammod;
		ND_UPSK pskmod;
		int nt = 2;
		for (int nb = 1; nb <= 2; nb++) {
			pammod.set_M(2*nt, 1<<nb);
			qammod.set_M(nt, 1<<2*nb);
      pskmod.set_M(nt, 1<<2*nb);

      bvec b=randb(2*nt*nb);
			cvec xqam = qammod.modulate_bits(b);
			cvec xpsk = pskmod.modulate_bits(b);
			cmat Hc = randn_c(nt, nt);
			cvec ec=sqrt(2*sigma2) * randn_c(nt);
			cvec yqam = Hc * xqam + ec;
			cvec ypsk = Hc * xpsk + ec;
			vec x=pammod.modulate_bits(b);
			mat H(2*nt,2*nt);
			H.set_submatrix(0, 0, real(Hc)/sqrt(2.0));
			H.set_submatrix(nt, 0, imag(Hc)/sqrt(2.0));
			H.set_submatrix(nt, nt, real(Hc)/sqrt(2.0));
			H.set_submatrix(0, nt, -imag(Hc)/sqrt(2.0));
			vec y=H*x+sqrt(sigma2) * randn(2*nt);
			
			QLLRvec LLR_ap = randi(2*nt*nb,-5000,5000);
			QLLRvec LLR;

			pammod.init_soft_demodulator(H, sigma2);
			qammod.init_soft_demodulator(Hc, 2*sigma2);
			pskmod.init_soft_demodulator(Hc, 2*sigma2);

			pammod.demodulate_soft_bits(y, LLR_ap, LLR);
      vec act = pammod.get_llrcalc().to_double(LLR);
      vec ref;
      switch (nb) {
          case 1:
          ref = "-61.2605 -21.5132 64.5134 -17.156";
          break;
          case 2:
          ref = "-9.55176 2.23926 -15.8718 2.30127 12.1917 2.27637 9.37671 -2.38306";
          break;
          default:
          ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			pammod.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR);
      act = pammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-61.2605 -21.5132 64.5134 -17.156";
        break;
        case 2:
        ref = "-9.55176 2.23926 -15.8718 2.30127 12.1917 2.27637 9.37671 -2.38306";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			
			qammod.demodulate_soft_bits(yqam, LLR_ap, LLR);
      act = qammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-26.7622 -31.0969 20.4077 -25.7842";
        break;
        case 2:
        ref = "-14.3232 3.77539 -38.3303 13.3435 44.6099 13.2417 3.77856 -11.3962";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			qammod.demodulate_soft_bits(yqam, Hc, 2*sigma2, LLR_ap, LLR);
      act = qammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-26.7622 -31.0969 20.4077 -25.7842";
        break;
        case 2:
        ref = "-14.3232 3.77539 -38.3303 13.3435 44.6099 13.2417 3.77856 -11.3962";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);

			pskmod.demodulate_soft_bits(ypsk, LLR_ap, LLR);
      act = pskmod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-28.6494 -28.8706 20.1487 -22.7725";
        break;
        case 2:
        ref = "-23.209 -2.04126 -5.97974 0.705566 5.83984 15.2581 2.95044 -0.491211";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			pskmod.demodulate_soft_bits(ypsk, Hc, 2*sigma2, LLR_ap, LLR);
      act = pskmod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:      
        ref = "-28.6494 -28.8706 20.1487 -22.7725";
        break;
        case 2:
        ref = "-23.209 -2.04126 -5.97974 0.705566 5.83984 15.2581 2.95044 -0.491211";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);

			pammod.demodulate_soft_bits(y, LLR_ap, LLR, ND_UPAM::FULL_ENUM_MAXLOG);
      act = pammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-61.2605 -21.5132 64.5134 -17.156";
        break;
        case 2:
        ref = "-9.38965 2.40625 -16.5156 2.302 12.1526 2.302 9.17822 -2.40625";
        break;        
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			LLR_calc_unit llrcalctmp=pammod.get_llrcalc();
			pammod.set_llrcalc(LLR_calc_unit(12, 0, 7));
			pammod.demodulate_soft_bits(y, LLR_ap, LLR);
      act = pammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-61.2605 -21.5132 64.5134 -17.156";
        break;
        case 2:
        ref = "-9.38965 2.40625 -16.5156 2.302 12.1526 2.302 9.17822 -2.40625";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);

			qammod.demodulate_soft_bits(yqam, LLR_ap, LLR, ND_UQAM::FULL_ENUM_MAXLOG);
      act = qammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-26.7622 -31.0969 20.4077 -26.2703";
        break;
        case 2:
        ref = "-14.3618 3.79297 -38.3418 13.7219 44.8179 13.626 3.79297 -11.3765";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);

			pskmod.demodulate_soft_bits(ypsk, LLR_ap, LLR, ND_UPSK::FULL_ENUM_MAXLOG);
      act = pskmod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-28.6494 -28.8706 20.1489 -22.7754";
        break;
        case 2:
        ref = "-22.5317 -1.95972 -5.72705 0.475098 5.49121 14.6724 2.47705 -0.475098";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);

			pammod.set_llrcalc(llrcalctmp);
			pammod.demodulate_soft_bits(y, diag(H), sigma2, LLR_ap, LLR);
      act = pammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "1.36011 6.94312 -0.0161133 -12.0933";
        break;
        case 2:
        ref = "2.62354 1.00073 -17.6406 8.08691 -1.66406 0.052002 2.45728 0.338379";
        break;        
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      assert_vec(ref, act);
			pammod.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::ZF_LOGMAP);
      act = pammod.get_llrcalc().to_double(LLR);
      switch (nb) {
        case 1:
        ref = "-19.4644 -7.35938 18.3091 -6.09521";
        break;
        case 2:
        ref = "-8.52222 1.0542 -15.8455 3.78247 12.3943 3.12109 5.45532 -2.86353";
        break;
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
      ivec zhat;
      ivec ref_i;
      pammod.sphere_decoding(y, H, 0.01, 10000, 2.0, zhat);
      switch (nb) {
        case 1:
        ref_i = "-1000 -1000 1000 -1000";
        break;
        case 2:
        ref_i = "-1000 1000 -1000 1000 1000 1000 1000 -1000";
        break;        
        default:
        ASSERT_TRUE(false) << "should not be here";
      }
			assert_vec(ref_i, zhat);
		}
	}

	{ // --- Test for M-PAM and M^2-QAM with variable constellations
		ND_UQAM qammod;
		ND_UPAM pammod;
		int nt = 3;
		
 		qammod.set_M(nt, "4 64 16");
		pammod.set_M(2*nt, "2 8 4 2 8 4");
		cmat Hc = randn_c(nt, nt);
		mat H(2*nt,2*nt);
		H.set_submatrix(0, 0, real(Hc)/sqrt(2.0));
		H.set_submatrix(nt, 0, imag(Hc)/sqrt(2.0));
		H.set_submatrix(nt, nt, real(Hc)/sqrt(2.0));
		H.set_submatrix(0, nt, -imag(Hc)/sqrt(2.0));
		QLLRvec LLR_ap = zeros_i(sum(qammod.get_k()));
		pammod.init_soft_demodulator(H, sigma2);
		qammod.init_soft_demodulator(Hc, 2*sigma2);	
		double pamerrs=0, qamerrs=0;
		int TRIALS=1000;
		for(int i=0; i<TRIALS; i++){
			bvec b=randb(sum(pammod.get_k()));
			cvec xqam = qammod.modulate_bits(b);
			cvec yqam = Hc * xqam + sqrt(2*sigma2) * randn_c(nt);
			vec x=pammod.modulate_bits(b);
			vec y=H*x+sqrt(sigma2) * randn(2*nt);
			
			QLLRvec LLR;
			pammod.demodulate_soft_bits(y, LLR_ap, LLR);
			pamerrs+=hamming_distance(LLR<0,b);
			qammod.demodulate_soft_bits(yqam, LLR_ap, LLR);
			qamerrs+=hamming_distance(LLR<0,b);
		}
		pamerrs/=TRIALS*sum(pammod.get_k());
		qamerrs/=TRIALS*sum(pammod.get_k());
		ASSERT_EQ(round_i(100*qamerrs), round_i(100*pamerrs));
		bvec b=randb(sum(pammod.get_k()));		
		cvec xqam = qammod.modulate_bits(b);
		cvec yqam = Hc * xqam + sqrt(2*sigma2) * randn_c(nt);
		vec x=pammod.modulate_bits(b);
		vec y=H*x+sqrt(sigma2) * randn(2*nt);
		QLLRvec LLR;
    pammod.demodulate_soft_bits(y, LLR_ap, LLR);
    vec act = pammod.get_llrcalc().to_double(LLR);
    vec ref = "34.6016 -2.96387 -30.106 11.395 16.1267 -3.16919 29.9153 50.3799 4.47412 -8.37354 23.3008 3.49634";
    assert_vec(ref, act);
		qammod.demodulate_soft_bits(yqam, LLR_ap, LLR);
    act = qammod.get_llrcalc().to_double(LLR);
    ref = "60.292 -58.8679 -98.0229 22.6648 6.70361 -107.042 22.8962 3.21558 12.7021 -7.41406 34.2944 6.74219";
    assert_vec(ref, act);
	}
}
