/*!
 * \file
 * \brief Vector ("MIMO") modulator classes test program
 * \author Erik G. Larsson, Adam Piatyszek, and Mirsad Cirkic
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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
#include <iomanip>

using namespace std;
using namespace itpp;


int main()
{
	cout << "========================================================" << endl;
	cout << "               Test of ND (MIMO) Modulators             " << endl;
	cout << "========================================================" << endl;

	cout.setf(ios::fixed);
	cout.precision(2);

	RNG_reset(12345);
	double sigma2 = 0.05;

	{ // --- Test for M-PAM, M^2-QAM, and M^2-PSK constellations
		ND_UPAM pammod;
		ND_UQAM qammod;
		ND_UPSK pskmod;
		int nt = 2;
		for (int nb = 1; nb <= 2; nb++) {
			cout << "================== " << (1 << nb) << "-PAM, " << (1<<2*nb) << "-QAM, and " << (1<<2*nb) << "-PSK ==================\n";
			pammod.set_M(2*nt, 1<<nb);
			qammod.set_M(nt, 1<<2*nb);
			pskmod.set_M(nt, 1<<2*nb);
			cout << pammod << endl;
			cout << qammod << endl;
			cout << pskmod << endl;
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
			cout << "PAM: full LLR: " << pammod.get_llrcalc().to_double(LLR) << endl;
			pammod.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR);
			cout << "             : " << pammod.get_llrcalc().to_double(LLR) << endl;
			
			qammod.demodulate_soft_bits(yqam, LLR_ap, LLR);
			cout << "QAM: full LLR: " << qammod.get_llrcalc().to_double(LLR) << endl;
			qammod.demodulate_soft_bits(yqam, Hc, 2*sigma2, LLR_ap, LLR);
			cout << "             : " << qammod.get_llrcalc().to_double(LLR) << endl;

			pskmod.demodulate_soft_bits(ypsk, LLR_ap, LLR);
			cout << "PSK: full LLR: " << pskmod.get_llrcalc().to_double(LLR) << endl;
			pskmod.demodulate_soft_bits(ypsk, Hc, 2*sigma2, LLR_ap, LLR);
			cout << "             : " << pskmod.get_llrcalc().to_double(LLR) << endl;

			pammod.demodulate_soft_bits(y, LLR_ap, LLR, ND_UPAM::FULL_ENUM_MAXLOG);
			cout << "PAM: Max-Log : " << pammod.get_llrcalc().to_double(LLR) << endl;
			LLR_calc_unit llrcalctmp=pammod.get_llrcalc();
			pammod.set_llrcalc(LLR_calc_unit(12, 0, 7));
			pammod.demodulate_soft_bits(y, LLR_ap, LLR);
			cout << "             : " << pammod.get_llrcalc().to_double(LLR) << endl;

			qammod.demodulate_soft_bits(yqam, LLR_ap, LLR, ND_UQAM::FULL_ENUM_MAXLOG);
			cout << "QAM: Max-Log : " << qammod.get_llrcalc().to_double(LLR) << endl;

			pskmod.demodulate_soft_bits(ypsk, LLR_ap, LLR, ND_UPSK::FULL_ENUM_MAXLOG);
			cout << "PSK: Max-Log : " << pskmod.get_llrcalc().to_double(LLR) << endl << endl;

			pammod.set_llrcalc(llrcalctmp);
			pammod.demodulate_soft_bits(y, diag(H), sigma2, LLR_ap, LLR);
			cout << "PAM: diag. model : " << pammod.get_llrcalc().to_double(LLR) << endl;
			pammod.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::ZF_LOGMAP);
			cout << "PAM: zero-forcing: " << pammod.get_llrcalc().to_double(LLR) << endl;
			ivec zhat;
			pammod.sphere_decoding(y, H, 0.01, 10000, 2.0, zhat);
			cout << "PAM: sphere-dec. : " << zhat << endl << endl;
		}
	}

	{ // --- Test for M-PAM and M^2-QAM with variable constellations
		ND_UQAM qammod;
		ND_UPAM pammod;
		int nt = 3;
		
 		qammod.set_M(nt, "4 64 16");
		pammod.set_M(2*nt, "2 8 4 2 8 4");
		cout << qammod << endl << pammod << endl;
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
		cout << itpp::round(100*pamerrs)/100 << " " << itpp::round(100*qamerrs)/100 << endl;
		if(round_i(100*qamerrs)==round_i(100*pamerrs)) cout << "BER passed" << endl;
		else cout << "BER failed" << endl;			
		bvec b=randb(sum(pammod.get_k()));		
		cvec xqam = qammod.modulate_bits(b);
		cvec yqam = Hc * xqam + sqrt(2*sigma2) * randn_c(nt);
		vec x=pammod.modulate_bits(b);
		vec y=H*x+sqrt(sigma2) * randn(2*nt);
		QLLRvec LLR;
		pammod.demodulate_soft_bits(y, LLR_ap, LLR);
		cout << "bit sequence : " << b << endl;
		cout << "PAM: full LLR: " << pammod.get_llrcalc().to_double(LLR) << endl;		
		qammod.demodulate_soft_bits(yqam, LLR_ap, LLR);
		cout << "QAM: full LLR: " << qammod.get_llrcalc().to_double(LLR) << endl;
	}

	return 0;
}
