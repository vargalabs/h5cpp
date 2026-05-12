/*!
 * \file
 * \brief test program for MIMO demapper
 * \author Bogdan Cristea
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

#include "itpp/itcomm.h"
#include <string>

using namespace std;
using namespace itpp;

int main()
{
    //parameters
    int const_size = 16;
    string demapper_method[] = {"Hassibi_maxlogMAP", "GA", "sGA", "mmsePIC", "zfPIC"};
    string code_name = "V-BLAST_MxN";
    int em_antennas = 2;
    int rec_antennas = 2;
    int channel_uses = 1;
    int perm_len = pow2i(6);//permutation length

    //QAM modulator class
    QAM mod(const_size);

    //Space-Time code parameters
    STC st_block_code(code_name, const_size, em_antennas, channel_uses);//generate matrices for LD code (following Hassibi's approach)

    //SISO blocks
    SISO siso;
    siso.set_constellation(mod.bits_per_symbol(), mod.get_symbols(),
    		mod.get_bits2symbols());
    siso.set_st_block_code(st_block_code.get_nb_symbols_per_block(),
    		st_block_code.get_1st_gen_matrix(), st_block_code.get_2nd_gen_matrix(),
    		rec_antennas);
    siso.set_noise(1e-1);

    //bits generation
    bvec bits = randb(perm_len);

    //QAM modulation
    cvec em = mod.modulate_bits(bits)/sqrt(double(em_antennas));//normalize emitted symbols

    //ST code
    cmat S = st_block_code.encode(em);

    //internal variables
    int symb_block = st_block_code.get_nb_symbols_per_block();
    int nb_symb = perm_len/mod.bits_per_symbol();//number of symbols at the modulator output
    int nb_subblocks = nb_symb/symb_block;//number of blocks of ST code emitted in an interleaver period
    int tx_duration = channel_uses*nb_subblocks;//transmission duration expressed in number of symbol periods

    //ideal channel
    if (em_antennas != rec_antennas)
    {
	    return EXIT_FAILURE;
    }
    cmat mimo_channel = eye_c(em_antennas);
    cmat rec(tx_duration,rec_antennas);
    for (int ns=0;ns<nb_subblocks;ns++)
    {
    	rec.set_submatrix(ns*channel_uses, 0,
    			S(ns*channel_uses, (ns+1)*channel_uses-1, 0, em_antennas-1)*
			mimo_channel);
    }    

    //first decoder
    vec demapper_apriori_data(perm_len);
    vec demapper_extrinsic_data(perm_len);
    demapper_apriori_data.zeros();//a priori information of emitted bits
    cmat ch_att(em_antennas*rec_antennas, tx_duration);
    ch_att = kron(reshape(mimo_channel, em_antennas*rec_antennas, 1), ones_c(1, tx_duration));
    siso.set_impulse_response(ch_att);
    BERC ber;
    for (unsigned int n = 0; n < sizeof(demapper_method)/sizeof(demapper_method[0]); ++n)
    {
        siso.set_demapper_method(demapper_method[n]);
        siso.demapper(demapper_extrinsic_data, rec, demapper_apriori_data);

        //show results
        ber.count((demapper_extrinsic_data > 0), bits);
        cout << demapper_method[n] << ", BER = " << ber.get_errorrate() << endl;
	demapper_extrinsic_data.zeros();
    }
}

