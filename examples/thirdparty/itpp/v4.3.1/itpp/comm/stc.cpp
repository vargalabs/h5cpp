/*!
 * \file
 * \brief Implementation of Space Time block Codes (STC) class
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2011  (see AUTHORS file for a list of contributors)
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

#include <itpp/comm/stc.h>

namespace itpp
{

void STC::Hassibi_block_code(void)
/* This function generates the A and B matrices needed for Space-Time block codes
 * generation following Hassibi's approach:
 * S = sum_{q=1}^symb_block (A_q alpha_q + jB_q beta_q),
 * where s_q = alpha_q+jbeta_q is the symbol after modulation
 * each A_q and B_q matrix has dimension TxM
 * different A_q and B_q matrices are stacked one below the other, e.g.
 * [A_1;A_2;...;A_Q]

 * input: code_name - code name whose generator matrices are to be generated
 *        const_size - constellation size (used in Damen code)
 * outputs: symb_block - number of symbols per block
 *         A, B - generator matrices
 * inputs/outputs: for some codes these are inputs for others they are
 * predefined, so they are outputs only
 *               em_antennas - number of emission antenna
 *               channel_uses - channel uses
 */
{
    if (code_name==Code_Names::V_BLAST_MxN)//classical V-BLAST
    {
    	it_assert(channel_uses > 0, "Channel uses should be strictly positive");
    	it_assert(em_antenna > 0, "Number of emission antenna should be strictly positive");
        symb_block = channel_uses*em_antenna;//number of symbols/block
        A.set_size(symb_block*channel_uses, em_antenna);
        A.zeros();
        itpp::mat temp(channel_uses, em_antenna);
        temp.zeros();
        register int tau,m;
        for (tau=0; tau<channel_uses; tau++)
        {
            for (m=0; m<em_antenna; m++)
            {
                temp(tau,m) = 1;
                A.set_submatrix(symb_block*tau+channel_uses*m, 0, itpp::to_cmat(temp));
                temp(tau,m) = 0;
            }
        }
        B = A;
    }
    else if (code_name == Code_Names::imp_V_BLAST_MxN)//improved V-BLAST (code (31) in Hassibi's paper)
    {
    	it_assert(em_antenna > 0, "Number of emission antenna should be strictly positive");
    	it_assert(channel_uses == em_antenna, "Channel uses and the number of emission antenna must be equal");
        symb_block = channel_uses*em_antenna;//number of symbols/block
        std::complex<double> j(0,1);
        itpp::cmat D = itpp::diag(exp(j*(2*itpp::pi/em_antenna)*
        		itpp::linspace(0, em_antenna-1, em_antenna)));
        itpp::mat P = itpp::diag(itpp::ones(em_antenna-1), -1);
        P(0,em_antenna-1) = 1;
        A.set_size(symb_block*channel_uses, em_antenna);
        A.zeros();
        register int k,l;
        for (k=0; k<channel_uses; k++)
        {
            for (l=0; l<em_antenna; l++)
            {
                A.set_submatrix(symb_block*k+l*channel_uses, 0,
                		diag_pow(D, k)*itpp::to_cmat(mat_pow(P, l))/
                		std::sqrt(double(em_antenna)));
            }
        }
        B = A;
    }
    else if (code_name==Code_Names::Alamouti_2xN)//Alamouti's orthogonal code
    {
        em_antenna = 2;//emission antenna
        channel_uses = 2;//channel uses
        symb_block = 2;//number of symbols/block

        A = "1  0;"
            "0  1;"
            "0  1;"
        	"-1  0";//A_1; A_2
        B = "1  0;"
            "0 -1;"
            "0  1;"
            "1  0";//B_1; B_2
    }
    else if (code_name==Code_Names::Switched_Alamouti_4xN)
    {
        em_antenna = 4;//emission antenna
        channel_uses = 4;//channel uses
        symb_block = 4;//number of symbols/block

        A = "1  0  0  0;"
            "0  1  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  1  0  0;"
            "-1 0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  1  0;"
            "0  0  0  1;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  1;"
            "0  0  -1 0";//A_1; A_2; A_3; A_4
        A *= std::sqrt(2.0);//normalization
        B = "1  0  0  0;"
            "0  -1 0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  1  0  0;"
            "1  0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  1  0;"
            "0  0  0  -1;"
            "0  0  0  0;"
            "0  0  0  0;"
            "0  0  0  1;"
            "0  0  1  0";//B_1; B_2; B_3; B_4
        B *= std::sqrt(2.0);
    }
    else if (code_name==Code_Names::Double_Alamouti_4xN)
    {
        em_antenna = 4;//emission antenna
        channel_uses = 2;//channel uses
        symb_block = 4;//number of symbols/block

        A = "1  0  0  0;"
            "0  1  0  0;"
            "0  0  1  0;"
            "0  0  0  1;"
            "0  1  0  0;"
            "-1 0  0  0;"
            "0  0  0  1;"
            "0  0  -1 0";//A_1; A_2; A_3; A_4
        B = "1  0  0  0;"
            "0  -1 0  0;"
            "0  0  1  0;"
            "0  0  0 -1;"
            "0  1  0  0;"
            "1  0  0  0;"
            "0  0  0  1;"
            "0  0  1  0";//B_1; B_2; B_3; B_4
    }
    else if (code_name== Code_Names::Jafarkhani_4xN)//Jafarkhani's quasi-orthogonal code
    {
        em_antenna = 4;//emission antenna
        channel_uses = 4;//channel uses
        symb_block = 4;//number of symbols/block

        A = "1  0  0  0;"
            "0  1  0  0;"
            "0  0  1  0;"
            "0  0  0  1;"
            "0  1  0  0;"
            "-1 0  0  0;"
            "0  0  0  1;"
            "0  0  -1 0;"
            "0  0  1  0;"
            "0  0  0  1;"
            "-1 0  0  0;"
            "0  -1 0  0;"
            "0  0  0  1;"
            "0  0  -1 0;"
            "0  -1 0  0;"
            "1  0  0  0";//A_1; A_2; A_3; A_4
        B = "1  0  0  0;"
            "0  -1  0  0;"
            "0  0  -1  0;"
            "0  0  0  1;"
            "0  1  0  0;"
            "1  0  0  0;"
            "0  0  0  -1;"
            "0  0  -1 0;"
            "0  0  1  0;"
            "0  0  0  -1;"
            "1  0  0  0;"
            "0  -1 0  0;"
            "0  0  0  1;"
            "0  0  1  0;"
            "0  1  0  0;"
            "1  0  0  0";//B_1; B_2; B_3; B_4
    }
    else if (code_name== Code_Names::Golden_2x2)//Golden code as proposed by Belfiore
    {
        em_antenna = 2;//emission antenna
        channel_uses = 2;//channel uses
        symb_block = 4;//number of symbols/block

        std::complex<double> theta((1+std::sqrt(5.0))/2,0);
        std::complex<double> theta_b((1-std::sqrt(5.0))/2,0);
        std::complex<double> j(0,1);
        std::complex<double> one(1,0);
        std::complex<double> alpha = one+j*(one-theta);
        std::complex<double> alpha_b = one+j*(one-theta_b);
        std::complex<double> gamma = j;
        A.set_size(8,2);
        A(0,0) = alpha/std::sqrt(5.0);
        A(0,1) = 0;
        A(1,0) = 0;
        A(1,1) =  alpha_b/std::sqrt(5.0);//A_1
        A(2,0) = alpha*theta/std::sqrt(5.0);
        A(2,1) = 0;
        A(3,0) = 0;
        A(3,1) = alpha_b*theta_b/std::sqrt(5.0);//A_2
        A(4,0) = 0;
        A(4,1) = gamma*alpha_b/std::sqrt(5.0);
        A(5,0) = alpha/std::sqrt(5.0);
        A(5,1) = 0;//A_3
        A(6,0) = 0;
        A(6,1) = gamma*alpha_b*theta_b/std::sqrt(5.0);
        A(7,0) = alpha*theta/std::sqrt(5.0);
        A(7,1) = 0;//A_4
        B = A;
    }
    else if (code_name== Code_Names::Damen_2x2)//ST code based on number theory as proposed by Damen
    {
        em_antenna = 2;//emission antenna
        channel_uses = 2;//channel uses
        symb_block = 4;//number of symbols/block

        double lambda;
        if (const_size==4)
            lambda = 0.5;
        else if (const_size==16)
            lambda = 0.521;
        else if (const_size>=256)
            lambda = itpp::pi/4;
        else
        {
            lambda = itpp::pi/4;
            std::cout << "STC::LDcode: Warning! For " << string_from_code_name(code_name) <<
            		" and const. size " << const_size << ", lambda has the "
            		"value " << lambda << std::endl;
        }
        std::complex<double> j(0,1);
        std::complex<double> phi = std::exp(j*lambda);
        std::complex<double> theta = std::exp(j*(lambda/2));
        A.set_size(8, 2);
        A(0,0) = 1/std::sqrt(2.0);
        A(0,1) = 0;
        A(1,0) = 0;
        A(1,1) = 1/std::sqrt(2.0);//A_1
        A(2,0) = phi/std::sqrt(2.0);
        A(2,1) = 0;
        A(3,0) = 0;
        A(3,1) = -phi/std::sqrt(2.0);//A_2
        A(4,0) = 0;
        A(4,1) = theta/std::sqrt(2.0);
        A(5,0) = theta/std::sqrt(2.0);
        A(5,1) = 0;//A_3
        A(6,0) = 0;
        A(6,1) = -theta*phi/std::sqrt(2.0);
        A(7,0) = theta*phi/std::sqrt(2.0);
        A(7,1) = 0;//A_4
        B = A;
    }
    else if (code_name== Code_Names::ortho34_3xN)//rate 3/4 orthogonal code (mutual information 5.13 bits/channel use at rho=20 dB)
    {
        em_antenna = 3;//emission antenna
        channel_uses = 4;//channel uses
        symb_block = 3;//number of symbols/block

        A = "1 0 0;"
            "0 1 0;"
            "0 0 1;"
            "0 0 0;"
            "0 1 0;"
            "-1 0 0;"
            "0 0 0;"
            "0 0 1;"
            "0 0 1;"
            "0 0 0;"
            "-1 0 0;"
            "0 -1 0";//A_1; A_2; A_3
        A /= std::sqrt(double(4)/double(3));
        B = "1 0 0;"
            "0 -1 0;"
            "0 0 -1;"
            "0 0 0;"
            "0 1 0;"
            "1 0 0;"
            "0 0 0;"
            "0 0 -1;"
            "0 0 1;"
            "0 0 0;"
            "1 0 0;"
            "0 1 0";//B_1; B_2; B_3
        B /= std::sqrt(double(4)/double(3));
    }
    else if (code_name== Code_Names::LD36_3xN)//(36) LD code with mutual info. 6.25bits/channel use at rho=20dB
    {
        em_antenna = 3;//emission antenna
        channel_uses = 4;//channel uses
        symb_block = 4;//number of symbols/block

        A.set_size(16, 3);
        A(0,0) = 1;
        A(0,1) = 0;
        A(0,2) = 0;
        A(1,0) = 1;
        A(1,1) = 1;
        A(1,2) = 0;
        A(2,0) = 0;
        A(2,1) = 0;
        A(2,2) = 1;
        A(3,0) = 0;
        A(3,1) = 0;
        A(3,2) = 0;//A_1
        A(4,0) = 0;
        A(4,1) = 1/std::sqrt(2.0);
        A(4,2) = 0;
        A(5,0) = -1/std::sqrt(2.0);
        A(5,1) = 0;
        A(5,2) = -1/std::sqrt(2.0);
        A(6,0) = 0;
        A(6,1) = 1/std::sqrt(2.0);
        A(6,2) = 0;
        A(7,0) = 1/std::sqrt(2.0);
        A(7,1) = 0;
        A(7,2) = -1/std::sqrt(2.0);//A_2
        A(8,0) = 1;
        A(8,1) = 0;
        A(8,2) = 0;
        A(9,0) = 0;
        A(9,1) = 0;
        A(9,2) = 0;
        A(10,0) = 0;
        A(10,1) = 0;
        A(10,2) = -1;
        A(11,0) = 0;
        A(11,1) = -1;
        A(11,2) = 0;//A_3
        A(12,0) = 0;
        A(12,1) = -1/std::sqrt(2.0);
        A(12,2) = 0;
        A(13,0) = 1/std::sqrt(2.0);
        A(13,1) = 0;
        A(13,2) = -1/std::sqrt(2.0);
        A(14,0) = 0;
        A(14,1) = 1/std::sqrt(2.0);
        A(14,2) = 0;
        A(15,0) = -1/std::sqrt(2.0);
        A(15,1) = 0;
        A(15,2) = -1/std::sqrt(2.0);//A_4
        B.set_size(16, 3);
        B(0,0) = 0;
        B(0,1) = -1/std::sqrt(2.0);
        B(0,2) = 0;
        B(1,0) = -1/std::sqrt(2.0);
        B(1,1) = 0;
        B(1,2) = 1/std::sqrt(2.0);
        B(2,0) = 0;
        B(2,1) = 1/std::sqrt(2.0);
        B(2,2) = 0;
        B(3,0) = 1/std::sqrt(2.0);
        B(3,1) = 0;
        B(3,2) = 1/std::sqrt(2.0);//B_1
        B(4,0) = 1/std::sqrt(2.0);
        B(4,1) = double(-1)/double(2);
        B(4,2) = 0;
        B(5,0) = double(-1)/double(2);
        B(5,1) = -1/std::sqrt(2.0);
        B(5,2) = double(-1)/double(2);
        B(6,0) = 0;
        B(6,1) = double(-1)/double(2);
        B(6,2) = 1/std::sqrt(2.0);
        B(7,0) = double(1)/double(2);
        B(7,1) = 0;
        B(7,2) = double(-1)/double(2);//B_2
        B(8,0) = 1/std::sqrt(2.0);
        B(8,1) = double(1)/double(2);
        B(8,2) = 0;
        B(9,0) = double(1)/double(2);
        B(9,1) = -1/std::sqrt(2.0);
        B(9,2) = double(1)/double(2);
        B(10,0) = 0;
        B(10,1) = double(1)/double(2);
        B(10,2) = 1/std::sqrt(2.0);
        B(11,0) = double(-1)/double(2);
        B(11,1) = 0;
        B(11,2) = double(1)/double(2);//B_3
        B(12,0) = 1;
        B(12,1) = 0;
        B(12,2) = 0;
        B(13,0) = 0;
        B(13,1) = 0;
        B(13,2) = 0;
        B(14,0) = 0;
        B(14,1) = 0;
        B(14,2) = -1;
        B(15,0) = 0;
        B(15,1) = 1;
        B(15,2) = 0;//B_4
    }
    else if (code_name== Code_Names::LD37_3xN)//(37) LD code 3-antenna LD code obtained from the symetrical concatenation of 3 2-antenna orthogonal design
    {
        em_antenna = 3;//emission antenna
        channel_uses = 6;//channel uses
        symb_block = 6;//number of symbols/block

        A = "1  0  0;"
            "0  1  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  1  0;"
            "-1 0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  1  0;"
            "0  0  1;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  1;"
            "0 -1  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "1  0  0;"
            "0  0  1;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  1;"
            "-1 0  0";//A_1; A_2; A_3; A_4; A_5; A_6
        A *= std::sqrt(double(3)/double(2));
        B = "1  0  0;"
            "0  -1  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  1  0;"
            "1 0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  1  0;"
            "0  0  -1;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  1;"
            "0  1  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "1  0  0;"
            "0  0  -1;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  0;"
            "0  0  1;"
            "1  0  0";//B_1; B_2; B_3; B_4; B_5; B_6
        B *= std::sqrt(double(3)/double(2));
    }
    else if (code_name== Code_Names::LD39_3xN)
    {
        em_antenna = 3;//emission antenna
        channel_uses = 6;//channel uses
        symb_block = 6;//number of symbols/block

        A.set_size(36, 3);
        A(0,0) = 1/std::sqrt(2.0);
        A(0,1) = 0;
        A(0,2) = 0;
        A(1,0) = 0;
        A(1,1) = 1/std::sqrt(2.0);
        A(1,2) = 0;
        A(2,0) = 0;
        A(2,1) = 1/std::sqrt(2.0);
        A(2,2) = 0;
        A(3,0) = 0;
        A(3,1) = 0;
        A(3,2) = 1/std::sqrt(2.0);
        A(4,0) = 1/std::sqrt(2.0);
        A(4,1) = 0;
        A(4,2) = 0;
        A(5,0) = 0;
        A(5,1) = 0;
        A(5,2) = 1/std::sqrt(2.0);//A_1
        A(6,0) = 0;
        A(6,1) = 1/std::sqrt(2.0);
        A(6,2) = 0;
        A(7,0) = -1/std::sqrt(2.0);
        A(7,1) = 0;
        A(7,2) = 0;
        A(8,0) = 0;
        A(8,1) = 0;
        A(8,2) = 1/std::sqrt(2.0);
        A(9,0) = 0;
        A(9,1) = -1/std::sqrt(2.0);
        A(9,2) = 0;
        A(10,0) = 0;
        A(10,1) = 0;
        A(10,2) = 1/std::sqrt(2.0);
        A(11,0) = -1/std::sqrt(2.0);
        A(11,1) = 0;
        A(11,2) = 0;//A_2
        A(12,0) = 1/std::sqrt(2.0);
        A(12,1) = 0;
        A(12,2) = 0;
        A(13,0) = 0;
        A(13,1) = 1/std::sqrt(2.0);
        A(13,2) = 0;
        A(14,0) = 0;
        A(14,1) = -1/(2*std::sqrt(2.0));
        A(14,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(15,0) = 0;
        A(15,1) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(15,2) = -1/(2*std::sqrt(2.0));
        A(16,0) = -1/(2*std::sqrt(2.0));
        A(16,1) = 0;
        A(16,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(17,0) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(17,1) = 0;
        A(17,2) = -1/(2*std::sqrt(2.0));//A_3
        A(18,0) = 0;
        A(18,1) = 1/std::sqrt(2.0);
        A(18,2) = 0;
        A(19,0) = -1/std::sqrt(2.0);
        A(19,1) = 0;
        A(19,2) = 0;
        A(20,0) = 0;
        A(20,1) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(20,2) = -1/(2*std::sqrt(2.0));
        A(21,0) = 0;
        A(21,1) = 1/(2*std::sqrt(2.0));
        A(21,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(22,0) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(22,1) = 0;
        A(22,2) = -1/(2*std::sqrt(2.0));
        A(23,0) = 1/(2*std::sqrt(2.0));
        A(23,1) = 0;
        A(23,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));//A_4
        A(24,0) = 1/std::sqrt(2.0);
        A(24,1) = 0;
        A(24,2) = 0;
        A(25,0) = 0;
        A(25,1) = 1/std::sqrt(2.0);
        A(25,2) = 0;
        A(26,0) = 0;
        A(26,1) = -1/(2*std::sqrt(2.0));
        A(26,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(27,0) = 0;
        A(27,1) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(27,2) = -1/(2*std::sqrt(2.0));
        A(28,0) = -1/(2*std::sqrt(2.0));
        A(28,1) = 0;
        A(28,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(29,0) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(29,1) = 0;
        A(29,2) = -1/(2*std::sqrt(2.0));//A_5
        A(30,0) = 0;
        A(30,1) = 1/std::sqrt(2.0);
        A(30,2) = 0;
        A(31,0) = -1/std::sqrt(2.0);
        A(31,1) = 0;
        A(31,2) = 0;
        A(32,0) = 0;
        A(32,1) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(32,2) = -1/(2*std::sqrt(2.0));
        A(33,0) = 0;
        A(33,1) = 1/(2*std::sqrt(2.0));
        A(33,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(34,0) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        A(34,1) = 0;
        A(34,2) = -1/(2*std::sqrt(2.0));
        A(35,0) = 1/(2*std::sqrt(2.0));
        A(35,1) = 0;
        A(35,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));//A_6
        B.set_size(36, 3);
        B(0,0) = 1/std::sqrt(2.0);
        B(0,1) = 0;
        B(0,2) = 0;
        B(1,0) = 0;
        B(1,1) = -1/std::sqrt(2.0);
        B(1,2) = 0;
        B(2,0) = 0;
        B(2,1) = 1/std::sqrt(2.0);
        B(2,2) = 0;
        B(3,0) = 0;
        B(3,1) = 0;
        B(3,2) = -1/std::sqrt(2.0);
        B(4,0) = 1/std::sqrt(2.0);
        B(4,1) = 0;
        B(4,2) = 0;
        B(5,0) = 0;
        B(5,1) = 0;
        B(5,2) = -1/std::sqrt(2.0);//B_1
        B(6,0) = 0;
        B(6,1) = 1/std::sqrt(2.0);
        B(6,2) = 0;
        B(7,0) = 1/std::sqrt(2.0);
        B(7,1) = 0;
        B(7,2) = 0;
        B(8,0) = 0;
        B(8,1) = 0;
        B(8,2) = 1/std::sqrt(2.0);
        B(9,0) = 0;
        B(9,1) = 1/std::sqrt(2.0);
        B(9,2) = 0;
        B(10,0) = 0;
        B(10,1) = 0;
        B(10,2) = 1/std::sqrt(2.0);
        B(11,0) = 1/std::sqrt(2.0);
        B(11,1) = 0;
        B(11,2) = 0;//B_2
        B(12,0) = 1/std::sqrt(2.0);
        B(12,1) = 0;
        B(12,2) = 0;
        B(13,0) = 0;
        B(13,1) = -1/std::sqrt(2.0);
        B(13,2) = 0;
        B(14,0) = 0;
        B(14,1) = -1/(2*std::sqrt(2.0));
        B(14,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(15,0) = 0;
        B(15,1) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(15,2) = 1/(2*std::sqrt(2.0));
        B(16,0) = -1/(2*std::sqrt(2.0));
        B(16,1) = 0;
        B(16,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(17,0) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(17,1) = 0;
        B(17,2) = 1/(2*std::sqrt(2.0));//B_3
        B(18,0) = 0;
        B(18,1) = 1/std::sqrt(2.0);
        B(18,2) = 0;
        B(19,0) = 1/std::sqrt(2.0);
        B(19,1) = 0;
        B(19,2) = 0;
        B(20,0) = 0;
        B(20,1) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(20,2) = -1/(2*std::sqrt(2.0));
        B(21,0) = 0;
        B(21,1) = -1/(2*std::sqrt(2.0));
        B(21,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(22,0) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(22,1) = 0;
        B(22,2) = -1/(2*std::sqrt(2.0));
        B(23,0) = -1/(2*std::sqrt(2.0));
        B(23,1) = 0;
        B(23,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));//B_4
        B(24,0) = 1/std::sqrt(2.0);
        B(24,1) = 0;
        B(24,2) = 0;
        B(25,0) = 0;
        B(25,1) = -1/std::sqrt(2.0);
        B(25,2) = 0;
        B(26,0) = 0;
        B(26,1) = -1/(2*std::sqrt(2.0));
        B(26,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(27,0) = 0;
        B(27,1) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(27,2) = 1/(2*std::sqrt(2.0));
        B(28,0) = -1/(2*std::sqrt(2.0));
        B(28,1) = 0;
        B(28,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(29,0) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(29,1) = 0;
        B(29,2) = 1/(2*std::sqrt(2.0));//B_5
        B(30,0) = 0;
        B(30,1) = 1/std::sqrt(2.0);
        B(30,2) = 0;
        B(31,0) = 1/std::sqrt(2.0);
        B(31,1) = 0;
        B(31,2) = 0;
        B(32,0) = 0;
        B(32,1) = -std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(32,2) = -1/(2*std::sqrt(2.0));
        B(33,0) = 0;
        B(33,1) = -1/(2*std::sqrt(2.0));
        B(33,2) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(34,0) = std::sqrt(3.0)/(2*std::sqrt(2.0));
        B(34,1) = 0;
        B(34,2) = -1/(2*std::sqrt(2.0));
        B(35,0) = -1/(2*std::sqrt(2.0));
        B(35,1) = 0;
        B(35,2) = -std::sqrt(3.0)/(2*std::sqrt(2.0));//B_6
    }
    else
    {
        it_assert(false, "Unknown code name.");
    }
}

itpp::cmat STC::encode(const itpp::cvec &symb)
//LD code generation (symb_block symbols go to an channel_uses x em_antennas matrix) following Hassibi's approach
{
    int nb_subblocks = symb.length()/symb_block;
    int tx_duration = channel_uses*nb_subblocks;
    itpp::cmat S(tx_duration,em_antenna);
    itpp::cmat temp(channel_uses,em_antenna);
    std::complex<double> j(0,1);
    register int ns,k;
    for (ns=0; ns<nb_subblocks; ns++)//encode block by block (symb_block symbols)
    {
        temp.zeros();
        for (k=0; k<symb_block; k++)//sum over all symb_block matrices
        {
            temp += (A(k*channel_uses,(k+1)*channel_uses-1,0,em_antenna-1)*
            		static_cast< std::complex<double> >(symb(k+ns*symb_block).real())+
                    j*B(k*channel_uses,(k+1)*channel_uses-1,0,em_antenna-1)*
                    static_cast< std::complex<double> >(symb(k+ns*symb_block).imag()));
        }
        S.set_submatrix(ns*channel_uses, 0, temp);
    }
    return S;
}

itpp::cmat STC::diag_pow(const itpp::cmat &in_mat, double in_exp)
//first input should be a diagonal square matrix with complex elements
{
    register int n;
    int dim = in_mat.rows();
    itpp::cmat out_mat(dim,dim);
    out_mat.zeros();
    for (n=0; n<dim; n++)
    {
        out_mat(n,n) = std::pow(in_mat(n,n), in_exp);
    }
    return out_mat;
}

itpp::mat STC::mat_pow(const itpp::mat &in_mat, int in_exp)
//square matrix power of integer exponent
{
    if (in_exp==0)
    {
        return itpp::eye(in_mat.rows());
    }
    itpp::mat out = in_mat;
    int abs_in_exp = std::abs(in_exp);
    register int n;
    for (n=1; n<abs_in_exp; n++)
    {
        out *= in_mat;
    }
    return (in_exp>0)?out:itpp::inv(out);
}

}
