/*!
 * \file
 * \brief Implementation of SISO modules for demappers
 * \author Bogdan Cristea
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

#include <itpp/comm/siso.h>
#include <limits>
#ifndef INFINITY
#define INFINITY std::numeric_limits<double>::infinity()
#endif

namespace itpp
{
void SISO::find_half_const(int &select_half, itpp::vec &re_part, itpp::bmat &re_bin_part, itpp::vec &im_part, itpp::bmat &im_bin_part)
/* finds real in imaginary parts of the constellation and its corresponding bits
 * this approach is used for equivalent channel according to Hassibi's model
 * the constellation must be quadratic and the number of bits per symbol must be a multiple of two
 */
{
    //values needed for initializations
    int const_size = itpp::pow2i(nb_bits_symb);//constellation size
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part
    //initialize output variables
    select_half = 0;
    re_part.set_size(half_len);
    re_bin_part.set_size(half_len, half_nb_bits_symb);
    re_part.zeros();
    re_part(0) = constellation(0).real();
    im_part.set_size(half_len);
    im_bin_part.set_size(half_len, half_nb_bits_symb);
    im_part.zeros();
    im_part(0) = constellation(0).imag();
    //select half for real (imaginary) to binary correspondence
    if (nb_bits_symb%2)
    {
        print_err_msg("SISO::find_half_const: number of bits per symbol must be a multiple of two");
        return;
    }
    const double min_diff = 1e-3;
    itpp::ivec idx = itpp::find(itpp::abs(itpp::real(constellation)-re_part(0))<min_diff);
    if (idx.length()!=half_len)
    {
        print_err_msg("SISO::find_half_const: the constellation must be quadratic");
        return;
    }
    itpp::bvec temp(nb_bits_symb);
    register int n;
    for (n=0; n<2; n++)
    {
        temp = bin_constellation.get_row(idx(n));
        re_bin_part.set_row(n,temp(0,half_nb_bits_symb-1));
    }
    select_half = (re_bin_part.get_row(0)==re_bin_part.get_row(1))?0:1;
    //algorithm
    double buffer;
    temp = bin_constellation.get_row(0);
    re_bin_part.set_row(0,temp(select_half*half_nb_bits_symb,(1+select_half)*half_nb_bits_symb-1));
    im_bin_part.set_row(0,temp((1-select_half)*half_nb_bits_symb,(2-select_half)*half_nb_bits_symb-1));
    int re_idx = 0;
    int im_idx = 0;
    for (n=1; n<const_size; n++)
    {
        temp = bin_constellation.get_row(n);
        buffer = constellation(n).real();
        if (fabs(itpp::prod(re_part-buffer))>min_diff)
        {
            re_idx++;
            re_part(re_idx) = buffer;
            re_bin_part.set_row(re_idx, temp(select_half*half_nb_bits_symb,(1+select_half)*half_nb_bits_symb-1));
        }
        buffer = constellation(n).imag();
        if (fabs(itpp::prod(im_part-buffer))>min_diff)
        {
            im_idx++;
            im_part(im_idx) = buffer;
            im_bin_part.set_row(im_idx, temp((1-select_half)*half_nb_bits_symb,(2-select_half)*half_nb_bits_symb-1));
        }
    }
}

void SISO::EquivRecSig(itpp::vec &x_eq, const itpp::cmat &rec_sig)
//finds equivalent received signal with real coefficients
//the equivalent received signal follows the model of Hassibi's paper
//ouput:
//x_eq - equivalent received signal with real coefficients
//inputs:
//rec_sig - received signal
{
    for (int k=0; k<nb_rec_ant; k++)
    {
        x_eq.set_subvector(k*2*block_duration, itpp::real(rec_sig.get_col(k)));
        x_eq.set_subvector(k*2*block_duration+block_duration, itpp::imag(rec_sig.get_col(k)));
    }
}

void SISO::EquivCh(itpp::mat &H_eq, const itpp::cvec &H)
//finds equivalent channel with real coefficients following the model of Hassibi's paper
//output:
//H_eq - equivalent channel
//input:
//H - channel matrix
{
    itpp::mat Aq(2*block_duration,2*nb_em_ant);
    itpp::mat Bq(2*block_duration,2*nb_em_ant);
    itpp::cmat temp(block_duration,nb_em_ant);
    itpp::vec h(2*nb_em_ant);
    itpp::mat AhBh(2*block_duration,2);
    register int n,k;
    for (k=0; k<symbols_block; k++)
    {
        temp = ST_gen1.get(k*block_duration,k*block_duration+block_duration-1,0,nb_em_ant-1);
        Aq.set_submatrix(0, 0, itpp::real(temp));
        Aq.set_submatrix(0, nb_em_ant, -itpp::imag(temp));
        Aq.set_submatrix(block_duration, 0, itpp::imag(temp));
        Aq.set_submatrix(block_duration, nb_em_ant, itpp::real(temp));
        temp = ST_gen2.get(k*block_duration,k*block_duration+block_duration-1,0,nb_em_ant-1);
        Bq.set_submatrix(0, 0, -itpp::imag(temp));
        Bq.set_submatrix(0, nb_em_ant, -itpp::real(temp));
        Bq.set_submatrix(block_duration, 0, itpp::real(temp));
        Bq.set_submatrix(block_duration, nb_em_ant, -itpp::imag(temp));
        for (n=0; n<nb_rec_ant; n++)
        {
            h.set_subvector(0, real(H.mid(n*nb_em_ant,nb_em_ant)));
            h.set_subvector(nb_em_ant, imag(H.mid(n*nb_em_ant,nb_em_ant)));
            AhBh.set_col(0, Aq*h);
            AhBh.set_col(1, Bq*h);
            H_eq.set_submatrix(2*block_duration*n, 2*k, AhBh);
        }
    }
}

void SISO::compute_symb_stats(itpp::vec &Es, itpp::vec &Vs,
		int ns, int select_half, const itpp::vec &apriori_data,
		const itpp::vec &re_part, const itpp::vec &im_part,
		const itpp::bmat &re_bin_part, const itpp::bmat &im_bin_part)
{
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part

    Es.zeros();
    Vs.zeros();
    for (int q = 0; q < symbols_block; q++)
    {
        int index = q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
        for (int k = 0; k < half_len; k++)
        {
            Es(2*q) += re_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                                               apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                               1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
            Es(1+2*q) += im_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                                               apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                               1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
            Vs(2*q) += itpp::sqr(re_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                           apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                           1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
            Vs(1+2*q) += itpp::sqr(im_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                           apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                           1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
        }
        Vs(2*q) -= itpp::sqr(Es(2*q));
        Vs(1+2*q) -= itpp::sqr(Es(1+2*q));
    }
}

void SISO::Hassibi_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//maxlogMAP algorithm for ST block codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks of ST matrices/int period
    double N0 = 2*sigma2;//noise DSP
    int nb_all_symb = itpp::pow2i(nb_bits_symb*symbols_block);//nb. of all possible input symbols as a binary vector
    double nom,denom;//nominator and denominator of extrinsic information
    itpp::bvec bin_frame(nb_bits_symb*symbols_block);//binary frame at channel input
    itpp::bmat mat_bin_frame(nb_bits_symb, symbols_block);
    itpp::vec symb_frame_eq(2*symbols_block);//frame of symbols at equivalent channel input
    double temp;
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);//equivalent channel matrix
    itpp::vec x_eq(2*block_duration*nb_rec_ant);//equivalent received signal
    register int ns,q,nb,n,k;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    //main loop
    for (ns=0; ns<nb_subblocks; ns++)//for each subblock
    {
        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        //compute the LLR of each bit in a frame of symbols_block symbols
        for (q=0; q<symbols_block; q++)//for each symbol in a subblock
        {
            for (nb=0; nb<nb_bits_symb; nb++)//for a given bit try all possible sollutions for the input symbol vector
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (n=0; n<nb_all_symb; n++)//all possible symbols
                {
                    bin_frame = itpp::dec2bin(nb_bits_symb*symbols_block, n);
                    mat_bin_frame = itpp::reshape(bin_frame, nb_bits_symb, symbols_block);
                    for (k=0; k<symbols_block; k++)
                    {
                        symb_frame_eq(2*k) = constellation(itpp::bin2dec(mat_bin_frame.get_col(k))).real();
                        symb_frame_eq(1+2*k) = constellation(itpp::bin2dec(mat_bin_frame.get_col(k))).imag();
                    }
                    temp = -itpp::sum_sqr(x_eq-H_eq*symb_frame_eq)/N0+\
                           itpp::to_vec(bin_frame)*apriori_data.mid(ns*nb_bits_symb*symbols_block,nb_bits_symb*symbols_block);
                    if (bin_frame(nb+q*nb_bits_symb))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                index = nb+q*nb_bits_symb+ns*nb_bits_symb*symbols_block;
                extrinsic_data(index) = (nom-denom)-apriori_data(index);
            }//bits/symbol
        }//symbols/subblock
    }//subblocks
}

void SISO::GA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
// Gaussian Approximation algorithm for ST codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

    //equivalent channel
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::vec Ey(2*block_duration*nb_rec_ant);
    itpp::mat Cy(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::mat Cy_inv(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::vec x_eq(2*block_duration*nb_rec_ant);
    itpp::vec EZeta(2*block_duration*nb_rec_ant);
    itpp::mat CZeta_inv(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    double nom,denom;
    double temp;
    register int ns,q,p,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0; ns<nb_subblocks; ns++)//subblock by subblock
    {
        //mean and variance of real and imaginary parts of emitted symbols
    	compute_symb_stats(Es, Vs, ns, select_half, apriori_data,
    			re_part, im_part, re_bin_part, im_bin_part);

        //find equivalent channel
        EquivCh(H_eq, c_impulse_response.get_col(ns));

        //compute E[y] and Cov[y]
        Ey.zeros();
        Cy = sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        for (q=0; q<symbols_block; q++)
        {
            //real & imaginary
            Ey += (H_eq.get_col(2*q)*Es(2*q)+H_eq.get_col(1+2*q)*Es(1+2*q));
            Cy += (itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Vs(2*q))+\
                   itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Vs(1+2*q)));
        }

        //inverse of Cov[y]
        Cy_inv = itpp::inv(Cy);

        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));

        //compute extrinsic information of coded bits
        for (q=0; q<symbols_block; q++)
        {
            //real part
            EZeta = Ey-H_eq.get_col(2*q)*Es(2*q);
            CZeta_inv = Cy_inv+itpp::outer_product(Cy_inv*\
                                                   ((Vs(2*q)/(1-(((H_eq.get_col(2*q)).transpose()*Cy_inv)*(H_eq.get_col(2*q)*Vs(2*q)))(0)))*\
                                                    H_eq.get_col(2*q)), Cy_inv.transpose()*H_eq.get_col(2*q));
            index = select_half*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0; p<half_nb_bits_symb; p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0; cs<half_len; cs++)
                {
                    temp = -0.5*((x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta).transpose()*CZeta_inv*(x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta))(0)+\
                           itpp::to_vec(re_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (re_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
            //imaginary part
            EZeta = Ey-H_eq.get_col(1+2*q)*Es(1+2*q);
            CZeta_inv = Cy_inv+itpp::outer_product(Cy_inv*\
                                                   ((Vs(1+2*q)/(1-(((H_eq.get_col(1+2*q)).transpose()*Cy_inv)*(H_eq.get_col(1+2*q)*Vs(1+2*q)))(0)))*\
                                                    H_eq.get_col(1+2*q)), Cy_inv.transpose()*H_eq.get_col(1+2*q));
            index = (1-select_half)*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0; p<half_nb_bits_symb; p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0; cs<half_len; cs++)
                {
                    temp = -0.5*((x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta).transpose()*CZeta_inv*(x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta))(0)+\
                           itpp::to_vec(im_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (im_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
        }
    }//subblock by subblock
}

void SISO::sGA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//simplified Gaussian Approximation algorithm for ST codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = (int)(rec_sig.rows()/block_duration);//number of subblocks
    int half_nb_bits_symb = (int)(nb_bits_symb/2);
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

    //equivalent channel
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);

    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::vec Ey(2*block_duration*nb_rec_ant);
    itpp::mat Cy(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::vec x_eq(2*block_duration*nb_rec_ant);
    itpp::vec EZeta(2*block_duration*nb_rec_ant);
    itpp::vec CZeta(2*block_duration*nb_rec_ant);
    double nom,denom;
    double temp;
    register int ns,q,p,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0; ns<nb_subblocks; ns++)//subblock by subblock
    {
        //mean and variance of real and imaginary parts of emitted symbols
    	compute_symb_stats(Es, Vs, ns, select_half, apriori_data,
    			re_part, im_part, re_bin_part, im_bin_part);

        //find equivalent channel
        EquivCh(H_eq, c_impulse_response.get_col(ns));

        //compute E[y] and Cov[y]
        Ey.zeros();
        Cy = sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        for (q=0; q<symbols_block; q++)
        {
            //real & imaginary
            Ey += (H_eq.get_col(2*q)*Es(2*q)+H_eq.get_col(1+2*q)*Es(1+2*q));
            Cy += (itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Vs(2*q))+\
                   itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Vs(1+2*q)));
        }

        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));

        //compute extrinsic INFINITYormation of coded bits
        for (q=0; q<symbols_block; q++)
        {
            //real part
            EZeta = Ey-H_eq.get_col(2*q)*Es(2*q);
            CZeta = diag(Cy-itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Vs(2*q)));
            index = select_half*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0; p<half_nb_bits_symb; p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0; cs<half_len; cs++)
                {
                    temp = -0.5*itpp::sum(itpp::elem_div(sqr(x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta), CZeta))+\
                           itpp::to_vec(re_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (re_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
            //imaginary part
            EZeta = Ey-H_eq.get_col(1+2*q)*Es(1+2*q);
            CZeta = itpp::diag(Cy-itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Vs(1+2*q)));
            index = (1-select_half)*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0; p<half_nb_bits_symb; p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0; cs<half_len; cs++)
                {
                    temp = -0.5*itpp::sum(itpp::elem_div(sqr(x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta), CZeta))+\
                           itpp::to_vec(im_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (im_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
        }
    }//subblock by subblock
}

void SISO::mmsePIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//MMSE Parallel Interference Canceller
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_const_len = itpp::pow2i(half_nb_bits_symb);
    int nb_bits_subblock = nb_bits_symb*symbols_block;//number of coded bits in an ST block
    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::mat K(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::mat K_inv(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::vec x_eq(2*nb_rec_ant*block_duration);
    itpp::vec interf(2*symbols_block);
    itpp::vec temp(2*nb_rec_ant*block_duration);
    itpp::vec w(2*nb_rec_ant*block_duration);//filter impulse response
    double s_tilde;
    double mu_res;
    double sigma2_res;
    double nom,denom;
    double tmp;
    register int ns,q,k,s;
    int index;

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);
    double part_var = 1/(double)(2*nb_em_ant);//real and imaginary part variance

    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0; ns<nb_subblocks; ns++)//compute block by block
    {
        //mean and variance of real and imaginary parts of emitted symbols
    	compute_symb_stats(Es, Vs, ns, select_half, apriori_data,
    			re_part, im_part, re_bin_part, im_bin_part);

        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //compute invariant inverse
        K = H_eq*diag(Vs)*H_eq.transpose()+sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        K_inv = itpp::inv(K);
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        for (q=0; q<symbols_block; q++)//symbols/block
        {
            //compute the extrinsic information of coded bits
            //real part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q);
            w = (part_var*temp.transpose())*(K_inv-itpp::outer_product(K_inv*(((part_var-Vs(2*q))/ \
                                             (1+((temp.transpose()*K_inv)*(temp*(part_var-Vs(2*q))))(0)))*temp), K_inv.transpose()*temp));
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = select_half*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0; k<half_nb_bits_symb; k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0; s<half_const_len; s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(re_part(s))-Vs(2*q)))))*w)(0)-itpp::sqr(re_part(s)*mu_res);//variance of the filtered signal
                    tmp = -itpp::sqr(s_tilde-mu_res*re_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(re_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (re_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end real part
            //imaginary part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q+1) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q+1);
            w = (part_var*temp.transpose())*(K_inv-itpp::outer_product(K_inv*(((part_var-Vs(1+2*q))/ \
                                             (1+((temp.transpose()*K_inv)*(temp*(part_var-Vs(1+2*q))))(0)))*temp), K_inv.transpose()*temp));
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = (1-select_half)*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0; k<half_nb_bits_symb; k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0; s<half_const_len; s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(im_part(s))-Vs(1+2*q)))))*w)(0)-itpp::sqr(im_part(s)*mu_res);//variance of the filtered signal
                    tmp = -itpp::sqr(s_tilde-mu_res*im_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(im_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (im_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end imaginary part
        }//symbols/block
    }//block by block
}

void SISO::zfPIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//ZF Parallel Interference Canceller
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_const_len = itpp::pow2i(half_nb_bits_symb);
    int nb_bits_subblock = nb_bits_symb*symbols_block;//number of coded bits in an ST block
    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::mat K(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::vec x_eq(2*nb_rec_ant*block_duration);
    itpp::vec interf(2*symbols_block);
    itpp::vec temp(2*nb_rec_ant*block_duration);
    itpp::vec w(2*nb_rec_ant*block_duration);//filter impulse response
    double s_tilde;
    double mu_res;
    double sigma2_res;
    double nom,denom;
    double tmp;
    register int ns,q,k,s;
    int index;

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0; ns<nb_subblocks; ns++)//compute block by block
    {
        //mean and variance of real and imaginary parts of emitted symbols
    	compute_symb_stats(Es, Vs, ns, select_half, apriori_data,
    			re_part, im_part, re_bin_part, im_bin_part);

        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //compute invariant inverse
        K = H_eq*itpp::diag(Vs)*H_eq.transpose()+sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        for (q=0; q<symbols_block; q++)//symbols/block
        {
            //compute the extrinsic information of coded bits
            //real part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q);
            w = temp/(temp*temp);//filter impulse response
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = select_half*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0; k<half_nb_bits_symb; k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0; s<half_const_len; s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(re_part(s))-Vs(2*q)))))*w)(0)-itpp::sqr(re_part(s)*mu_res);
                    tmp = -itpp::sqr(s_tilde-mu_res*re_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(re_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (re_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end real part
            //imaginary part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q+1) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q+1);
            w = temp/(temp*temp);//filter impulse response
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = (1-select_half)*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0; k<half_nb_bits_symb; k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0; s<half_const_len; s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(im_part(s))-Vs(1+2*q)))))*w)(0)-itpp::sqr(im_part(s)*mu_res);
                    tmp = -itpp::sqr(s_tilde-mu_res*im_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(im_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (im_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end imaginary part
        }//symbols/block
    }//block by block
}

void SISO::Alamouti_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//maxlogMAP algorithm for Alamouti ST code
{
    //matched filter
    int int_len = apriori_data.length();//interleaver length
    int nb_symb = (int)(int_len/nb_bits_symb);//number of symbols/block
    itpp::cvec comb_sig(nb_symb);
    comb_sig.zeros();
    itpp::cmat conj_H = itpp::conj(c_impulse_response);
    itpp::cmat conj_X = itpp::conj(rec_sig);
    register int nr,n,cs;
    for (nr=0; nr<nb_rec_ant; nr++)
    {
        for (n=0; n<(nb_symb/2); n++)
        {
            comb_sig(2*n) += (conj_H(2*nr,n)*rec_sig(2*n,nr)+c_impulse_response(1+2*nr,n)*conj_X(1+2*n,nr));
            comb_sig(1+2*n) += (conj_H(1+2*nr,n)*rec_sig(2*n,nr)-c_impulse_response(2*nr,n)*conj_X(1+2*n,nr));
        }
    }

    //extrinsic information of coded bits
    int const_size = itpp::pow2i(nb_bits_symb);//constellation size
    double buffer;
    double nom,denom;
    double temp;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for (n=0; n<nb_symb; n++)
    {
        buffer = itpp::sum_sqr(itpp::abs(c_impulse_response.get_col(n/2)));
        for (nr=0; nr<nb_bits_symb; nr++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for (cs=0; cs<const_size; cs++)
            {
                temp = -itpp::sqr(comb_sig(n)-buffer*constellation(cs))/(2*buffer*sigma2)+ \
                       itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(n*nb_bits_symb, nb_bits_symb);
                if (bin_constellation(cs,nr))
                    nom = std::max(nom, temp);
                else
                    denom = std::max(denom, temp);
            }
            index = n*nb_bits_symb+nr;
            extrinsic_data(index) = (nom-denom)-apriori_data(index);//extrinsic information
        }
    }
}

void SISO::demodulator_logMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data)
/// logMAP demodulator
{
    int nb_symb = rec_sig.length();
    int const_size = itpp::pow2i(nb_bits_symb);
    double nom,denom,temp;
    register int k,i,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for (k=0; k<nb_symb; k++)
    {
        for (i=0; i<nb_bits_symb; i++)
        {
            nom = 0;
            denom = 0;
            for (cs=0; cs<const_size; cs++)
            {
                temp = -itpp::sqr(rec_sig(k)-c_impulse_response(0,k)*constellation(cs))/(2*sigma2)+\
                       itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(k*nb_bits_symb, nb_bits_symb);
                if (bin_constellation(cs,i))
                    nom += std::exp(temp);
                else
                    denom += std::exp(temp);
            }
            index = k*nb_bits_symb+i;
            extrinsic_data(index) = std::log(nom/denom)-apriori_data(index);//extrinsic information
        }
    }
}

void SISO::demodulator_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data)
/// maxlogMAP demodulator
{
    int nb_symb = rec_sig.length();
    int const_size = itpp::pow2i(nb_bits_symb);
    double nom,denom,temp;
    register int k,i,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for (k=0; k<nb_symb; k++)
    {
        for (i=0; i<nb_bits_symb; i++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for (cs=0; cs<const_size; cs++)
            {
                temp = -itpp::sqr(rec_sig(k)-c_impulse_response(0,k)*constellation(cs))/(2*sigma2)+\
                       itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(k*nb_bits_symb, nb_bits_symb);
                if (bin_constellation(cs,i))
                    nom = std::max(nom, temp);
                else
                    denom = std::max(denom, temp);
            }
            index = k*nb_bits_symb+i;
            extrinsic_data(index) = (nom-denom)-apriori_data(index);//extrinsic information
        }
    }
}
}//end namespace tr
