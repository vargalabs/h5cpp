/*!
 * \file
 * \brief Implementation of SISO modules for descrambler and MUDs
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
void SISO::descrambler(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
  inputs:
  intrinsic_coded - intrinsic information of coded bits (repetition code output)
  apriori_data - a priori information of informational bits (repetition code input)
  outputs:
  extrinsic_coded - extrinsic information of coded bits
  extrinsic_data - Logarithm of Likelihood Ratio of informational bits
*/
{
    //get parameters
    int nb_bits = apriori_data.length();
    int Nc = scrambler_pattern.length();
    //implementation
    extrinsic_data.set_size(nb_bits);
    extrinsic_coded.set_size(nb_bits*Nc);
    register int n,k;
#pragma omp parallel for private(n,k)
    for (k=0; k<nb_bits; k++)
    {
        extrinsic_data[k] = 0;//apriori_data[k];//add a priori info
        for (n=0; n<Nc; n++)
        {
            extrinsic_data[k] += (1-2*double(scrambler_pattern[n]))*intrinsic_coded[n+k*Nc];
        }
        for (n=0; n<Nc; n++)
        {
            extrinsic_coded[n+k*Nc] = (1-2*double(scrambler_pattern[n]))*extrinsic_data[k]-intrinsic_coded[n+k*Nc];
        }
    }
}

void SISO::zpFIRfilter(itpp::vec &filt, const itpp::vec &h, const itpp::vec &sig)
//FIR filter for a zero padded signal (L zeros are added at the end of the signal)
{
    //get parameters
    int L = h.length()-1;
    int N = sig.length();
    //implementation
    register int n,l;
#pragma omp parallel for private(n,l)
    for (n=0; n<(N+L); n++)
    {
        filt[n] = 0;
        for (l=0; l<=L; l++)
        {
            if ((n-l)<0)
            {
                break;//channel has state 0 at the beginning
            }
            if ((n-l)>=N)
            {
                continue;//channel has state 0 in the end
            }
            filt[n] += (h[l]*sig[n-l]);
        }
    }
}

void SISO::gen_hyperTrellis(void)
/* generates channel hyper trellis for binary symbols
 * the channel is a MISO system
 * BPSK mapping: 0->+1, 1->-1
 */
{
    //get parameters
    int nb_usr = impulse_response.rows();
    int ch_order = impulse_response.cols()-1;
    int p_order = prec_gen.length()-1;
    int max_order = std::max(ch_order, p_order);

    //initialize hypertrellis
    chtrellis.numInputSymbols = itpp::pow2i(nb_usr);
    int mem_len = nb_usr*max_order;
    if (mem_len>=(int)(8*sizeof(int)))
    {
        std::string msg = "SISO::gen_hyperTrellis: memory length of the hyperchannel should be at most ";
        msg += itpp::to_str(8*sizeof(int)-1);
        msg += ". Try to lower the number of users, channel order or the order of the precoding polynomial (if any)";
        print_err_msg(msg);
        return;
    }
    chtrellis.stateNb = itpp::pow2i(mem_len);
    try
    {
        unsigned int len =  static_cast<unsigned int>(chtrellis.stateNb)*static_cast<unsigned int>(chtrellis.numInputSymbols);
        chtrellis.nextState = new int[len];
        chtrellis.prevState = new int[len];
        chtrellis.output = new double[len];
        chtrellis.input = new int[len];
    } catch (std::bad_alloc)
    {
        std::string msg = "SISO::gen_hyperTrellis: not enough memory for the channel trellis variables. The number of states is ";
        msg += itpp::to_str(chtrellis.stateNb);
        msg += " and the number of input symbols ";
        msg += itpp::to_str(chtrellis.numInputSymbols);
        print_err_msg(msg);
        return;
    }
    itpp::ivec index(chtrellis.stateNb);
    index.zeros();
    itpp::bvec hyper_ch_mem(mem_len);
    itpp::bvec hyper_ch_in(nb_usr);
    itpp::bvec hyper_states(mem_len);
    itpp::bin feedback;

    //create hypertrellis
    register int n,k,p,r;
    int buffer;
    double hyper_ch_out;
    for (k=0; k<chtrellis.stateNb; k++)
    {
        hyper_ch_mem = itpp::dec2bin(mem_len, k);//initial state
        for (n=0; n<chtrellis.numInputSymbols; n++)
        {
            hyper_ch_in = itpp::dec2bin(nb_usr, n);//MISO channel input
            hyper_ch_out = 0;
            for (r=0; r<nb_usr; r++)
            {
                buffer = r*max_order;
                //precoder
                feedback = hyper_ch_in[r];
                for (p=1; p<=p_order; p++)
                {
                    if (prec_gen(p))
                    {
                        feedback ^= hyper_ch_mem[p-1+buffer];//xor
                    }
                }
                //FIR channel output
                hyper_ch_out += (1-2*double(feedback))*impulse_response(r,0);
                for (p=0; p<ch_order; p++)
                {
                    hyper_ch_out += (1-2*double(hyper_ch_mem[p+buffer]))*impulse_response(r,p+1);//hyper channel output for user r
                }
                //(equivalent) channel next state
                hyper_states[buffer] = feedback;
                for (p=0; p<(max_order-1); p++)
                {
                    hyper_states[p+buffer+1] = hyper_ch_mem[p+buffer];//next hyper state for user r
                }
            }
            chtrellis.output[k+n*chtrellis.stateNb] = hyper_ch_out;
            buffer = itpp::bin2dec(hyper_states);//next state from an initial state and a given input
            chtrellis.nextState[k+n*chtrellis.stateNb] = buffer;
            chtrellis.prevState[buffer+index[buffer]*chtrellis.stateNb] = k;
            chtrellis.input[buffer+index[buffer]*chtrellis.stateNb] = n;
            index[buffer]++;
        }
    }
}

/// Maximum A Posteriori algorithm for Multi-User Detection in IDMA systems
/** uses max log MAP algorithm
 * use with care for large number of users and/or FIR channel order
 */
void SISO::mud_maxlogMAP(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* output:
 * extrinsic_data - extrinsic information for the chips (usr_nb x block_len)
 * inputs:
 * rec_sig - received signal (1 x block_len)
 * apriori_data - a priori information for the chips (usr_nb x block_len)
 */
{
    //get parameters
    int nb_usr = apriori_data.rows();
    int block_len = apriori_data.cols();

    //init trellis
    gen_hyperTrellis();

    //initial conditions for A = log(alpha) and B = log(beta)
    double *A = NULL,*B = NULL;
    try
    {
        A = new double[chtrellis.stateNb*(block_len+1)];
        B = new double[chtrellis.stateNb*(block_len+1)];
    } catch (std::bad_alloc)
    {
        std::string msg = "SISO::mud_maxlogMAP: Not enough memory for alphas and betas. The number of states is ";
        msg += itpp::to_str(chtrellis.stateNb);
        msg += " and the block length ";
        msg += itpp::to_str(block_len);
        print_err_msg(msg);
    }
    register int n;
    A[0] = 0;
    B[block_len*chtrellis.stateNb] = 0;
    double buffer = (tail?-INFINITY:0);
#pragma omp parallel for private(n)
    for (n=1; n<chtrellis.stateNb; n++)
    {
        A[n] = -INFINITY;
        B[n+block_len*chtrellis.stateNb] = buffer;//if tail==false the final state is not known
    }

    //compute log(alpha) (forward recursion)
    register int s,k;
    int sp,i;
    itpp::bvec in_chips(nb_usr);
#pragma omp parallel sections private(n,buffer,s,k,sp,in_chips)
    {
        for (n=1; n<=block_len; n++)
        {
            buffer = -INFINITY;//normalization factor
            for (s=0; s<chtrellis.stateNb; s++)
            {
                A[s+n*chtrellis.stateNb] = -INFINITY;
                for (k=0; k<chtrellis.numInputSymbols; k++)
                {
                    sp = chtrellis.prevState[s+k*chtrellis.stateNb];
                    i = chtrellis.input[s+k*chtrellis.stateNb];
                    in_chips = itpp::dec2bin(nb_usr, i);
                    A[s+n*chtrellis.stateNb] = std::max(A[s+n*chtrellis.stateNb], \
                                                        A[sp+(n-1)*chtrellis.stateNb]-itpp::sqr(rec_sig[n-1]-chtrellis.output[sp+i*chtrellis.stateNb])/(2*sigma2)+\
                                                        itpp::to_vec(in_chips)*apriori_data.get_col(n-1));
                }
                buffer = std::max(buffer, A[s+n*chtrellis.stateNb]);
            }
            //normalization
            for (s=0; s<chtrellis.stateNb; s++)
            {
                A[s+n*chtrellis.stateNb] -= buffer;
            }
        }

        //compute log(beta) (backward recursion)
#pragma omp section
        for (n=block_len-1; n>=0; n--)
        {
            buffer = -INFINITY;//normalization factor
            for (s=0; s<chtrellis.stateNb; s++)
            {
                B[s+n*chtrellis.stateNb] = -INFINITY;
                for (k=0; k<chtrellis.numInputSymbols; k++)
                {
                    sp = chtrellis.nextState[s+k*chtrellis.stateNb];
                    in_chips = itpp::dec2bin(nb_usr, k);
                    B[s+n*chtrellis.stateNb] = std::max(B[s+n*chtrellis.stateNb], \
                                                        B[sp+(n+1)*chtrellis.stateNb]-itpp::sqr(rec_sig[n]-chtrellis.output[s+k*chtrellis.stateNb])/(2*sigma2)+\
                                                        itpp::to_vec(in_chips)*apriori_data.get_col(n));
                }
                buffer = std::max(buffer, B[s+n*chtrellis.stateNb]);
            }
            //normalization
            for (s=0; s<chtrellis.stateNb; s++)
            {
                B[s+n*chtrellis.stateNb] -= buffer;
            }
        }
    }

    //compute extrinsic information
    double nom, denom;
    extrinsic_data.set_size(nb_usr,block_len);
    register int u;
#pragma omp parallel for private(u,n,s,k,nom,denom,in_chips,buffer)
    for (u=0; u<nb_usr; u++)
    {
        for (n=1; n<=block_len; n++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for (s=0; s<chtrellis.stateNb; s++)
            {
                for (k=0; k<chtrellis.numInputSymbols; k++)
                {
                    in_chips = itpp::dec2bin(nb_usr, k);
                    buffer = A[s+(n-1)*chtrellis.stateNb]+B[chtrellis.nextState[s+k*chtrellis.stateNb]+n*chtrellis.stateNb]-\
                             itpp::sqr(rec_sig[n-1]-chtrellis.output[s+k*chtrellis.stateNb])/(2*sigma2)+\
                             itpp::to_vec(in_chips)*apriori_data.get_col(n-1);
                    if (in_chips[u])
                    {
                        nom = std::max(nom, buffer);
                    }
                    else
                    {
                        denom = std::max(denom, buffer);
                    }
                }
            }
            extrinsic_data(u,n-1) = (nom-denom)-apriori_data(u,n-1);
        }
    }
    //free memory
    delete[] chtrellis.nextState;
    delete[] chtrellis.prevState;
    delete[] chtrellis.output;
    delete[] chtrellis.input;
    delete[] A;
    delete[] B;
}

/// Gaussian Chip Detector for IDMA systems
/** Use with care for large size of interleavers.
 */
void SISO::GCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* Gaussian Chip Detector
 * output:
 * extrinsic_data - extrinsic information of emitted chips
 * inputs:
 * rec_sig - received signal
 * apriori_data - a priori information of emitted chips
 */
{
    //get parameters
    int N = apriori_data.cols();//emitted frames of non-zero samples
    int K = apriori_data.rows();//number of users
    int L = impulse_response.cols()-1;//channel order
    //other parameters
    register int n,k;

    //mean and variance of each chip (NxK)
    itpp::mat Ex = -itpp::tanh(apriori_data/2.0);//take into account BPSK mapping
    itpp::mat Vx = 1.0-sqr(Ex);

    //expectation and variance of the received signal
    itpp::vec Er(N+L);
    Er.zeros();
    itpp::mat Covr;
    try
    {
        Covr.set_size(N+L,N+L);
    } catch (std::bad_alloc)
    {
        std::string msg = "SISO::GCD: not enough memory when allocating space for the covariance matrix. The interleaver length is ";
        msg += itpp::to_str(N);
        msg += ". Use sGCD instead or try to reduce the interleaver length";
        print_err_msg(msg);
        return;
    }
    //create eye function
    Covr.zeros();
    for (n=0; n<(N+L); n++)
        Covr(n,n) = 1;
    itpp::vec col(N+L);
    col.zeros();
    itpp::vec row(N);
    row.zeros();
    itpp::mat h_eq(N+L,N);
    for (k=0; k<K; k++)
    {
        col.replace_mid(0, impulse_response.get_row(k));
        row(0) = impulse_response(k,0);
        h_eq = itpp::toeplitz(col, row);
        Er += h_eq*Ex.get_row(k);
        Covr += (h_eq*itpp::diag(Vx.get_row(k)))*h_eq.T();
    }

    //extrinsic information
    double nom,denom;
    Er = rec_sig-Er;
    itpp::mat inv_Covr(N+L,N+L);
    inv_Covr = itpp::inv(Covr);
    itpp::mat inv_cov_zeta(N+L,N+L);
    itpp::vec rec_chip(N+L);
    extrinsic_data.set_size(K,N);
    for (k=0; k<K; k++)
    {
#pragma omp parallel for private(n,inv_cov_zeta,rec_chip,nom,denom) firstprivate(col)
        for (n=0; n<N; n++)
        {
            col.replace_mid(n, impulse_response.get_row(k));
            inv_cov_zeta = inv_Covr+itpp::outer_product(inv_Covr*col, inv_Covr.T()*(col*Vx(k,0)))/(1-(col*Vx(k,0))*(inv_Covr*col));
            rec_chip = Er-col*(+1-Ex(k,n));
            nom = -0.5*rec_chip*(inv_cov_zeta*rec_chip);
            rec_chip = Er-col*(-1-Ex(k,n));
            denom = -0.5*rec_chip*(inv_cov_zeta*rec_chip);
            extrinsic_data(k,n) = denom-nom;//take into account BPSK mapping
            col(n) = 0;
        }
    }
}

/// simplified Gaussian Chip Detector for IDMA systems
/** This algorithm is simplified and uses much less memory than its counterpart, the GCD.
 * Recommended to use in most cases.
 */
void SISO::sGCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* simplified GCD
 * output:
 * extrinsic_data - extrinsic information of emitted chips
 * inputs:
 * rec_sig - received signal
 * apriori_data - a priori information of emitted chips
 */
{
    //get parameters
    int N = apriori_data.cols();//emitted frames of non-zero samples
    int K = apriori_data.rows();//number of users
    int L = impulse_response.cols()-1;//channel order
    //other parameters
    register int n,k;

    //mean and variance of each chip (NxK)
    itpp::mat Ex = -itpp::tanh(apriori_data/2.0);//take into account BPSK mapping
    itpp::mat Vx = 1.0-itpp::sqr(Ex);

    //mean and variance for the samples of the received signal
    itpp::vec Er(N+L);
    Er.zeros();
    itpp::vec Vr = sigma2*itpp::ones(N+L);
    itpp::vec buffer(N+L);
    for (k=0; k<K; k++)
    {
        zpFIRfilter(buffer, impulse_response.get_row(k), Ex.get_row(k));
        Er += buffer;
        zpFIRfilter(buffer, itpp::sqr(impulse_response.get_row(k)), Vx.get_row(k));
        Vr += buffer;
    }

    //extrinsic information for the samples of the received signal
    Er = rec_sig-Er;
    itpp::vec ch(L+1);
    extrinsic_data.set_size(K,N);
    for (k=0; k<K; k++)
    {
        ch = impulse_response.get_row(k);
#pragma omp parallel for private(n)
        for (n=0; n<N; n++)
        {
            extrinsic_data(k,n) = -2*itpp::elem_div(ch, Vr.mid(n,L+1)-sqr(ch)*Vx(k,n))*(Er.mid(n,L+1)+ch*Ex(k,n));//take into account BPSK mapping
        }
    }
}
}//end namespace tr
