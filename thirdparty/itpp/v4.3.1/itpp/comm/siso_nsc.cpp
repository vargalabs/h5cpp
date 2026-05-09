/*!
 * \file
 * \brief Implementation of SISO modules for NSC codes
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
#include <itpp/base/itcompat.h>
#include <limits>
#ifndef INFINITY
#define INFINITY std::numeric_limits<double>::infinity()
#endif

namespace itpp
{
void SISO::gen_nsctrellis(void)
/*
 generate trellis for non systematic non recursive convolutional codes
 r - number of outputs of the CC
 mem_len - memory length of the CC
 */
{
    //get parameters
    int r = gen.rows();
    int mem_len = gen.cols()-1;
    //other parameters
    register int n,k,j,p;
    itpp::bin inputs[] = {0,1};
    int index;

    nsctrellis.stateNb = (1<<mem_len);
    nsctrellis.output = new double[nsctrellis.stateNb*2*r];
    nsctrellis.nextState = new int[nsctrellis.stateNb*2];
    nsctrellis.prevState = new int[nsctrellis.stateNb*2];
    nsctrellis.input = new int[nsctrellis.stateNb*2];

    itpp::bvec enc_mem(mem_len);
    itpp::bin out;
#pragma omp parallel for private(n,k,j,p,enc_mem,out)
    for (n=0; n<nsctrellis.stateNb; n++) //initial state
    {
        enc_mem = itpp::dec2bin(mem_len, n);
        //output
        for (k=0; k<2; k++)
        {
            for (j=0; j<r; j++)
            {
                out = inputs[k]*gen(j,0);
                for (p=1; p<=mem_len; p++)
                {
                    out ^= (enc_mem[p-1]*gen(j,p));//0 or 1
                }
                nsctrellis.output[n+k*nsctrellis.stateNb+j*nsctrellis.stateNb*2] = double(out);
            }
        }
        //next state
        for (j=(mem_len-1); j>0; j--)
        {
            enc_mem[j] = enc_mem[j-1];
        }
        for (k=0; k<2; k++)
        {
            enc_mem[0] = inputs[k];
            nsctrellis.nextState[n+k*nsctrellis.stateNb] = itpp::bin2dec(enc_mem, true);
        }
    }

#pragma omp parallel for private(n,k,j,index)
    for (j=0; j<nsctrellis.stateNb; j++)
    {
        index = 0;
        for (n=0; n<nsctrellis.stateNb; n++)
        {
            for (k=0; k<2; k++)
            {
                if (nsctrellis.nextState[n+k*nsctrellis.stateNb]==j)
                {
                    nsctrellis.prevState[j+index*nsctrellis.stateNb] = n;
                    nsctrellis.input[j+index*nsctrellis.stateNb]  = k;//0 or 1
                    index++;
                }
            }
        }
    }
}

void SISO::nsc_logMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 * generalized decoder for NSC codes (after the NSC code a scrambler of pattern phi is used) using log MAP alg.
 * extrinsic_coded - extrinsic information of coded bits (output if the shaping filter)
 * extrinsic_data - extrinsic information of data bits
 * intrinsic_coded - intrinsic information of coded bits
 * apriori_data - a priori information of data bits
 */
{
    //get parameters
    int N = apriori_data.length();
    int Nc = scrambler_pattern.length();
    int r = gen.rows();//number of outputs of the CC
    //other parameters
    register int n,k,m,mp,j,i;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    int index;

    //initialize trellis
    gen_nsctrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[nsctrellis.stateNb*(N+1)];
    double* B = new double[nsctrellis.stateNb*(N+1)];
    A[0] = 0;
    B[N*nsctrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
#pragma omp parallel for private(n)
    for (n=1; n<nsctrellis.stateNb; n++)
    {
        A[n] = -INFINITY;
        B[n+N*nsctrellis.stateNb] = sum;//if tail==false the final state is not known
    }

#pragma omp parallel sections private(n,sum,m,k,i,j,C)
    {
        //forward recursion
        for (n=1; n<(N+1); n++)
        {
            sum = 0;//normalization factor
            for (m=0; m<nsctrellis.stateNb; m++) //final state
            {
                for (k=0; k<2; k++)
                {
                    pstates[k] = nsctrellis.prevState[m+nsctrellis.stateNb*k];//determine previous states
                    inputs[k] = nsctrellis.input[m+nsctrellis.stateNb*k];//determine input
                    C[k] = (inputs[k]?apriori_data[n-1]:0);//compute log of gamma
                }
                for (i=0; i<2; i++)//for each C[i]
                {
                    for (k=0; k<r; k++)
                    {
                        for (j=0; j<Nc; j++)
                        {
                            C[i] += nsctrellis.output[pstates[i]+inputs[i]*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+(n-1)*Nc*r];
                        }
                    }
                }
                A[m+n*nsctrellis.stateNb] = itpp::log_add(A[pstates[0]+(n-1)*nsctrellis.stateNb]+C[0], A[pstates[1]+(n-1)*nsctrellis.stateNb]+C[1]);
                sum += std::exp(A[m+n*nsctrellis.stateNb]);
            }
            //normalization
            sum = std::log(sum);
            if (std::isinf(sum))
            {
                continue;
            }
            for (m=0; m<nsctrellis.stateNb; m++)
            {
                A[m+n*nsctrellis.stateNb] -= sum;
            }
        }

        //backward recursion
#pragma omp section
        for (n=N-1; n>=0; n--)
        {
            sum = 0;//normalisation factor
            for (m=0; m<nsctrellis.stateNb; m++) //initial state
            {
                for (k=0; k<2; k++)
                {
                    nstates[k] = nsctrellis.nextState[m+k*nsctrellis.stateNb];//determine next states
                    C[k] = (k?apriori_data[n]:0);//compute log of gamma
                }
                for (i=0; i<2; i++)
                {
                    for (k=0; k<r; k++)
                    {
                        for (j=0; j<Nc; j++)
                        {
                            C[i] += nsctrellis.output[m+i*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                        }
                    }
                }
                B[m+n*nsctrellis.stateNb] = itpp::log_add(B[nstates[0]+(n+1)*nsctrellis.stateNb]+C[0], B[nstates[1]+(n+1)*nsctrellis.stateNb]+C[1]);
                sum += std::exp(B[m+n*nsctrellis.stateNb]);
            }
            //normalisation
            sum = std::log(sum);
            if (std::isinf(sum))
            {
                continue;
            }
            for (m=0; m<nsctrellis.stateNb; m++)
            {
                B[m+n*nsctrellis.stateNb] -= sum;
            }
        }
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
#pragma omp parallel for private(n,k,sum,sumbis)
    for (n=1; n<(N+1); n++)
    {
        sum = 0;
        sumbis = 0;
        for (k=0; k<(nsctrellis.stateNb/2); k++)
        {
            sum += std::exp(A[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]+B[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]);//nominator
            sumbis += std::exp(A[k+n*nsctrellis.stateNb]+B[k+n*nsctrellis.stateNb]);//denominator
        }
        extrinsic_data[n-1] = std::log(sum/sumbis)-apriori_data[n-1];
    }

    //compute extrinsic_coded
    double *sum0 = new double[r];
    double *sum1 = new double[r];
    extrinsic_coded.set_size(N*Nc*r);
    for (n=0; n<N; n++)
    {
        for (k=0; k<r; k++)
        {
            sum0[k] = 0;
            sum1[k] = 0;
        }
        for (mp=0; mp<nsctrellis.stateNb; mp++) //initial state
        {
            for (i=0; i<2; i++)
            {
                m = nsctrellis.nextState[mp+i*nsctrellis.stateNb];//final state
                //compute log of sigma
                index = (m>=(nsctrellis.stateNb/2));//0 if input is 0, 1 if input is 1
                C[0] = (index?apriori_data[n]:0);
                for (k=0; k<r; k++)
                {
                    for (j=0; j<Nc; j++)
                    {
                        C[0] += nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                    }
                }
                C[1] = A[mp+n*nsctrellis.stateNb]+C[0]+B[m+(n+1)*nsctrellis.stateNb];//this is only a buffer
                //compute sums
                for (k=0; k<r; k++)
                {
                    if (nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2])
                    {
                        sum1[k] += std::exp(C[1]);
                    }
                    else
                    {
                        sum0[k] += std::exp(C[1]);
                    }
                }
            }
        }
        for (k=0; k<r; k++)
        {
            for (j=0; j<Nc; j++)
            {
                index = j+k*Nc+n*r*Nc;
                extrinsic_coded[index] = (1-2*double(scrambler_pattern[j]))*std::log(sum1[k]/sum0[k])-intrinsic_coded[index];
            }
        }
    }

    //free memory
    delete[] nsctrellis.output;
    delete[] nsctrellis.nextState;
    delete[] nsctrellis.prevState;
    delete[] nsctrellis.input;
    delete[] A;
    delete[] B;
    delete[] sum0;
    delete[] sum1;
}

void SISO::nsc_maxlogMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 * generalized decoder for NSC codes (after the NSC code a scrambler of pattern phi is used) using max log MAP alg.
 * extrinsic_coded - extrinsic information of coded bits (output if the shaping filter)
 * extrinsic_data - extrinsic information of data bits
 * intrinsic_coded - intrinsic information of coded bits
 * apriori_data - a priori information of data bits
 */
{
    //get parameters
    int N = apriori_data.length();
    int Nc = scrambler_pattern.length();
    int r = gen.rows();//number of outputs of the CC
    //other parameters
    register int n,k,m,mp,j,i;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    int index;

    //initialize trellis
    gen_nsctrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[nsctrellis.stateNb*(N+1)];
    double* B = new double[nsctrellis.stateNb*(N+1)];
    A[0] = 0;
    B[N*nsctrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
#pragma omp parallel for private(n)
    for (n=1; n<nsctrellis.stateNb; n++)
    {
        A[n] = -INFINITY;
        B[n+N*nsctrellis.stateNb] = sum;//if tail==false the final state is not known
    }

    //forward recursion
#pragma omp parallel sections private(n,sum,m,k,i,j,C)
    {
        for (n=1; n<(N+1); n++)
        {
            sum = -INFINITY;//normalisation factor
            for (m=0; m<nsctrellis.stateNb; m++) //final state
            {
                for (k=0; k<2; k++)
                {
                    pstates[k] = nsctrellis.prevState[m+nsctrellis.stateNb*k];//determine previous states
                    inputs[k] = nsctrellis.input[m+nsctrellis.stateNb*k];//determine input
                    C[k] = (inputs[k]?apriori_data[n-1]:0);//compute log of gamma
                }
                for (i=0; i<2; i++)//for each C[i]
                {
                    for (k=0; k<r; k++)
                    {
                        for (j=0; j<Nc; j++)
                        {
                            C[i] += nsctrellis.output[pstates[i]+inputs[i]*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+(n-1)*Nc*r];
                        }
                    }
                }
                A[m+n*nsctrellis.stateNb] = std::max(A[pstates[0]+(n-1)*nsctrellis.stateNb]+C[0], A[pstates[1]+(n-1)*nsctrellis.stateNb]+C[1]);
                sum = std::max(sum, A[m+n*nsctrellis.stateNb]);
            }
            //normalization
            if (std::isinf(sum))
            {
                continue;
            }
            for (m=0; m<nsctrellis.stateNb; m++)
            {
                A[m+n*nsctrellis.stateNb] -= sum;
            }
        }

        //backward recursion
#pragma omp section
        for (n=N-1; n>=0; n--)
        {
            sum = -INFINITY;//normalisation factor
            for (m=0; m<nsctrellis.stateNb; m++) //initial state
            {
                for (k=0; k<2; k++)
                {
                    nstates[k] = nsctrellis.nextState[m+k*nsctrellis.stateNb];//determine next states
                    C[k] = (k?apriori_data[n]:0);//compute log of gamma
                }
                for (i=0; i<2; i++)
                {
                    for (k=0; k<r; k++)
                    {
                        for (j=0; j<Nc; j++)
                        {
                            C[i] += nsctrellis.output[m+i*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                        }
                    }
                }
                B[m+n*nsctrellis.stateNb] = std::max(B[nstates[0]+(n+1)*nsctrellis.stateNb]+C[0], B[nstates[1]+(n+1)*nsctrellis.stateNb]+C[1]);
                sum = std::max(sum, B[m+n*nsctrellis.stateNb]);
            }
            //normalisation
            if (std::isinf(sum))
            {
                continue;
            }
            for (m=0; m<nsctrellis.stateNb; m++)
            {
                B[m+n*nsctrellis.stateNb] -= sum;
            }
        }
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
#pragma omp parallel for private(n,k,sum,sumbis)
    for (n=1; n<(N+1); n++)
    {
        sum = -INFINITY;
        sumbis = -INFINITY;
        for (k=0; k<(nsctrellis.stateNb/2); k++)
        {
            sum = std::max(sum, A[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]+B[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]);//nominator
            sumbis = std::max(sumbis, A[k+n*nsctrellis.stateNb]+B[k+n*nsctrellis.stateNb]);//denominator
        }
        extrinsic_data[n-1] = (sum-sumbis)-apriori_data[n-1];
    }

    //compute extrinsic_coded
    double *sum0 = new double[r];
    double *sum1 = new double[r];
    extrinsic_coded.set_size(N*Nc*r);
    for (n=0; n<N; n++)
    {
        for (k=0; k<r; k++)
        {
            sum0[k] = -INFINITY;
            sum1[k] = -INFINITY;
        }
        for (mp=0; mp<nsctrellis.stateNb; mp++) //initial state
        {
            for (i=0; i<2; i++)
            {
                m = nsctrellis.nextState[mp+i*nsctrellis.stateNb];//final state
                //compute log of sigma
                index = (m>=(nsctrellis.stateNb/2));//0 if input is 0, 1 if input is 1
                C[0] = (index?apriori_data[n]:0);
                for (k=0; k<r; k++)
                {
                    for (j=0; j<Nc; j++)
                    {
                        C[0] += nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                    }
                }
                C[1] = A[mp+n*nsctrellis.stateNb]+C[0]+B[m+(n+1)*nsctrellis.stateNb];//this is only a buffer
                //compute sums
                for (k=0; k<r; k++)
                {
                    if (nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2])
                    {
                        sum1[k] = std::max(sum1[k], C[1]);
                    }
                    else
                    {
                        sum0[k] = std::max(sum0[k], C[1]);
                    }
                }
            }
        }
        for (k=0; k<r; k++)
        {
            for (j=0; j<Nc; j++)
            {
                index = j+k*Nc+n*r*Nc;
                extrinsic_coded[index] = (1-2*double(scrambler_pattern[j]))*(sum1[k]-sum0[k])-intrinsic_coded[index];
            }
        }
    }

    //free memory
    delete[] nsctrellis.output;
    delete[] nsctrellis.nextState;
    delete[] nsctrellis.prevState;
    delete[] nsctrellis.input;
    delete[] A;
    delete[] B;
    delete[] sum0;
    delete[] sum1;
}
}//end namespace tr
