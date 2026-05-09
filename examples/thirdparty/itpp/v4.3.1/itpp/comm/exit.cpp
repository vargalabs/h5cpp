/*!
 * \file
 * \brief Implementation of EXtrinsic Information Transfer (EXIT) chart class
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

#include <itpp/comm/exit.h>
#include <itpp/stat/histogram.h> //histogram class for mutual information computation
#include <itpp/base/itcompat.h>

namespace itpp
{
double EXIT::Gaussian_Fct::operator()(double x) const
{
  return (1.0/std::sqrt(_sigma*itpp::m_2pi))*std::exp(-itpp::sqr(x-(_sigma/2.0))/(2.0*_sigma))*::log2(1+std::exp(-x));
}

double EXIT::extrinsic_mutual_info(const itpp::vec &obs, const itpp::bvec &cond, const int &N)
{
    //initialize histogram
    itpp::Histogram<double> hist(itpp::min(obs), itpp::max(obs), N);//common definition interval for both PDFs

    //conditional PDF knowing that a bit of 0 was emitted
    itpp::ivec idx = itpp::find(cond==itpp::bin(0));
    itpp::vec cond_obs = obs(idx);
    hist.reset();//start counting
    hist.update(cond_obs);
    itpp::vec left_pdf = hist.get_pdf();//the pdf is computed without taking into account the interval length (step)
    itpp::ivec left_int = itpp::find(left_pdf!=0);//integration interval for the left PDF

    //conditional PDF knowing that a bit of 1 was emitted
    idx = itpp::find(cond==itpp::bin(1));
    cond_obs = obs(idx);
    hist.reset();//restart counting
    hist.update(cond_obs);
    itpp::vec right_pdf = hist.get_pdf();
    itpp::ivec right_int = itpp::find(right_pdf!=0);//integration interval for the right PDF

    //mutual extrinsic information
    itpp::vec left_half = itpp::elem_mult(left_pdf(left_int), itpp::log2(itpp::elem_div(2.0*left_pdf(left_int), left_pdf(left_int)+right_pdf(left_int))));
    double IE = itpp::sum(left_half)-0.5*(left_half(0)+left_half(left_half.length()-1));//numerical integration without taking into account the inteval length (see conditional PDF computation)
    itpp::vec right_half = itpp::elem_mult(right_pdf(right_int), itpp::log2(itpp::elem_div(2.0*right_pdf(right_int), left_pdf(right_int)+right_pdf(right_int))));
    IE += itpp::sum(right_half)-0.5*(right_half(0)+right_half(right_half.length()-1));//numerical integration
    IE *= 0.5;

    return IE;
}

}//namespace itpp

