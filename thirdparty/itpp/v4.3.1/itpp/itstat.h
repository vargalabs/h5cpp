/*!
 * \file
 * \brief Include file for the IT++ statistics module
 * \author Adam Piatyszek and Conrad Sanderson
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

#ifndef ITSTAT_H
#define ITSTAT_H

/*!
 * \defgroup stat Statistics Module
 * @{
 */

//! \defgroup histogram Histogram
//! \defgroup statistics Miscellaneous Statistics Functions

/*! \defgroup MOG Mixture of Gaussians (MOG)
    \brief Classes and functions for modelling multivariate data as a Mixture of Gaussians
    \author Conrad Sanderson

    The following example shows how to model data:
    \code
    Array<vec> X;
    // ... fill X with vectors ...
    int K = 3;     // specify the number of Gaussians
    int D = 10;    // specify the dimensionality of vectors
    MOG_diag model(K,D);
    MOG_diag_kmeans(model, X, 10, 0.5, true, true); // initial optimisation using 10 iterations of k-means
    MOG_diag_ML(model, X, 10, 0.0, 0.0, true);      // final optimisation using 10 iterations of ML version of EM
    double avg = model.avg_log_lhood(X);            // find the average log likelihood of X
    \endcode

    See also the tutorial section for a more elaborate example.
*/

/*!
 * @}
 */


#include <itpp/itbase.h>
#include <itpp/stat/histogram.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/stat/mog_generic.h>
#include <itpp/stat/mog_diag.h>
#include <itpp/stat/mog_diag_kmeans.h>
#include <itpp/stat/mog_diag_em.h>

#endif // #ifndef ITSTAT_H
