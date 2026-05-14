/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#pragma once

// xtensor-blas does not introduce new container types;
// it operates on xt::xarray and xt::xtensor which are handled by H5Mxtensor.hpp.
// This header ensures xtensor-blas compatibility.
#if defined(XTENSOR_BLAS_HPP) || defined(H5CPP_USE_XTENSOR_BLAS)
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
// No additional specializations needed — H5Mxtensor.hpp covers the types.
#endif
