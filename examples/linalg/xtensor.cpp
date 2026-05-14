/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#include <xtensor.hpp>
#include <h5cpp/all>

int main(){
	{ // CREATE - WRITE xtensor
		xt::xtensor<double,2> M = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
		h5::fd_t fd = h5::create("xtensor.h5",H5F_ACC_TRUNC);
		h5::ds_t ds = h5::create<double>(fd,"create then write"
				,h5::current_dims{10,20}
				,h5::max_dims{10,H5S_UNLIMITED}
				,h5::chunk{2,3} | h5::fill_value<double>{3} |  h5::gzip{9}
		);
		h5::write( ds,  M, h5::offset{2,2}, h5::stride{1,3}  );
	}
	{
		xt::xarray<double> V = xt::xarray<double>::from_shape({8});
		for(size_t i=0; i<8; ++i) V(i) = static_cast<double>(i+1);
		// simple one shot write that computes current dimensions and saves vector
		h5::write( "xtensor.h5", "one shot create write",  V);
	}
	{ // READ xtensor
		xt::xtensor<double,2> M = h5::read<xt::xtensor<double,2>>("xtensor.h5","create then write");
	}
	{ // READ xarray
		xt::xarray<double> V = h5::read<xt::xarray<double>>("xtensor.h5","one shot create write");
	}
}
