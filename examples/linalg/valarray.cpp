/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#include <iostream>
#include <valarray>
#include <h5cpp/all>

int main(){
	h5::fd_t fd = h5::create("valarray.h5",H5F_ACC_TRUNC);
	{ // CREATE - WRITE
		std::valarray<double> v(10);
		for (size_t i = 0; i < v.size(); ++i)
			v[i] = static_cast<double>(i);
		h5::write(fd, "std::valarray", v);
	}
	{ // READ
		std::valarray<double> v = h5::read<std::valarray<double>>(fd, "std::valarray");
		for (size_t i = 0; i < v.size(); ++i)
			std::cout << v[i] << " ";
		std::cout << std::endl;
	}
}
