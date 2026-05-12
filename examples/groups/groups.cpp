/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 * Author: Varga, Steven <steven@vargaconsulting.ca>
 */
#include <h5cpp/all>
#include <cstddef>


int main() {

	auto fd = h5::create("001.h5", H5F_ACC_TRUNC);
	{	// create intermediate groups
		// the default h5::default_lcpl will create immediate path and the encoding is set to utf8
		h5::gr_t gr{H5Gcreate(fd, "my-group/sub/path", static_cast<hid_t>(h5::default_lcpl), H5P_DEFAULT, H5P_DEFAULT)};

		h5::lcpl_t lcpl = h5::char_encoding{H5T_CSET_UTF8} | h5::create_intermediate_group{1};
		H5Gcreate(fd, "/mygroup", static_cast<hid_t>(lcpl), H5P_DEFAULT, H5P_DEFAULT);
	}
	{ // exceptions
		h5::mute(); // mute CAPI error handler
		try { // this group already exists, will throw exception
			H5Gcreate(fd, "/mygroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}catch (  const h5::error::any& e){
			std::cout << e.what() <<"\n";
		}

		// catching all exceptions 
		try { // this group already exists, will throw exception
			H5Gcreate(fd, "/mygroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}catch (  const h5::error::any& e){
			std::cout << e.what() <<"\n";
		}
		h5::unmute(); // re-enable CAPI error handling
	}
	{ // opening a group and adding attributes
		h5::gr_t gr{H5Gopen(fd, "/mygroup", H5P_DEFAULT)};
		h5::awrite(gr, "temperature", 42.0);
		h5::awrite(gr, "unit", "C");
		h5::awrite(gr, "vector of ints", std::vector<int>({1,2,3,4,5}));
		h5::awrite(gr, "initializer list", std::initializer_list({1,2,3,4,5}));
		h5::awrite(gr, "strings", std::initializer_list({"first", "second", "third","..."}));
	}
}
