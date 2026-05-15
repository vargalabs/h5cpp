# Copyright (c) 2018-2026 Steven Varga, Toronto,ON Canada
if(H5CPP_HAVE_ROS3_VFD)
    add_executable(s3 s3.cpp)
    set_target_properties(s3 PROPERTIES
        CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
    target_link_libraries(s3 PRIVATE h5cpp)
    message(STATUS "H5CPP: building S3/ROS3 example")
else()
    message(STATUS "H5CPP: skipping S3/ROS3 example (ROS3 VFD not available)")
endif()
