/*
 * Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
 */
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <h5cpp/all>

int main() {
    auto fd = h5::create("example_xtensor_blas.h5", H5F_ACC_TRUNC);

    xt::xtensor<double, 2> A = {{1, 2}, {3, 4}};
    xt::xtensor<double, 2> B = {{5, 6}, {7, 8}};

    auto C = xt::linalg::dot(A, B);
    h5::write(fd, "xtensor_blas/dot_product", C);

    auto C_read = h5::read<xt::xtensor<double, 2>>(fd, "xtensor_blas/dot_product");

    return 0;
}
