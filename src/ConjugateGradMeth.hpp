#pragma once
#include "csr_matrix.hpp"
#include <cmath>

std::vector<double> CGM(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const size_t & N, const double & tolerable_error)
{
    std::vector<double> x = x0;
    std::vector<double> r = A * x - b;
    std::vector<double> d = r;
    double alpha = 0;
    double beta = 0;

    for(size_t i = 0; i < N; ++i)
    {
        alpha = (r * r)/(d * (A * d));
        x = x - alpha * d;
        std::vector<double> prev_r = r;
        r = A * x - b;
        beta = (r * r)/(prev_r * prev_r);
        d = r + beta * d;
    }
    return x;
}