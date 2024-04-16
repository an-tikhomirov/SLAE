#pragma once
#include "csr_matrix.hpp"
#include <cmath>

std::vector<double> FGD(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const size_t & N, const double & tolerable_error)
{ 
    std::vector<double> x = x0;
    std::vector<double> r_i = A * x - b;
    double tau_i;
    for(size_t i = 0; i < N; ++i)
    {
        x = x - (r_i * r_i * r_i)/(r_i * (A * r_i));
        r_i = A * x - b;
        if (vlength(r_i) < tolerable_error)
        {
            break;
        }
    }

    return x;
}