#pragma once
#include "csr_matrix.hpp"
#include <cmath>

std::vector<double> SGS(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const size_t & N, const double & tolerable_error)
{
    std::size_t n = x0.size();
    std::vector<double> x = x0;
    for(int k = 0; k < N; ++k)
    {
        for(int i = 0; i < n; ++i)
        {
            double sum = 0;
            for(int j = 0; j < n; ++j)
            {
                if (i != j) {sum += A(i,j) * x[j];}
            }
            x[i] = (b[i] - sum) / A(i,i);
        }
        for(int i = n - 1; i >= 0; --i)
        {
            double sum = 0;
            for(int j = 0; j < n; ++j)
            {
                if (i != j) {sum += A(i,j) * x[j];}
            }
            x[i] = (b[i] - sum) / A(i,i);
        }
    }

    return x;
}