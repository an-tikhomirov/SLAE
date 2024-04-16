#include "SIMethods.hpp"

std::vector<double> SIM(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & tau, const size_t & N, const double & tolerable_error)
{
    std::vector<double> x = x0;
    std::vector<double> temp = A * x - b;
    for(size_t i = 0; i < N; ++i)
    {
        x = x - tau*temp;
        temp = A * x - b;
        if (vlength(temp) < tolerable_error)
        {
            break;
        }
    }

    return x;
} 