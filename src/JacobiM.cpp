#include "SIMethods.hpp"

std::vector<double> JacobiM(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & tau, const size_t & N, const double & tolerable_error)
{
    std::vector<double> x = x0;
    std::vector<double> temp = b - A*x;
    for(size_t i = 0; i < N; ++i)
    {
        for(size_t j = 0; j < b.size(); ++j)
        {
            x[j] = (temp[j] + x[j])/A(j,j);
        }
        
        temp = b - A*x;
        
        if (vlength(temp) < tolerable_error)
        {
            break;
        }
        
    }
    return x;

}