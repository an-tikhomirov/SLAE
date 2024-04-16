#pragma once
#include "csr_matrix.hpp"

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

// int main()
// {
//     std::vector<double> x0 = {1., 2., 4.};
//     std::vector<double> v = {7., 1., 5., 1.};
//     std::vector<std::size_t> c = {0, 2, 1, 2};
//     std::vector<std::size_t> r = {0, 2, 3, 4};
//     CSR_Matrix A{v, c, r};
//     A.show();
//     std::vector<double> b = {6., 1., 1.};
//     std::vector<double> res = SIM(A, x0, b, 0.25, 7, 0.000001);
//     for(size_t  i = 0; i < b.size(); i++)
//     {
//         std::cout<<res[i]<<std::endl;
//     }
// }