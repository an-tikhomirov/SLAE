#pragma once
#include "csr_matrix.hpp"
#include <cmath>

std::vector<double> Chebyshev(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & lambda_min, const double & lambda_max, const size_t & r, const double & tolerable_error)
{
    std::vector<double> x = x0;
    std::vector<double> temp = A * x - b;
    size_t N = pow(2, r); // degree of the Chebyshev polynomial, number of its roots, number of MPI steps
    std::vector <size_t> indx(N, 0);
    std::vector <double> taus(N, 0);
    const double incos = M_PI/2/N;
    indx[0] = 0;
    indx[N/2] = 1;
    size_t mid = N/2;

    for (size_t i = 2; i < r + 1; ++i)
    {
        for(size_t j = 0; j < N; j += mid)
        {
            indx[j + mid/2] = pow(2, i) - 1 - indx[j];
        }
        mid /=2;

    } 
    // std::cout<<"elements of indx"<<std::endl;
    // for(size_t k = 0; k < N; ++k)
    // {
    //     std::cout<<indx[k];
    // }
    
    for(size_t i = 0; i < N; ++i)
    {
        taus[indx[i]] = 2/(lambda_min + lambda_max + (lambda_min - lambda_max) * std::cos(incos * (2 * i + 1)));
    }

    for(size_t i = 0; i < N; ++i)
    {
        x = x - taus[i]*temp;
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
//     std::vector<double> res = Chebyshev(A, x0, b, 1., 7., 7, 0.000001);
//     for(size_t  i = 0; i < b.size(); i++)
//     {
//         std::cout<<res[i]<<std::endl;
//     }
// }