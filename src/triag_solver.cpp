#include <iostream>
#include "matrix.h"


std::vector<double> thomas_meth(tridiag_matrix& matrix, const std::vector<double> &d)
{
    size_t N = d.size();
    std::vector<double> p(N, 0.0);
    std::vector<double> q(N, 0.0);
    std::vector<double> result(N, 0.0);

    for (int i = 1; i < N; i++)
    {
        p[i] = -(matrix.c[i-1]) / (matrix.a[i-1] * p[i-1] + matrix.b[i-1]);
        q[i] = (d[i-1] - matrix.a[i-1] * q[i-1]) / (matrix.a[i-1] * p[i-1] + matrix.b[i-1]);    
    }
    for (int i = N - 1 ; i > -1; i-- ) {
        if ( i != N - 1) {
            result[i] = p[i+1] * result[i+1] + q[i+1];
        }
        else {
            result[i] = (d[i] - matrix.a[i] * q[i]) / (matrix.a[i] * p[i] + matrix.b[i]);
        }
    }
    return result;
}   
