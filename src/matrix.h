#include <iostream>
#include <vector>

class tridiag_matrix {
    public: 
        std::vector<double> a;
        std::vector<double> b;
        std::vector<double> c;

        tridiag_matrix(const std::vector<double> &_a, const std::vector<double> &_b, const std::vector<double> &_c) : 
        a(_a), b(_b), c(_c) {};
};


std::vector<double> thomas_meth(tridiag_matrix &matrix, const std::vector<double> &d);




