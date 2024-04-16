#pragma once
#include "csr_matrix.hpp"

std::vector<double> SIM(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & tau, const size_t & N, const double & tolerable_error); 

std::vector<double> JacobiM(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & tau, const size_t & N, const double & tolerable_error);

// std::vector<double> GaussSeidel(const CSR_Matrix & A, const std::vector<double> & x0, const std::vector<double> & b, const double & tau, const size_t & N, const double & tolerable_error);