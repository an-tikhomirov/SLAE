#include <gtest/gtest.h>
#include "../src/SymGaussSeid.hpp"

std::vector<double> V = {8, 3, 5, 1, 4, 7};
std::vector<size_t> col_index = {0, 1, 3, 1, 2, 3};
std::vector<size_t> row_index = {0, 3, 4, 5, 6};
std::vector<double> x0 = {-1, 3, 0.5, 1};
std::vector<double> b = {1, 4, 0, 2};
double N = 40;
double tolerable_error = 0.0001;
CSR_Matrix A(V, col_index, row_index);
std::vector<double> result = SGS(A, x0, b, N, tolerable_error);
std::vector<double> expected_result {-87.0/56.0, 4.0, 0.0, 2.0/7.0};
TEST(FGD, x)
{
    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_NEAR(result[i], expected_result[i], 0.01); }
}