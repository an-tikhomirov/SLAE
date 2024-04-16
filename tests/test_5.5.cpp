#include <iostream>
#include <gtest/gtest.h>
#include "../src/triag_solver.cpp"

std::vector<double> d = {2.0, 2.0, 2.0, 2.0, 2.0};

TEST(random_5x5, HandlesAnswer) {

    std::vector<double> b = {15.0, 12.0, 19.0, 10.0, 2.0};
    std::vector<double> a = {0.0, 4.5, 3.0, 2.0, -1.0};
    std::vector<double> c = {2.0, -2.0, 5.5, 7.0, 0.0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result {3236.0/28797.0, 1509.0/9599.0, 1882.0/9599.0, (-3834.0)/9599.0, 7682.0/9599.0};

    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_NEAR(result[i], expected_result[i], 0.000001); }

};
