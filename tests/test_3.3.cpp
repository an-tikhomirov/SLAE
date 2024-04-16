#include <iostream>
#include <gtest/gtest.h>
#include "../src/triag_solver.cpp"

std::vector<double> d = {7.0, 7.0, 7.0};

TEST(diagonal_3x3, HandlesAnswer) {

    std::vector<double> a = {0, 0, 0};
    std::vector<double> b = {9, 9, 9};
    std::vector<double> c = {0, 0, 0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = {7.0 / 9.0, 7.0 / 9.0, 7.0 / 9.0};
    
    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_DOUBLE_EQ(result[i], expected_result[i]); }
    
};


TEST(symmetric_3x3, HandlesAnswer) {

    std::vector<double> b = {8.0, 8.0, 8.0};
    std::vector<double> a = {0, 2.0, 2.0};
    std::vector<double> c = {2.0, 2.0, 0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = {3.0 / 4.0, 1.0 / 2.0, 3.0 / 4.0};
    
    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_DOUBLE_EQ(result[i], expected_result[i]); }

};

TEST(random_3x3, HandlesAnswer) {
    
    std::vector<double> b = {8.0, 8.0, -9.0};
    std::vector<double> a = {0, -1.0, 6.0};
    std::vector<double> c= {4.0, 3.0, 0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = {7.0 / 18.0, 35.0 / 36.0, -7.0 / 54.0};

    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_DOUBLE_EQ(result[i], expected_result[i]); }
};

TEST(random_3x3_harder, HandlesAnswer) {
    
    std::vector<double> b = {8.2, 8.3, 9.2};
    std::vector<double> a = {0, 1.6, 5.7};
    std::vector<double> c = {4.6, 2.9, 0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = { 0.510648, 0.611453, 0.382034};

    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_NEAR(result[i], expected_result[i], 0.000001); }

};

TEST(random_3x3_BigNumbers, HandlesAnswer) {

    std::vector<double> b = {88878.0, 83765.0, 925000.0};
    std::vector<double> a = {0.0, -168.0, 57965.0};
    std::vector<double> c = {4679.0, 2.9, 0.0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = { 0.0000743524, 0.0000837162, 0.00000232232};

    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_NEAR(result[i], expected_result[i], 0.000001); }

};


TEST(random_3x3_IncorrectMatrix, HandlesAnswer) {

    std::vector<double> b = {8.0, 8.0, 1.0};
    std::vector<double> a = {0.0, -1.0, 6.0};
    std::vector<double> c = {4.0, 3.0, 0.0};

    tridiag_matrix matrix(a, b, c);

    std::vector<double> result = thomas_meth(matrix, d);
    std::vector<double> expected_result = { 35.0 / 106.0, 231.0 / 212.0, -49.0 / 106.0};

    for (std::size_t i = 0; i < result.size(); i++) { EXPECT_NE(result[i], expected_result[i]); }

};

