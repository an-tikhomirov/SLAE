#include <gtest/gtest.h>
#include "../src/csr_matrix.hpp"

std::vector <double> V = {7, 5, 1, 3, 2};               // 0 0 0 7
std::vector <size_t> col_index = {3, 0, 1, 2, 1};       // 5 1 0 0
std::vector <size_t> row_index = {0, 1, 3, 4, 5};       // 0 0 3 0 - matrix
CSR_Matrix Cmatrix(V, col_index, row_index);            // 0 2 0 0
TEST(CSR, matrix_mult_num)
{
    double num = 2;
    CSR_Matrix RES = Cmatrix * num;

    std::vector <double> V_new = {14, 10, 2, 6, 4};
    CSR_Matrix CURRENT_RES(V_new, col_index, row_index);

    for (std::size_t i = 0; i < 4; i++)
	{
		for (std::size_t j = 0; j < 4; j++)
		{
			ASSERT_NEAR(RES(i, j), CURRENT_RES(i, j), 0.1);
		}
	}
}

TEST(CSR, matrix_mult_vec)
{
    std::vector <double> V = {1,1};
    std::vector <double> a = {1,1};
    std::vector <size_t> col_index = {0, 1};
    std::vector <size_t> row_index = {0, 1, 2};
    CSR_Matrix Cmatrix1(V, col_index, row_index);
    std::vector <double> RES = Cmatrix1 * a;
    std::vector <double> CURRENT_RES = {1,1};

    for(std::size_t i = 0; i < 2; i++)
    {
        ASSERT_NEAR(RES[i], CURRENT_RES[i], 0.1);
    }

}