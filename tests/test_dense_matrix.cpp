#include <gtest/gtest.h>
#include "../src/dense_matrix.hpp"

TEST(Dense, matrix_mult_num)
{
	std::vector <double> V1 = {50, 3, 7, 99, 8, 21, 37, 1, 9, 11, 31, 17};
	double num = 2;
	std::vector <double> V2 = {100, 6, 14, 198, 16, 42, 74, 2, 18, 22, 62, 34};

	Dense_Matrix M1(V1, 3, 4);
	Dense_Matrix CURRENT_RES(V2, 3, 4);

	Dense_Matrix RES = M1 * num;;

	for (std::size_t i = 0; i < 3; i++)
	{
		for (std::size_t j = 0; j < 4; j++)
		{
			ASSERT_NEAR(RES(i, j), CURRENT_RES(i, j), 0.1);
		}
	}
}

TEST(Dense, matrix_plus_another_matrix)
{
	std::vector <double> V1 = {50, 3, 7, 99, 8, 21, 37, 1, 9, 11, 31, 17};
	std::vector <double> V2 = {50, 97, 93, 1, 92, 79, 63, 99, 91, 89, 69, 83};
	std::vector <double> V3 = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};

	Dense_Matrix M1(V1, 3, 4);
	Dense_Matrix M2(V2, 3, 4);
	Dense_Matrix CURRENT_RES(V3, 3, 4);

	Dense_Matrix RES = M1 + M2;

	for (std::size_t i = 0; i < 3; i++)
	{
		for (std::size_t j = 0; j < 4; j++)
		{
			ASSERT_NEAR(RES(i, j), CURRENT_RES(i, j), 0.1);
		}
	}
}
TEST(Dense, matrix_mult_another_matrix)
{
	std::vector <double> V1 = {3, 7, 1, 5};
	std::vector <double> V2 = {-1, 0, 3, 2};
	std::vector <double> V3 = {18, 14, 14, 10};

	Dense_Matrix M1(V1, 2, 2);
	Dense_Matrix M2(V2, 2, 2);
	Dense_Matrix CURRENT_RES(V3, 2, 2);

	Dense_Matrix RES = M1 * M2;

	for (std::size_t i = 0; i < 2; i++)
	{
		for (std::size_t j = 0; j < 2; j++)
		{
			ASSERT_NEAR(RES(i, j), CURRENT_RES(i, j), 0.1);
		}
	}
}

TEST(Dense, matrix_mult_vec)
{
	std::vector <double> V1 = {3, 7, 1, 5};
	std::vector <double> vec = {1, 0};
	std::vector <double> CURRENT_RES = {3, 1};

	std::vector RES = Dense_Matrix(V1, 2, 2) * vec;

	for(std::size_t i = 0; i < 2; i++)
	{
		ASSERT_NEAR(RES[i], CURRENT_RES[i], 0.5);
	}
}