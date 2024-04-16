#include "dense_matrix.hpp"

double Dense_Matrix::operator()(std::size_t i, std::size_t j) const
{
    return V[cols * i + j];
}

Dense_Matrix Dense_Matrix::operator+(const Dense_Matrix & another) const
{
    std::vector<double>result;
    result.reserve(rows * cols);
    for (std::size_t i = 0; i < rows; i++)
    {
        for(std::size_t j = 0; j < cols; j++)
        {
            result.push_back(V[cols * i + j] + another(i, j));
        }
    }
    Dense_Matrix result_Matrix(result, cols, rows);
    return result_Matrix;
}

Dense_Matrix Dense_Matrix::operator*(const Dense_Matrix & another) const
{   

    std::vector<double>result;
    result.reserve(rows * another.GetCols());
    for (std::size_t i = 0; i < rows; i++)
    {
        for(std::size_t j = 0; j < another.GetCols(); j++)
        {
            double sum = 0;

            for(std::size_t k = 0; k < cols; k++)
            {
                sum += V[cols * i + k] * another(k,j);
            }
            result.push_back(sum);
        }
    }
        
    Dense_Matrix result_Matrix(result, rows, another.GetCols());
    return result_Matrix; 
}

Dense_Matrix Dense_Matrix::operator*(const double & num) const
{
    std::vector<double>result;
    result.reserve(rows * cols);
    result = V * num;
    Dense_Matrix result_Matrix(result, rows, cols);
    return result_Matrix;
}

std::vector<double>Dense_Matrix::operator*(const std::vector<double>& vec) const
{
    std::vector <double> result;
    result.reserve(rows);
    for(std::size_t i = 0; i < rows; i++)
    {
        double sum = 0;
        for(std::size_t j = 0; j < cols; j++)
        {
            sum += V[cols * i + j] * vec[j];
        }
        result.push_back(sum);
    }   
    return result;
}

void Dense_Matrix::show() const
{
    for(std::size_t i = 0; i < rows; i++)
    {
        for(std::size_t j = 0; j < cols; j++)
        {
            std::cout << V[cols * i + j] << " ";
        }
        std::cout<< " "<<std::endl;
    }
    std::cout<< " "<<std::endl;
}

size_t Dense_Matrix::GetRows() const
{
    return rows;
}

size_t Dense_Matrix::GetCols() const
{
    return cols;
}
