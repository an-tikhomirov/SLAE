#include "csr_matrix.hpp"


double CSR_Matrix::operator()(std::size_t i, std::size_t j) const{
    for(std::size_t k = this->row_index[i]; k < this->row_index[i+1]; k++){
        if (col_index[k] == j){
            return this->V[k];
        }
    }
    return 0;
}

CSR_Matrix CSR_Matrix::operator*(const double & num) const
{
    std::vector<double> result = V * num;
    CSR_Matrix result_Matrix(result, col_index, row_index);
    return result_Matrix;
}

std::vector <double> CSR_Matrix::operator*(const std::vector <double> & vec) const
{
    std::vector <double> result;
    result.reserve(row_index.size() - 1);
    for (std::size_t i = 0; i < row_index.size() - 1; i++)
    {
        double sum  = 0;
        for (std::size_t j = row_index[i]; j < row_index[i + 1]; j++)
        {
            sum += V[j] * vec[col_index[j]];
        }
        result.push_back(sum);
    }
    return result;
}

void CSR_Matrix::show() const {
        std::cout << "non-zero values: ";
        for(std::size_t i = 0; i <V.size(); i++) std::cout << V[i] << " ";
        std::cout << std::endl;
        std::cout << "the column in which the corresponding entry V is located: ";
        for(std::size_t i = 0; i < col_index.size(); i++) std::cout << col_index[i] << " ";
        std::cout << std::endl;
        std::cout << "the number of non-zero elements in each string and row_index sums in every string: ";
        for(std::size_t i = 0; i < row_index.size(); i++) std::cout << row_index[i] << " ";
        std::cout << std::endl;    
}

// int main()
// {
//     std::vector <double> V = {1,1};
//     std::vector <double> a = {1,1};
//     std::vector <size_t> col_index = {0, 1};
//     std::vector <size_t> row_index = {0, 1, 2};
//     CSR_Matrix C1(V, col_index, row_index);
//     std::vector <double> result = C1 * a;
//     for(std::size_t i = 0; i < 2; i++)
//     {
//         std::cout << result[i] << std::endl;
//     }


// }