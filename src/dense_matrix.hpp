#pragma once
#include "vector_operations.hpp"

class Dense_Matrix
{
    private:
    std::vector <double> V; // values
    std::size_t rows, cols; // size of matrix

    public:
    Dense_Matrix(const std::vector <double> V, std::size_t rows, std::size_t cols): V(V), rows(rows), cols(cols) {} // constructor

    double operator()(std::size_t i, std::size_t j) const; //  access to i and j elements
    
    Dense_Matrix operator+(const Dense_Matrix & another) const; // matrix + another matrix

    Dense_Matrix operator*(const Dense_Matrix & another) const; // matrix * another matrix

    Dense_Matrix operator*(const double & num) const; // matrix * number

    std::vector <double> operator*(const std::vector<double> & vec) const; // matrix * vector

    void show() const; // just shows elements of dense matrix

    size_t GetRows() const; // Getter for i

    size_t GetCols() const; // Getter for j




};