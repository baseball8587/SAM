#pragma once
#include "Eigen/Sparse"
int size(const Eigen::SparseMatrix<double>& mat, int dim) {
    if (dim == 1) return mat.rows();  // Number of rows
    else if (dim == 2) return mat.cols();  // Number of columns
    else throw std::invalid_argument("Invalid dimension. Use 1 for rows and 2 for columns.");
}