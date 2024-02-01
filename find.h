#pragma once
#include "Eigen/Sparse"
Eigen::VectorXi find(const Eigen::SparseMatrix<double>& mat, int col) {
    std::vector<int> nonZeros;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, col); it; ++it) {
        nonZeros.push_back(it.row());
    }
    return Eigen::VectorXi::Map(nonZeros.data(), nonZeros.size());
}