#pragma once
#include "Eigen/Sparse"
Eigen::SparseMatrix<bool> convertToLogical(const Eigen::SparseMatrix<double>& findA) {
    Eigen::SparseMatrix<bool> PP1(findA.rows(), findA.cols());
    PP1.reserve(findA.nonZeros());

    // Iterate over the non-zero elements of findA
    for (int k = 0; k < findA.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(findA, k); it; ++it) {
            // Simply set the corresponding entry in PP1 to 'true'
            PP1.coeffRef(it.row(), it.col()) = true;
        }
    }

    return PP1;
}