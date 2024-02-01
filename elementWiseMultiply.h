#pragma once
#include <Eigen/Sparse>
Eigen::SparseMatrix<double> elementWiseMultiply(const Eigen::SparseMatrix<bool>& PP1, const Eigen::SparseMatrix<double>& msk) {
    if (PP1.rows() != msk.rows() || PP1.cols() != msk.cols()) {
        throw std::runtime_error("Dimensions of matrices do not match for element-wise multiplication.");
    }

    Eigen::SparseMatrix<double> result(PP1.rows(), PP1.cols());

    // Iterate over the non-zero elements of the boolean matrix
    for (int k = 0; k < PP1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<bool>::InnerIterator it1(PP1, k); it1; ++it1) {
            double value = msk.coeff(it1.row(), it1.col());
            if (value != 0) {
                result.coeffRef(it1.row(), it1.col()) = value;
            }
        }
    }

    return result;
}