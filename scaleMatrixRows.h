#pragma once

#include <Eigen/Sparse>


void scaleMatrixRows(Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& vec) {
    // Scale rows of mat
    for (int i = 0; i < mat.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
            it.valueRef() *= vec(it.row());
        }
    }
}