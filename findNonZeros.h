#pragma once
#include "Eigen/Sparse"
void findNonZeros(const Eigen::SparseMatrix<double>& As, std::vector<int>& I, std::vector<int>& J) {
    for (int k = 0; k < As.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(As, k); it; ++it) {
            I.push_back(it.row());  // Row index of non-zero entry
            J.push_back(it.col());  // Column index of non-zero entry
        }
    }
}