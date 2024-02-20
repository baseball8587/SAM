#pragma once
Eigen::SparseVector<double> getDiagonalSparse(const Eigen::SparseMatrix<double>& mat) {
    Eigen::SparseVector<double> diagonal(mat.rows());

    // Define the limit based on matrix dimensions
    Eigen::Index limit = (mat.rows() <= mat.cols()) ? mat.rows() : mat.cols();

    for (int k = 0; k < limit; ++k) {
        double value = mat.coeff(k, k);
        if (value != 0) {  // Only insert non-zero values
            diagonal.insert(k) = value;
        }
    }

    return diagonal;
}