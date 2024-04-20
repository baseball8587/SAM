#pragma once
void multiplyMatrixByVectorElementWise(Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& vec) {
    // Square every element of vec
    Eigen::VectorXd squaredVec = vec.array().square();

    // Multiply mat by squaredVec element-wise for rows
    for (int i = 0; i < mat.outerSize(); ++i)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
            it.valueRef() *= squaredVec(it.row());
}