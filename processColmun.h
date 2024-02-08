#pragma once
#include "Eigen/Sparse"
#include <iostream>
#include <unordered_map>
#include "createSparseSubMatrix.h"

void processColumn(int k,
    const std::vector<Eigen::VectorXi>& nz_LS,
    const std::vector<Eigen::VectorXi>& nz_M,
    const Eigen::SparseMatrix<double>& As,
    Eigen::VectorXd& rowM,
    Eigen::VectorXd& colM,
    Eigen::VectorXd& valM,
    int& cntrM) {
    if (nz_LS[k].size() == 0 || nz_M[k].size() == 0) {
        std::cerr << "Error: Empty row or column vectors for subMatrix at column " << k << std::endl;
        return;
    }

    std::unordered_map<int, int> rows_map;
    for (size_t i = 0; i < nz_LS[k].size(); ++i) {
        rows_map[nz_LS[k][i]] = i;
    }

    Eigen::SparseMatrix<double> G = createSparseSubMatrix(As, nz_LS[k], nz_M[k]);
    Eigen::VectorXd As0_col(G.rows());
    int idx = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(As, k); it; ++it) {
        if (rows_map.count(it.row())) {
            As0_col(idx++) = it.value();
        }
    }

    if (G.rows() != As0_col.size()) {
        std::cerr << "Error: Dimension mismatch before solving linear system at column " << k << std::endl;
        return;
    }

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(G);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed at column " << k << std::endl;
        return;
    }

    Eigen::VectorXd M = solver.solve(As0_col);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed at column " << k << std::endl;
        return;
    }

    for (int i = 0; i < M.size(); ++i) {
        rowM(cntrM) = static_cast<double>(nz_M[k](i));
        colM(cntrM) = static_cast<double>(k);
        valM(cntrM) = M(i);
        cntrM++;
    }
}