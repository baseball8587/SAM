#pragma once
#include <Eigen/Sparse>

Eigen::SparseMatrix<double> computeMask(const std::vector<int>& rows, const std::vector<int>& cols, int nz_cnt) {
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nz_cnt);

    int maxRow = 0;
    int maxCol = 0;

    for (int i = 0; i < nz_cnt; ++i) {
        tripletList.emplace_back(rows[i], cols[i], 1.0);
        if (rows[i] > maxRow) maxRow = rows[i];
        if (cols[i] > maxCol) maxCol = cols[i];
    }

    Eigen::SparseMatrix<double> msk(maxRow + 1, maxCol + 1);
    msk.setFromTriplets(tripletList.begin(), tripletList.end());

    return msk;
}