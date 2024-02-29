#pragma once
#include "Eigen/Sparse"
Eigen::SparseMatrix<double> createSparseMatrix(const std::vector<int>& I, const std::vector<int>& J) {
    Eigen::SparseMatrix<double> findA;

    // Create a list of triplets for setting sparse matrix values
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(I.size());

    for (size_t i = 0; i < I.size(); ++i) {
        tripletList.push_back(Eigen::Triplet<double>(I[i], J[i], 1.0));
    }

    int numRows = *std::max_element(I.begin(), I.end()) + 1;
    int numCols = *std::max_element(J.begin(), J.end()) + 1;

    findA.resize(numRows, numCols);

    // Fill the sparse matrix
    findA.setFromTriplets(tripletList.begin(), tripletList.end());
 
    return findA;
}