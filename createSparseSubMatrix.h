#pragma once
Eigen::SparseMatrix<double> createSparseSubMatrix(
    const Eigen::SparseMatrix<double>& matrix,
    const Eigen::VectorXi& rows,
    const Eigen::VectorXi& cols) {

    std::vector<Eigen::Triplet<double>> triplets;

    // Reserve space for non-zero elements to avoid reallocations
    triplets.reserve(rows.size() * cols.size());

    // Maps to keep track of the new indices in the submatrix
    std::unordered_map<int, int> rows_map, cols_map;
    for (int i = 0; i < rows.size(); ++i) {
        rows_map[rows[i]] = i;
    }
    for (int i = 0; i < cols.size(); ++i) {
        cols_map[cols[i]] = i;
    }

    // Iterate over the specified columns
    for (int k = 0; k < cols.size(); ++k) {
        // Iterate over all non-zero elements in the current column
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, cols[k]); it; ++it) {
            // If the current row is one of the specified rows, add it to the triplets
            if (rows_map.count(it.row())) {
                triplets.emplace_back(rows_map[it.row()], cols_map[it.col()], it.value());
            }
        }
    }

    // Create the submatrix with the correct size
    Eigen::SparseMatrix<double> submatrix(rows.size(), cols.size());
    submatrix.setFromTriplets(triplets.begin(), triplets.end());

    return submatrix;
}
