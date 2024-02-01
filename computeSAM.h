#pragma once
#pragma once
#include "Eigen/Sparse"
#include <iostream>
#include "find.h"
#include <unordered_map>
#include "createSparseSubMatrix.h"

Eigen::SparseMatrix<double> computeSAM(const Eigen::SparseMatrix<double>& As,
    const Eigen::SparseMatrix<double>& PP,
    const Eigen::SparseMatrix<double>& PP2,
    int N) {

    std::cout << "Starting computeSAM function." << std::endl;

    // 1. Check Input Matrices
    if (As.rows() == 0 || As.cols() == 0 || PP.rows() == 0 || PP.cols() == 0 || PP2.rows() == 0 || PP2.cols() == 0) {
        std::cerr << "Error: One or more input matrices are empty." << std::endl;
        return Eigen::SparseMatrix<double>();
    }

    std::vector<Eigen::VectorXi> nz_M(N), nz_LS(N);
    Eigen::VectorXi nnz_M(N), nnz_LS(N);

    int nnzMM = PP2.nonZeros();
    Eigen::VectorXd rowM = Eigen::VectorXd::Zero(2 * nnzMM);
    Eigen::VectorXd colM = Eigen::VectorXd::Zero(2 * nnzMM);
    Eigen::VectorXd valM = Eigen::VectorXd::Zero(2 * nnzMM);

    std::cout << "Processing columns..." << std::endl;
    // printFirst10Elements(PP, "PP");
     //error might be here
    for (int k = 1; k < N; ++k) {

        nz_M[k] = find(PP, k);
        nz_LS[k] = find(PP2, k);



        // 2. Check Non-Zero Elements After find
        if (nz_M[k].size() == 0 || nz_LS[k].size() == 0) {
            std::cerr << "Error: Empty vector at column " << k << std::endl;
            continue;
        }

        nnz_M(k) = nz_M[k].size();
        nnz_LS(k) = nz_LS[k].size();
    }

    int max_col = nnz_M.maxCoeff();
    int max_row = nnz_LS.maxCoeff();


    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(max_row, max_col);
    Eigen::VectorXd M = Eigen::VectorXd::Zero(max_row);
    int cntrM = 0;

    for (int k = 1; k < N; ++k) {
        // Check Submatrix Operations

        if (nz_LS[k].size() == 0 || nz_M[k].size() == 0) {
            std::cerr << "Error: Empty row or column vectors for subMatrix at column " << k << std::endl;
            continue;
        }
        // Declare rows_map at the beginning of the loop
        std::unordered_map<int, int> rows_map;

        // Populate rows_map with the mapping from nz_LS to local indices
        for (int i = 0; i < nz_LS[k].size(); ++i) {
            rows_map[nz_LS[k][i]] = i;
        }

        // Create the submatrix G
        Eigen::SparseMatrix<double> G = createSparseSubMatrix(As, nz_LS[k], nz_M[k]);

        // Create the right-hand side vector for the current column of As
        Eigen::VectorXd As0_col(G.rows());
        int idx = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(As, k); it; ++it) {
            if (rows_map.count(it.row())) {  // Ensure that the row is part of the submatrix
                As0_col(idx++) = it.value();
            }
        }

        // Check Before Solving Linear System
        if (G.rows() != As0_col.size()) {
            std::cerr << "Error: Dimension mismatch before solving linear system at column " << k << std::endl;
            continue;
        }

        // Solve the linear system using SparseQR
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(G);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Decomposition failed at column " << k << std::endl;
            continue;
        }

        Eigen::VectorXd M = solver.solve(As0_col);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solving failed at column " << k << std::endl;
            continue;
        }
        //std::cout << "size " << M.size() << std::endl;
        // Insert the results into the appropriate places in the rowM, colM, and valM arrays
        for (int i = 0; i < M.size(); ++i) {
            rowM(cntrM) = static_cast<double>(nz_M[k](i));
            colM(cntrM) = static_cast<double>(k);
            valM(cntrM) = M(i);
            cntrM++;
        }
    }


    std::cout << "Finished processing columns. Constructing result matrix..." << std::endl;

    Eigen::SparseMatrix<double> result(As.rows(), N);
    for (int i = 0; i < cntrM; ++i) {
        result.insert(rowM(i), colM(i)) = valM(i);
    }

    result.makeCompressed(); // This is important for efficiency

    std::cout << "Finished computeSAM function." << std::endl;

    return result;
}