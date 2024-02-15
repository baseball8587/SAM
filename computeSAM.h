#pragma once
#pragma once
#include "Eigen/Sparse"
#include <iostream>
#include "find.h"
#include <unordered_map>
#include "createSparseSubMatrix.h"
#include <unordered_map>
#include <omp.h>
#include <numeric>

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
    std::vector<Eigen::Triplet<double>> allTriplets; // Container for all triplets
    int maxThreads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> threadLocalTriplets(maxThreads);

#pragma omp parallel
    {
        std::vector<Eigen::Triplet<double>> localTriplets; // Local storage for triplets
        int threadId = omp_get_thread_num();

#pragma omp for nowait // Schedule loop iterations without waiting
        for (int k = 1; k < N; ++k) {
            if (nz_LS[k].size() == 0 || nz_M[k].size() == 0) {
                std::cerr << "Error: Empty row or column vectors for subMatrix at column " << k << std::endl;
                continue;
            }

            std::unordered_map<int, int> rows_map;
            for (int i = 0; i < nz_LS[k].size(); ++i) {
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

            /*Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
            solver.compute(G);
            if (solver.info() == Eigen::Success) {
                Eigen::VectorXd M = solver.solve(As0_col);
                if (solver.info() == Eigen::Success) {*/
                    for (int i = 0; i < As0_col.size(); ++i) {
                        localTriplets.emplace_back(nz_LS[k][i], k, As0_col(i)); // Construct triplet for thread-local matrix
                    }
               // }
            //}
        }

#pragma omp critical
        {
            threadLocalTriplets[threadId] = std::move(localTriplets); // Move localTriplets to the thread-specific vector
        }
    }

    // Sequentially merge the results from each thread
    for (const auto& threadTriplets : threadLocalTriplets) {
        allTriplets.insert(allTriplets.end(), threadTriplets.begin(), threadTriplets.end());
    }

    // Continue with the rest of your code to construct the result matrix
    Eigen::SparseMatrix<double> result(As.rows(), N);
    result.setFromTriplets(allTriplets.begin(), allTriplets.end()); // Construct the result matrix from triplets
    result.makeCompressed(); // Compress the result matrix for efficient storage and computation

    // Optional: Print some triplets for verification
    for (int i = 4326; i <4340 ; i++) {
        const Eigen::Triplet<double>& triplet = allTriplets[i];
        std::cout << "Triplet at index " << i << ": Row: " << triplet.row()
            << ", Column: " << triplet.col() << ", Value: " << triplet.value() << std::endl;
    }

    std::cout << "Finished computeSAM function." << std::endl;
    return result;
}