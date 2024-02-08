#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <chrono>
#include <windows.h>
#include <algorithm>
#include <cstdio>
#include <string>
#include <fstream>
#include <exception>
#include <tuple>
#include <omp.h>
#include <unordered_map>
#include <iomanip>  // For std::setw and std::setprecision
#include <deque>
#include "readMatrix.h"
#include "GetDiagnoalSparse.h"
#include "multiplyMatrixByVectorElementWise.h"
#include "scaleMatrixRows.h"
#include "computeIndices.h"
#include "computeMask.h"
#include "findNonZeros.h"
#include "createSparseMatrix.h"
#include "elementWiseMultiply.h"
#include "convertToLogical.h"
#include "find.h"
#include "createSparseSubMatrix.h"
#include "computeSAM.h"
#include "memMap.h"
#include "size.h"





int main() {
    auto start_time_total = std::chrono::high_resolution_clock::now();

    // Open a binary file using Windows-specific API
    HANDLE hFile = CreateFileW(L"C:\\Users\\baseb\\source\\repos\\SAM\\files\\system_00100.bin", GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return -1;
    }

    // Create a file mapping for the opened file
    HANDLE hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!hMapFile) {
        CloseHandle(hFile);
        std::cerr << "Error: Could not create file mapping." << std::endl;
        return -1;
    }

    // Map the file contents into memory
    char* dataPtr = static_cast<char*>(MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0));
    if (!dataPtr) {
        CloseHandle(hMapFile);
        CloseHandle(hFile);
        std::cerr << "Error: Could not map view of file." << std::endl;
        return -1;
    }

    // Read two matrices from the mapped memory and measure the time taken for each
    double elapsed_A, elapsed_B;
    Eigen::SparseMatrix<double> mat_A = readMatrix(dataPtr, elapsed_A);
    Eigen::SparseMatrix<double> mat_B = readMatrix(dataPtr, elapsed_B);

    auto start_copy = std::chrono::high_resolution_clock::now();
    Eigen::SparseMatrix<double> A0 = mat_A;
    auto end_copy = std::chrono::high_resolution_clock::now();
    double elapsed_copy = std::chrono::duration<double>(end_copy - start_copy).count();


    // Clean up mapped memory and file handles
    UnmapViewOfFile(dataPtr);
    CloseHandle(hMapFile);
    CloseHandle(hFile);

    // Print the times taken to read each matrix
    std::cout << "Time taken to read and process matrix A: " << elapsed_A << " seconds" << std::endl;
    std::cout << "Time taken to read and process matrix B: " << elapsed_B << " seconds" << std::endl;
    std::cout << "Time taken to copy: " << elapsed_copy << " seconds" << std::endl;


    // Timing for 'rows' operation
    auto start_time_N = std::chrono::high_resolution_clock::now();
    int N = mat_A.rows();
    auto end_time_N = std::chrono::high_resolution_clock::now();
    double elapsed_time_N = std::chrono::duration<double>(end_time_N - start_time_N).count();

    // Timing for 'diagonal' operation
    auto start_time_d = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd d = mat_A.diagonal();
    auto end_time_d = std::chrono::high_resolution_clock::now();
    double elapsed_time_d = std::chrono::duration<double>(end_time_d - start_time_d).count();
   

    // Timing for 'Dsr1' computation
    auto start_time_Dsr1 = std::chrono::high_resolution_clock::now();

    Eigen::VectorXd Dsr1 = d.array().sqrt().inverse();
 
    auto end_time_Dsr1 = std::chrono::high_resolution_clock::now();
    double elapsed_time_Dsr1 = std::chrono::duration<double>(end_time_Dsr1 - start_time_Dsr1).count();

    // Convert Dsr1 to a diagonal matrix
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagMatrix_Dsr1(Dsr1);


    // Print the elapsed times.
    std::cout << "Time taken for mat_A.rows(): " << elapsed_time_N << " seconds" << std::endl;
    std::cout << "Time taken for mat_A.diagonal(): " << elapsed_time_d << " seconds" << std::endl;
    std::cout << "Time taken for Dsr1 computation: " << elapsed_time_Dsr1 << " seconds" << std::endl;

    auto start_time_Dsr1_Dig = std::chrono::high_resolution_clock::now();
    Eigen::SparseVector<double> dSparse = getDiagonalSparse(mat_A);
    auto end_time_Dsr1_Dig = std::chrono::high_resolution_clock::now();
    double elapsed_time_Dsr1_Dig = std::chrono::duration<double>(end_time_Dsr1_Dig - start_time_Dsr1_Dig).count();
    std::cout << "Time taken for Dsr1 Dig computation: " << elapsed_time_Dsr1_Dig << " seconds" << std::endl;

    auto start_As = std::chrono::high_resolution_clock::now();
    multiplyMatrixByVectorElementWise(mat_A, Dsr1);// mat_A now is As0 so we don't need to copy it
    auto end_time_As = std::chrono::high_resolution_clock::now();
    double elapsed_time_As = std::chrono::duration<double>(end_time_As - start_As).count();
    std::cout << "Time taken for Finding As: " << elapsed_time_As << " seconds" << std::endl;
   
    auto start_Bs = std::chrono::high_resolution_clock::now();

    scaleMatrixRows(mat_B, Dsr1);
    auto end_time_Bs = std::chrono::high_resolution_clock::now();
    double elapsed_time_Bs = std::chrono::duration<double>(end_time_Bs - start_Bs).count();
    std::cout << "Time taken for Finding Bs: " << elapsed_time_Bs << " seconds" << std::endl;


    Eigen::SparseMatrix<double> PP2;// early initilization 
    Eigen::SparseMatrix<double> PP;
    for (int i = 105; i <= 110; i = i + 5) {
        if (i == 105){
        auto matrices = memMap(i);
        std::vector<int> patt = {
       -6603, -6602, -6601, -6600, -6303, -6302, -6301, -6300, -6003,
       -6002, -6001, -6000, -303, -302, -301, -300, -3, -2, -1, 0,
       1, 2, 3, 4, 300, 301, 302, 303, 304, 6000, 6001, 6002, 6003,
       6004, 6300, 6301, 6302, 6303, 6304, 6600, 6601, 6602, 6603, 6604
        };
        auto& matA = matrices.first;
        auto& matB = matrices.second;

        // Use computeIndices to get the IndicesResult
        IndicesResult indicesResult = computeIndices(N, patt);

        // Use the values from the IndicesResult struct to compute the mask
        Eigen::SparseMatrix<double> msk = computeMask(indicesResult.rows, indicesResult.cols, indicesResult.nz_count);

        Eigen::SparseMatrix<double> As;
        std::vector<int> I, J;
        findNonZeros(mat_A, I, J);
        Eigen::SparseMatrix<double> findA = createSparseMatrix(I, J);

        Eigen::SparseMatrix<bool> PP1 = convertToLogical(findA);
        PP = elementWiseMultiply(PP1, msk);
        PP2 = PP;
    }
       else if (i > 105) {
            auto matrices = memMap(i);
            auto& matA = matrices.first;
            auto& matB = matrices.second;

            int N = size(matA, 1);
            Eigen::VectorXd D = matA.diagonal();
            Eigen::VectorXd Dsr_1 = D.array().sqrt().inverse();
            Eigen::SparseVector<double> d_Sparse = getDiagonalSparse(matA);

            multiplyMatrixByVectorElementWise(matA, Dsr_1);
            scaleMatrixRows(matB, Dsr_1);



  

        }
        Eigen::SparseMatrix<double> MM = computeSAM(mat_A, PP, PP2, N);
   

    }
   
    auto end_time_total = std::chrono::high_resolution_clock::now();
    double elapsed_time_total = std::chrono::duration<double>(end_time_total - start_time_total).count();
    std::cout << "Total runtime of main(): " << elapsed_time_total << " seconds" << std::endl;


    return 0;
}



