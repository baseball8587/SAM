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



struct IndicesResult {
    std::vector<int> cols;
    std::vector<int> rows;
    int nz_count;
};
// Function to read and return a SparseMatrix from a binary data pointer
Eigen::SparseMatrix<double> readMatrix(char*& dataPtr, double& elapsed_time) {
    auto start = std::chrono::high_resolution_clock::now();

    int len = *reinterpret_cast<int*>(dataPtr);
    dataPtr += sizeof(int);  // Move the data pointer by the size of the integer

    int* i_data = reinterpret_cast<int*>(dataPtr);    // Row indices
    dataPtr += len * sizeof(int);
    int* j_data = reinterpret_cast<int*>(dataPtr);    // Column indices
    dataPtr += len * sizeof(int);
    double* s_data = reinterpret_cast<double*>(dataPtr);  // Values
    dataPtr += len * sizeof(double);

    int max_row = *std::max_element(i_data, i_data + len) + 1; // +1 to convert to size
    int max_col = *std::max_element(j_data, j_data + len) + 1;

    Eigen::SparseMatrix<double> mat(max_row, max_col);
    mat.reserve(len); // Reserve space for upcoming inserts
    for (int index = 0; index < len; ++index) {
        mat.insert(i_data[index], j_data[index]) = s_data[index];
    }

    auto end = std::chrono::high_resolution_clock::now();
    elapsed_time = std::chrono::duration<double>(end - start).count();

    return mat;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> memMap(int i) {
    std::wstring fileName;
    if (i < 100) {
        fileName = L"C:\\Users\\baseb\\source\\repos\\SAM\\files\\system_000" + std::to_wstring(i) + L".bin";
    }
    else {
        fileName = L"C:\\Users\\baseb\\source\\repos\\SAM\\files\\system_00" + std::to_wstring(i) + L".bin";
    }

    std::wcout << L"Attempting to load file: " << fileName << std::endl;

    Eigen::SparseMatrix<double> mat_A0, mat_B0; // Initialize empty matrices

    try {
        HANDLE hFile = CreateFileW(fileName.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (hFile == INVALID_HANDLE_VALUE) {
            std::cerr << "Error: Could not open the file." << std::endl;
            return { mat_A0, mat_B0 }; // Return empty matrices in case of an error
        }

        // Create a file mapping for the opened file
        HANDLE hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
        if (!hMapFile) {
            CloseHandle(hFile);
            std::cerr << "Error: Could not create file mapping." << std::endl;
            return { mat_A0, mat_B0 }; // Return empty matrices in case of an error
        }

        // Map the file contents into memory
        char* dataPtr = static_cast<char*>(MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0));
        if (!dataPtr) {
            CloseHandle(hMapFile);
            CloseHandle(hFile);
            std::cerr << "Error: Could not map view of file." << std::endl;
            return { mat_A0, mat_B0 }; // Return empty matrices in case of an error
        }

        double elapsed_C, elapsed_D;
        mat_A0 = readMatrix(dataPtr, elapsed_C);
        mat_B0 = readMatrix(dataPtr, elapsed_D);

        // Clean up mapped memory and file handles
        UnmapViewOfFile(dataPtr);
        CloseHandle(hMapFile);
        CloseHandle(hFile);

        std::wcout << L"Successfully loaded and processed: " << fileName << std::endl;
    }
    catch (const std::exception& e) {
        std::wcerr << L"Error processing file: " << fileName << L" - " << e.what() << std::endl;
        return { mat_A0, mat_B0 }; // Return empty matrices in case of an error
    }

    return { mat_A0, mat_B0 };
}

// Function to print the first 10 elements of a matrix
void printFirst10Elements(const Eigen::SparseMatrix<double>& mat, const std::string& name) {
    int count = 0;
    std::cout << "First 10 elements from " << name << ":\n";

    // Iterate over the matrix outer dimension (e.g., columns for a column-major matrix)
    for (int k = 0; k < mat.outerSize() && count < 100; ++k) {
        // Inner iterator for matrix entries within the current outer index
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it && count < 100; ++it) {
            std::cout << "(" << it.row()  << "," << it.col()  << "): " << it.value() << "\n";
            count++;
        }
    }
    std::cout << std::endl;
}
int size(const Eigen::SparseMatrix<double>& mat, int dim) {
    if (dim == 1) return mat.rows();  // Number of rows
    else if (dim == 2) return mat.cols();  // Number of columns
    else throw std::invalid_argument("Invalid dimension. Use 1 for rows and 2 for columns.");
}

Eigen::VectorXd diag(const Eigen::SparseMatrix<double>& mat) {
    return mat.diagonal();
}


Eigen::SparseVector<double> getDiagonalSparse(const Eigen::SparseMatrix<double>& mat) {
    Eigen::SparseVector<double> diagonal(mat.rows());

    // Define the limit based on matrix dimensions
    int limit = (mat.rows() <= mat.cols()) ? mat.rows() : mat.cols();

    for (int k = 0; k < limit; ++k) {
        double value = mat.coeff(k, k);
        if (value != 0) {  // Only insert non-zero values
            diagonal.insert(k) = value;
        }
    }

    return diagonal;
}
Eigen::MatrixXd multiplyByDiagonalMatrix(const Eigen::MatrixXd& A, const Eigen::SparseVector<double>& diag) {
    if (A.rows() != diag.size()) {
        std::cerr << "Dimensions do not match for multiplication." << std::endl;
        return Eigen::MatrixXd();  // Return an empty matrix
    }

    Eigen::MatrixXd result(A);
    for (Eigen::SparseVector<double>::InnerIterator it(diag); it; ++it) {
        result.row(it.index()) *= it.value();
    }

    return result;
}


void multiplyMatrixByVectorElementWise(Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& vec) {
    // Square every element of vec
    Eigen::VectorXd squaredVec = vec.array().square();

    // Multiply mat by squaredVec element-wise for rows
    for (int i = 0; i < mat.outerSize(); ++i)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
            it.valueRef() *= squaredVec(it.row());
}




void scaleMatrixRows(Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& vec) {
    // Scale rows of mat
    for (int i = 0; i < mat.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
            it.valueRef() *= vec(it.row());
        }
    }
}
IndicesResult computeIndices(const int N, const std::vector<int>& patt) {
    std::vector<int> cols;
    std::vector<int> rows;
    int nz_cnt = 0;

    // Estimate the maximum size for cols and reserve space
    cols.reserve(N * patt.size());
    rows.reserve(N * patt.size());  // Assuming each pattern entry could theoretically be valid

    for (int k = 1; k <= N; k++) {
        for (size_t i = 0; i < patt.size(); i++) {
            int val = k + patt[i];
            if (val > 0 && val <= N) {
                nz_cnt++;
                cols.push_back(k - 1);
                rows.push_back(val - 1);
            }
        }
    }

    return { cols, rows, nz_cnt };
}


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



void printLastNonZeros(const Eigen::SparseMatrix<double>& msk) {
    std::vector<std::tuple<int, int, double>> allElements;

    // Iterate through non-zero elements
    for (int k = 0; k < msk.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(msk, k); it; ++it) {
            allElements.push_back(std::make_tuple(it.row(), it.col(), it.value()));
        }
    }

}

void findNonZeros(const Eigen::SparseMatrix<double>& As, std::vector<int>& I, std::vector<int>& J) {
    for (int k = 0; k < As.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(As, k); it; ++it) {
            I.push_back(it.row());  // Row index of non-zero entry
            J.push_back(it.col());  // Column index of non-zero entry
        }
    }
}
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

Eigen::SparseMatrix<bool> convertToLogical(const Eigen::SparseMatrix<double>& findA) {
    Eigen::SparseMatrix<bool> PP1(findA.rows(), findA.cols());
    PP1.reserve(findA.nonZeros());

    // Iterate over the non-zero elements of findA
    for (int k = 0; k < findA.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(findA, k); it; ++it) {
            // Simply set the corresponding entry in PP1 to 'true'
            PP1.coeffRef(it.row(), it.col()) = true;
        }
    }

    return PP1;
}

Eigen::SparseMatrix<double> elementWiseMultiply(const Eigen::SparseMatrix<bool>& PP1, const Eigen::SparseMatrix<double>& msk) {
    if (PP1.rows() != msk.rows() || PP1.cols() != msk.cols()) {
        throw std::runtime_error("Dimensions of matrices do not match for element-wise multiplication.");
    }

    Eigen::SparseMatrix<double> result(PP1.rows(), PP1.cols());

    // Iterate over the non-zero elements of the boolean matrix
    for (int k = 0; k < PP1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<bool>::InnerIterator it1(PP1, k); it1; ++it1) {
            double value = msk.coeff(it1.row(), it1.col());
            if (value != 0) {
                result.coeffRef(it1.row(), it1.col()) = value;
            }
        }
    }

    return result;
}


void printFirst10Elements(const Eigen::VectorXd& vec) {
  
    for (int i = 0; i < 100; ++i) {
        std::cout << "Element " << i + 1 << ": " << vec[i] << std::endl;
    }
}

void printFirst10Elements(const Eigen::VectorXi& vec) {

    for (int i = 0; i < 100; ++i) {
        std::cout << "Element " << i + 1 << ": " << vec[i] << std::endl;
    }
}

Eigen::MatrixXd extractSubmatrix(const Eigen::SparseMatrix<double>& mat,
    const Eigen::VectorXi& rows,
    const Eigen::VectorXi& cols) {
    Eigen::MatrixXd submat(rows.size(), cols.size());
    for (int i = 0; i < rows.size(); ++i) {
        for (int j = 0; j < cols.size(); ++j) {
            submat(i, j) = mat.coeff(rows(i), cols(j));
        }
    }
    return submat;
}
// Helper functions:
Eigen::VectorXi find(const Eigen::SparseMatrix<double>& mat, int col) {
    std::vector<int> nonZeros;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, col); it; ++it) {
        nonZeros.push_back(it.row());
    }
    return Eigen::VectorXi::Map(nonZeros.data(), nonZeros.size());
}

Eigen::MatrixXd subMatrix(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXi& rows, const Eigen::VectorXi& cols) {
    Eigen::MatrixXd result(rows.size(), cols.size());
    for (int i = 0; i < rows.size(); ++i) {
        for (int j = 0; j < cols.size(); ++j) {
            result(i, j) = mat.coeff(rows(i), cols(j));
        }
    }
    return result;
}
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
//this is the sparse one
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

//this is the dense one
//Eigen::SparseMatrix<double> computeSAM(const Eigen::SparseMatrix<double>& As,
//    const Eigen::SparseMatrix<double>& PP,
//    const Eigen::SparseMatrix<double>& PP2,
//    int N) {
//
//    std::cout << "Starting computeSAM function." << std::endl;
//
//    int nnzMM = PP2.nonZeros();
//    std::cout << nnzMM << std::endl;
//    Eigen::VectorXd rowM = Eigen::VectorXd::Zero(2*nnzMM);
//    Eigen::VectorXd colM = Eigen::VectorXd::Zero(2*nnzMM);
//    Eigen::VectorXd valM = Eigen::VectorXd::Zero(2*nnzMM);
//    std::vector<Eigen::VectorXi> nz_M(N), nz_LS(N);
//    Eigen::VectorXi nnz_M(N), nnz_LS(N);
//
//    int cntrM = 0;  // Declaration moved here
//
//    std::cout << "Initialized variables. Processing columns..." << std::endl;
//
//
//    for (int k = 1; k < N; ++k) {
//        std::vector<int> nonZeros_PP, nonZeros_PP2;
//        for (Eigen::SparseMatrix<double>::InnerIterator it(PP, k); it; ++it) {
//            nonZeros_PP.push_back(it.row());
//        }
//        for (Eigen::SparseMatrix<double>::InnerIterator it(PP2, k); it; ++it) {
//            nonZeros_PP2.push_back(it.row());
//        }
//        //this part is wrong
//        nz_M[k] = Eigen::VectorXi::Map(nonZeros_PP.data(), nonZeros_PP.size());
//        nz_LS[k] = Eigen::VectorXi::Map(nonZeros_PP2.data(), nonZeros_PP2.size());
//        nnz_M(k) = static_cast<int>(nonZeros_PP.size());
//        nnz_LS(k) = static_cast<int>(nonZeros_PP2.size());
//
//        if (nnz_LS[k] <= 0 || nnz_M[k] <= 0) {
//            std::cout << "Incompatible dimensions at k = " << k << ". Skipping this iteration." << std::endl;
//            continue;
//        }
//    }
//
//
//    for (int k = 1; k < N; ++k) {
//
//        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(nnz_LS[k], nnz_M[k]);
//        for (int i = 0; i < nnz_LS[k]; ++i) {
//            for (int j = 0; j < nnz_M[k]; ++j) {
//                G(i, j) = As.coeff(nz_LS[k](i), nz_M[k](j));
//            }
//        }
//
//
//        Eigen::VectorXd As0_col = Eigen::VectorXd::Zero(As.rows());
//        for (Eigen::SparseMatrix<double>::InnerIterator it(As, k); it; ++it) {
//            As0_col(it.row()) = it.value();
//        }
//
//        Eigen::VectorXd sub_vector;
//        try {
//            sub_vector = G.colPivHouseholderQr().solve(As0_col.head(nnz_LS[k]));
//        }
//        catch (const std::exception& e) {
//            std::cout << "Error during matrix operations at k = " << k << ": " << e.what() << ". Skipping this iteration." << std::endl;
//            continue;
//        }
//
//        for (int i = 0; i < nnz_M[k]; ++i) {
//            rowM(cntrM) = static_cast<double>(nz_M[k](i));
//            colM(cntrM) = static_cast<double>(k);
//            valM(cntrM) = sub_vector(i);
//            cntrM++;
//        }
//    }
//
//    std::cout << "Finished processing columns. Constructing result matrix..." << std::endl;
//
//    Eigen::SparseMatrix<double> result(As.rows(), N);
//    for (int i = 0; i < cntrM; ++i) {
//        result.insert(static_cast<int>(rowM(i)), static_cast<int>(colM(i))) = valM(i);
//    }
//
//    std::cout << "Finished computeSAM function." << std::endl;
//
//    return result;
//}
//


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



            Eigen::SparseMatrix<double> MM = computeSAM(mat_A, PP, PP2, N);

            // Set the output formatting
            std::cout << std::fixed << std::setprecision(4);

            //// Iterate over each element in the matrix
            //for (int i = 0; i < MM.rows(); ++i) {
            //    for (int j = 0; j <2; ++j) {
            //        double value = MM.coeff(i, j);
            //        std::cout << "(" << i << "," << j << ") " << std::setw(16) << value << std::endl;
            //    }
            //}


           // printFirst10Elements(MM, "MM");

        }

    }
    auto end_time_total = std::chrono::high_resolution_clock::now();
    double elapsed_time_total = std::chrono::duration<double>(end_time_total - start_time_total).count();
    std::cout << "Total runtime of main(): " << elapsed_time_total << " seconds" << std::endl;


    return 0;
}



