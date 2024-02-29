#pragma once
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


    return mat;
}