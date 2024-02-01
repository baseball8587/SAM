#pragma once
#include <Eigen/Sparse>
#include <iostream>
#include <windows.h>
#include "readMatrix.h"

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