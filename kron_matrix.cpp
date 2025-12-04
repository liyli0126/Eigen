
//  kron_matrix.cpp
//  2D_Possion_Equ
//
//  Created by Yingli Li on 7/8/24.
//

#include "kron_matrix.hpp"


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

// Compute the Kronecker product of two matrices
Eigen::MatrixXd kroneckerProduct(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    int rowsA = static_cast<int>(A.rows());
    int colsA = static_cast<int>(A.cols());
    int rowsB = static_cast<int>(B.rows());
    int colsB = static_cast<int>(B.cols());

    Eigen::MatrixXd result(rowsA * rowsB, colsA * colsB);

    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsA; ++j) {
            result.block(i * rowsB, j * colsB, rowsB, colsB) = A(i, j) * B;
        }
    }

    return result;
}
// Function to generate 1D Poisson matrix
Eigen::MatrixXd generate1DPoissonMatrix(int n) {
    // Initialize subblock matrix T1D (1D Poisson matrix)
    Eigen::MatrixXd T1D(n, n);
    T1D.setZero(); // Initialize with zeros
    T1D.diagonal().setConstant(2); // Set diagonal elements to 2
    T1D.diagonal(1).setConstant(-1); // Set super-diagonal elements to -1
    T1D.diagonal(-1).setConstant(-1); // Set sub-diagonal elements to -1
    
    // Initialize identity matrix I1D
    Eigen::MatrixXd I1D = Eigen::MatrixXd::Identity(n, n);
    
    std::cout << "Whole block matrix T1D:" << std::endl << T1D << std::endl;
    
    return T1D;
}
    
// Function to generate 2D Poisson matrix
Eigen::SparseMatrix<double> generate2DPoissonMatrix(int n) {
    // Initialize subblock matrix T1D (1D Poisson matrix)
    Eigen::MatrixXd T1D(n, n);
    T1D.setZero(); // Initialize with zeros
    T1D.diagonal().setConstant(2); // Set diagonal elements to 2
    T1D.diagonal(1).setConstant(-1); // Set super-diagonal elements to -1
    T1D.diagonal(-1).setConstant(-1); // Set sub-diagonal elements to -1

    // Initialize identity matrix I1D
    Eigen::MatrixXd I1D = Eigen::MatrixXd::Identity(n, n);

   
    // Construct the whole block matrix T2D
    Eigen::MatrixXd denseT2D = kroneckerProduct(T1D, I1D) + kroneckerProduct(I1D, T1D);
    
    Eigen::SparseMatrix<double> T2D = denseT2D.sparseView();

   // std::cout << "Whole block matrix T2D:" << std::endl << T2D << std::endl;

    return T2D;
}

Eigen::SparseMatrix<double> generate3DPoissonMatrix(int n) {
    // Initialize subblock matrix Tn
    Eigen::MatrixXd T1D(n, n);
    T1D.setZero(); // Initialize with zeros
    T1D.diagonal().setConstant(2); // Set diagonal elements to 2
    T1D.diagonal(1).setConstant(-1); // Set super-diagonal elements to -1
    T1D.diagonal(-1).setConstant(-1); // Set sub-diagonal elements to -1

    // Initialize identity matrix In
    Eigen::MatrixXd I1D = Eigen::MatrixXd::Identity(n, n);

    // Construct the whole block matrix T2D
    Eigen::MatrixXd T2D = kroneckerProduct(T1D, I1D) + kroneckerProduct(I1D, T1D);
    
    // Initialize identity matrix I2D
    Eigen::MatrixXd I2D = Eigen::MatrixXd::Identity(n*n, n*n);
    
    // Construct the whole block matrix T3D
        Eigen::MatrixXd denseT3D = kroneckerProduct(T2D, I1D) + kroneckerProduct(I2D, T1D);
        
        // Convert dense matrix to sparse matrix
        Eigen::SparseMatrix<double> T3D = denseT3D.sparseView();

     
   // std::cout << "Whole block matrix T3D:" << std::endl << T3D << std::endl;


    return T3D;
   
}

