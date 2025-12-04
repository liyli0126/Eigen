//
//  kron_matrix.hpp
//  2D_Possion_Equ
//
//  Created by Yingli Li on 7/8/24.
//

#ifndef kron_matrix_hpp
#define kron_matrix_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

// Function prototype for Kronecker product
Eigen::MatrixXd kroneckerProduct(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

// Function to generate 2D Poisson matrix
Eigen::SparseMatrix<double> generate2DPoissonMatrix(int n);

// Function to generate 3D Poisson matrix
Eigen::SparseMatrix<double> generate3DPoissonMatrix(int n);


#endif /* kron_matrix_hpp */


