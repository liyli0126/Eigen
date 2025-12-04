// LinSys.h
// Solvers for dense/full and sparse linear systems
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#ifndef LINSYS_H
#define LINSYS_H

#include <vector>

#include "matrix.h"
#include "vector.h"

// namespace LinLite {
  Vector slvFullLinSysGEPP(const FullMatrix &A, const Vector &b);

  Vector slvFullLowerTrigSys(const FullMatrix &L, const Vector &b);

  Vector slvFullUpperTrigSys(const FullMatrix &U, const Vector &b);

  Vector slvFullSpdSysCholesky(const FullMatrix &A, const Vector &b);

  void slvFullSpdSysCG(Vector &x, const FullMatrix &A, const Vector &b, 
    int &itr, int maxitr, double threshold, double tol);

  Vector slvDiagSys(const DiagMatrix &A, const Vector &b);

  int slvSpaLinSysBiCGStab(Vector &x, const SparseMatrix &A, 
    const Vector &b, const DiagMatrix &B, 
    int &maxitr, double &tol, double atol, int printit);

  void slvSpaSpdSysCG(Vector &x, const SparseMatrix &A, const Vector &b, 
    int &itr, int maxitr, double threshold, double tol);

  void slvSpaBlkSpdSysCG(Vector &x, const SparseBlockMatrix &A, const Vector &b,
    int &itr, int maxitr, double threshold, double tol);

  void slvBlkDiagSpdSysCG(Vector &x, const BlockDiagMatrix &A, 
    const Vector &b);

  void slvBlkDiagSpdSysCholesky(Vector &x, const BlockDiagMatrix &A, 
    const Vector &b);

  void slvSpaBlkSpdSysPCG(Vector &x, const SparseBlockMatrix &A, 
    const Vector &b, const BlockDiagMatrix &B, 
    int &itr, int maxitr, double threshold, double tol);

  void slvBlkDiagSysGEPP(Vector &x, const BlockDiagMatrix &A, 
    const Vector &b);

  int slvSpaBlkLinSysBiCGStab(Vector &x, const SparseBlockMatrix &A, 
    const Vector &b, const BlockDiagMatrix &B, 
    int &maxitr, double &tol, double atol, int printit);

  int slvSpaBlkSysPreCondGMRES(Vector &x, const Vector &b, 
    const SparseBlockMatrix &A, const BlockDiagMatrix &B, 
    int m, int &maxitr, double &atol, double &tol, int printit);

  // int PreCondGMRES(Vector &x, const SparseBlockMatrix &A, 
  //   const Vector &b, const BlockDiagMatrix &B, 
  //   int &maxitr, int m, double &tol, double &atol, int printit);

  // Solving a linear system based on Block-Diagonal-Schur-complement & CG
  // Assuming compatible patterns, mainly for weak Galerkin FEMs
  int slvBDSchurCG(Vector &x0, Vector &x1,
                   const BlockDiagMatrix &A00, const SparseBlockMatrix &A01,
                   const SparseBlockMatrix &A10, const SparseBlockMatrix &A11,
                   const Vector &b0, const Vector &b1,
                   int &itr, int maxitr, double threshold, double tol);
// }

#endif  // LINSYS_H
