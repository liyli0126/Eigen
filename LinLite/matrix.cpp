// matrix.cpp
// Source code for (double) full, sparse, and sparse block matrices
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "vector.h"
using namespace std;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// FullMatrix: Another constructor

FullMatrix::FullMatrix(int numRows, int numCols)
{
  m = numRows;
  n = numCols;
  p = new double[m*n];
  for (int k=0; k<m*n; ++k)  p[k]=0;
}


// FullMatrix: Copy constructor

FullMatrix::FullMatrix(const FullMatrix &A) 
{
  m = A.m;
  n = A.n;
  p = new double[m*n];
  for (int k=0; k<m*n; ++k)  p[k] = A.p[k];
}


// FullMatrix: Copy assignment

FullMatrix &FullMatrix::operator=(const FullMatrix &A) 
{
  if (this!=&A) {
    delete[] p;
    m = A.m;
    n = A.n;
    p = new double[m*n];
    for (int k=0; k<m*n; ++k)  p[k] = A.p[k];
  }
  return *this;
}


// FullMatrix: Mathmatical style access

double &FullMatrix::operator()(int i, int j) const 
{
  return p[(i-1)*n+(j-1)];  // Row-wise storage 
}


// FullMatrix: resize to a zero matrix of a given size

void FullMatrix::resize(int numRows, int numCols) 
{
  delete[] p;
  m = numRows;
  n = numCols;
  p = new double[m*n];
  for (int k=0; k<m*n; ++k)  p[k]=0;
  return;
}


// FullMatrix: zero out a specified row

void FullMatrix::zeroutRow(int rowPos) 
{
  for (int j=1; j<=n; ++j) 
    p[(rowPos-1)*n+(j-1)] = 0;
  return; 
}


// FullMatrix: zero out a specified column

void FullMatrix::zeroutCol(int colPos) 
{
  for (int i=1; i<=m; ++i) 
    p[(i-1)*n+(colPos-1)] = 0;
  return; 
}


// FullMatrix: save to a data file 

void FullMatrix::save2file(char *filename) const 
{
  FILE *fp;

  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit.");
    exit(-1);
  }

  fprintf(fp, "%3d  %3d\n", m, n);

  for (int i=1; i<=m; ++i) 
    for (int j=1; j<=n; ++j) 
      fprintf(fp, "%3d  %3d  %+12.6f\n", i, j, p[(i-1)*n+(j-1)]);

  fclose(fp);
  return;
}


// FullMatrix: cout to screen 

std::ostream &operator<<(std::ostream &strm, const FullMatrix &A) 
{
  strm << A.m << "   " << A.n << "\n";
  for (int i=1; i<=A.m; ++i) {
    for (int j=1; j<=A.n; ++j) {
      // strm << A(i,j) << "   ";
      printf("%+9.6G  ", A(i,j));
    }
    // cout << "\n";
    printf("\n");
  }
  return strm;
}


// FullMatrix: submatrix 

FullMatrix submatrix(const FullMatrix &A, int bgnRow, int bgnCol, 
  int numRows, int numCols) 
{
  FullMatrix B(numRows, numCols);

  for (int i=1; i<=numRows; ++i) 
    for (int j=1; j<=numCols; ++j) 
      B(i,j) = A(bgnRow+i-1, bgnCol+j-1);

  return B;
}


// FullMatrix: transpose 

FullMatrix transpose(const FullMatrix &A) 
{
  FullMatrix B(A.n, A.m);

  for (int i=1; i<=A.m; ++i) 
    for (int j=1; j<=A.n; ++j) 
      B(j,i) = A(i,j);

  return B;
}


// Inverting a full matrix by GEPP
// Finished by GH;  Revised by JKL

FullMatrix inverseGEPP(const FullMatrix &A)
{
  int i, im, j, k, m, n;
  double cmax, lier, tmp;
  double *p;
  
  m = A.rowSize();
  n = A.columnSize();

  if (m!=n) {
    std::cout << "Not a square matrix! Exit!\n";
    exit(EXIT_FAILURE);
  }
  
  // Initializing the inverse matrix and setting it to the identity matrix
  FullMatrix B(m,n);
  for (i=1; i<=n; ++i)  B(i,i) = 1;
  
  // Filling vector p with entries from A (p is in format [row*n+col])
  p = new double[n*n];
  for (i=1; i<=n; ++i)
    for (j=1; j<=n; ++j)
      p[(i-1)*n+(j-1)] = A(i,j);
  
  
  for (j=1; j<=(n-1); ++j) {
    cmax = fabs(p[(j-1)*n+(j-1)]);
    im = j;
    
    for (i=j+1; i<=n; ++i) {
      if (fabs(p[(i-1)*n+(j-1)])>cmax) {
        // im is the index of the maximal column element
        im = i;
        // Updating cmax now
        // GH20160802 This needs to be updated in LinSys GEPP code also
        cmax = p[(i-1)*n+(j-1)];
      }
    }
    
    // if the column absolute max is 0 then the system cannot be solved
    if (cmax==0) {
      std::cout << "Singular matrix!  Exit!\n";
      std::cout << "Column " << j << "\n";
      exit(EXIT_FAILURE);
    }
    
    for (k=j; k<=n; ++k) {  // Pivot rows in A
      tmp = p[(j-1)*n+(k-1)];
      p[(j-1)*n+(k-1)] = p[(im-1)*n+(k-1)];
      p[(im-1)*n+(k-1)] = tmp;
    }
    
    for (k=j; k<=n; ++k) {  // Pivot rows in B to mirror A
      tmp = B(j,k);
      B(j,k) = B(im,k);
      B(im,k) = tmp;
    }
    
    for (i=j+1; i<=n; ++i) {  // Performing elimination on rows
      lier = p[(i-1)*n+(j-1)]/p[(j-1)*n+(j-1)];
      
      
      for(k=1; k<=n; ++k) {  // Element by element row operations
        B(i,k) -= lier*B(j,k);
        p[(i-1)*n+(k-1)] -= lier*p[(j-1)*n+(k-1)];
      }
    }
  }
  
  // Note: At this point p should be the original matrix
  // in the reduced echelon form (REF)
  
  for (i=n; i>1; --i) {  // i counts columns backwards (zeroing column i)
    for (j=1; j<i; ++j) {  // j counts rows
      lier = p[(j-1)*n + (i-1)]/p[(i-1)*n+(i-1)];
      
      for(k=1; k<=n; ++k) {  // Performing row operation on matrix B
        B(j,k) -= lier*B(i,k);
      }
      p[(j-1)*n + (i-1)] = 0;
    }
    
    for(k=1; k<=n; ++k) {  // Dividing row by constant to get identity
      B(i,k) = B(i,k)/p[(i-1)*n+(i-1)];
    }
    p[(i-1)*n+(i-1)]=1;
  }
  
  for(k=1; k<=n; ++k) {  // Divide row 1 by the last constant to finish
    B(1,k) = B(1,k)/p[0];
  }
  p[0] = 1;
  
  delete[] p;
  
  return B;
}


// FullMatrix addition

FullMatrix operator+(const FullMatrix &A, const FullMatrix &B) 
{
  FullMatrix C(A.m, A.n);

  if (A.m!=B.m || A.n!=B.n) {
    std::cout << "Bad matrix sizes!\n";
    exit(EXIT_FAILURE);
  }

  for (int i=1; i<=A.m; ++i) 
    for (int j=1; j<=A.n; ++j) 
      C(i,j) = A(i,j) + B(i,j);

  return C;
}


// FullMatrix multiplication

FullMatrix operator*(const FullMatrix &A, const FullMatrix &B)
{
  FullMatrix C(A.m, B.n);
  
  if (A.n!=B.m) {
    std::cout << "Matrix sizes mismatch!\n";
    exit(EXIT_FAILURE);
  }
  
  for (int i=1; i<=A.m; ++i) {
    for (int j=1; j<=B.n; ++j) {
      C(i,j) = 0;
      for (int k=1; k<=B.m; ++k) {
        C(i,j) = C(i,j) + A(i,k)*B(k,j);
      }
    }
  }
  
  return C;
}


// FullMatrix tensor product: the same results as kronecker product

FullMatrix tensorProduct(const FullMatrix &A, const FullMatrix &B)
{
  int i, ia, ib, j, ja, jb;
  FullMatrix C((A.m)*(B.m), (A.n)*(B.n));
  
  for (ia=1; ia<=A.m; ++ia) {
    for (ib=1; ib<=B.m; ++ib) {
      i = (ia-1)*(B.m) + ib;
      for (ja=1; ja<=A.n; ++ja) {
        for (jb=1; jb<=B.n; ++jb) {
          j = (ja-1)*(B.n) + jb;
          C(i,j) = A(ia,ja) * B(ib,jb);
        }
      }
    }
  }
  
  return C;
}


// FullMatrix kronecker product: the same result as the above tensor product

FullMatrix kron(const FullMatrix &A, const FullMatrix &B)
{
  int i, ia, ib, j, ja, jb;
  FullMatrix C((A.m)*(B.m), (A.n)*(B.n));
  
  for (ia=1; ia<=A.m; ++ia) {
    for (ib=1; ib<=B.m; ++ib) {
      i = (ia-1)*(B.m) + ib;
      for (ja=1; ja<=A.n; ++ja) {
        for (jb=1; jb<=B.n; ++jb) {
          j = (ja-1)*(B.n) + jb;
          C(i,j) = A(ia,ja) * B(ib,jb);
        }
      }
    }
  }
  
  return C;
}


// Scalar-FullMatrix multiplication

FullMatrix operator*(double a, const FullMatrix &A) 
{
  FullMatrix B(A.m, A.n);

  for (int i=1; i<=A.m; ++i) 
    for (int j=1; j<=A.n; ++j) 
      B(i,j) = a*A(i,j);

  return B;
}


// FullMatrix: Taking out a row vector

Vector rowVector(const FullMatrix &A, int i)
{
  Vector vec(A.columnSize());
  for (int j=1; j<=A.columnSize(); ++j)  vec(j) = A(i,j);
  return vec;
}


// FullMatrix: Taking out a column vector

Vector columnVector(const FullMatrix &A, int j)
{
  Vector vec(A.rowSize());
  for (int i=1; i<=A.rowSize(); ++i)  vec(i) = A(i,j);
  return vec;
}


// * overloaded: FullMatrix-Vector multiplication

Vector operator*(const FullMatrix&A, const Vector&v) 
{
  if (A.n!=v.size()) {
    std::cout << "Matrix vector sizes do not match!\n";
    exit(EXIT_FAILURE);
  }

  Vector u(A.m);

  for (int i=1; i<=A.m; ++i) 
    for (int j=1; j<=A.n; ++j) 
      u(i) += A(i,j)*v(j);

  return u;
}


// FullMatrix: colon product
// Finished by RC

double colonProduct(const FullMatrix &A, const FullMatrix &B)
{
  if (A.m!=B.m || A.n!=B.n) {
    std::cout << "Matrix sizes do not match!\n";
    exit(EXIT_FAILURE);
  }
  
  double dp = 0.0;

  for (int i=1; i<=A.m; ++i) 
    for (int j=1; j<=A.n; ++j) 
      dp += A(i,j)*B(i,j);

  return dp;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// DiagMatrix: Another constructor

DiagMatrix::DiagMatrix(int order)
{
  n = order;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i]=0;
}


// DiagMatrix: Copy constructor

DiagMatrix::DiagMatrix(const DiagMatrix &A) 
{
  n = A.n;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i]=A.p[i];
}


// DiagMatrix: Copy assignment

DiagMatrix &DiagMatrix::operator=(const DiagMatrix &A) 
{
  if (this!=&A) {
    delete[] p;
    n = A.n;
    p = new double[n];
    for (int i=0; i<n; ++i)  p[i] = A.p[i];
  }
  return *this;
}


// DiagMatrix: resize to the zero matrix of a given size
// almost like a constructor

void DiagMatrix::resize(int order)
{
  delete[] p;
  n = order;
  p = new double[n];
  for (int i=0; i<n; ++i)  p[i] = 0;
  return;
}


// DiagMatrix: determinant

double DiagMatrix::det() const 
{
  double d = 1.0;
  for (int i=0; i<n; ++i)  d *= p[i];
  return d;
}


// DiagMatrix: save to a data file 

void DiagMatrix::save2file(char *filename) const 
{
  FILE *fp;

  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit!");
    exit(-1);
  }

  fprintf(fp, "%8d\n", n);

  for (int i=0; i<n; ++i)
    fprintf(fp, "%8d  %21.9f\n", i+1, p[i]);

  fclose(fp);
  return;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// SparseMatrix: A constructor with not-much info

SparseMatrix::SparseMatrix(int numRows, int numCols, int numNonZeros)
{
  m = numRows;
  n = numCols;
  nnz = numNonZeros;

  rowoff = new int[m+1];
  colind = new int[nnz];
  val = new double[nnz];

  for (int k=0; k<nnz; ++k)  val[k] = 0;
}


// SparseMatrix: Copy constructor

SparseMatrix::SparseMatrix(const SparseMatrix &A) 
{
  m = A.m;  n = A.n;  nnz= A.nnz;

  rowoff = new int[m+1];
  colind = new int[nnz];
  val = new double[nnz];

  for (int i=0; i<=m; ++i)  rowoff[i] = A.rowoff[i];

  for (int k=0; k<nnz; ++k) {
    colind[k] = A.colind[k];
    val[k] = A.val[k];
  }
}


// SparseMatrix: Copy assignment

SparseMatrix &SparseMatrix::operator=(const SparseMatrix &A) 
{
  if (this!=&A) {
    delete[] rowoff, colind, val;

    m = A.m;  n = A.n;  nnz= A.nnz;

    rowoff = new int[m+1];
    colind = new int[nnz];
    val = new double[nnz];

    for (int i=0; i<=m; ++i)  rowoff[i] = A.rowoff[i];

    for (int k=0; k<nnz; ++k) {
      colind[k] = A.colind[k];
      val[k] = A.val[k];
    }
  }

  return *this;
}


// SparseMatrix: Mathematical style access: entry (i,j), 1<=i<=m, 1<=j<=n

double SparseMatrix::entry(int i, int j) const
{
  double a=0;
  for (int k=rowoff[i-1]; k<rowoff[i]; ++k) {
    if (j==colind[k]) {
      a = val[k];
      break;
    }
  }
  return a;
}


// SparseMatrix: set entry (i,j) to a, 1<=i<=m, 1<=j<=n

void SparseMatrix::setEntry(int i, int j, double a)
{
  for (int k=rowoff[i-1]; k<rowoff[i]; ++k) {
    if (j==colind[k]) {
      val[k] = a;
      break;
    }
  }
  return;
}


// SparseMatrix: add a to entry (i,j), 1<=i<=m, 1<=j<=n

void SparseMatrix::addEntry(int i, int j, double a)
{
  for (int k=rowoff[i-1]; k<rowoff[i]; ++k) {
    if (j==colind[k]) {
      val[k] += a;
      break;
    }
  }
  return;
}


// SparseMatrix: unitfy a row (1<=i<=m)

void SparseMatrix::unitfyRow(int i)
{
  for (int k=rowoff[i-1]; k<rowoff[i]; ++k) {
    if (i==colind[k])  {val[k] = 1.0;}
    else  {val[k]=0.0;}
  }
  return;
}


// SparseMatrix: zero out a row (1<=i<=m)

void SparseMatrix::zeroutRow(int i)
{
  for (int k=rowoff[i-1]; k<rowoff[i]; ++k)  val[k] = 0;
  return;
}


// SparseMatrix: resize to a zero matrix of specified sparsity pattern
// almost like a constructor

void SparseMatrix::resize(int numRows, int numCols,
                          int *numEntrsInRow, int **colIdx)
{
  int i, j, k;
  delete[] rowoff;
  delete[] colind;
  delete[] val;

  m = numRows;
  n = numCols;

  nnz = 0;
  for (i=0; i<m; ++i)  nnz += numEntrsInRow[i];
 
  rowoff = new int[m+1];
  colind = new int[nnz];
  val = new double[nnz];

  k = 0;
  for (i=0; i<m; ++i) {
    rowoff[i] = k;
    for (j=0; j<numEntrsInRow[i]; ++j) {
       colind[k] = colIdx[i][j];
       val[k] = 0;
       k++;
    }
  }
  rowoff[m] = nnz;

  return;
}


// SparseMatrix: save to a data file 

void SparseMatrix::save2file(char *filename) const 
{
  FILE *fp;
  int i, j, k;

  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit.");
    exit(-1);
  }

  fprintf(fp, "%d  %d  %d\n", m, n, nnz);
  fprintf(fp, "\n");

  for (i=1; i<=m; ++i) {
    for (k=rowoff[i-1]; k<rowoff[i]; ++k) {
      j = colind[k];
      fprintf(fp, "%d  %d  %d  %+18.12f\n", i, j, k+1, val[k]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return;
}


// + overloaded: addition of two sparse matrices of the same sparsity pattern

SparseMatrix operator+(const SparseMatrix &A, const SparseMatrix &B) 
{
  SparseMatrix C(A);
  for (int k=0; k<A.nnz; ++k)  C.val[k] = A.val[k] + B.val[k];
  return C;
}


// * overloaded: scalar-SparseMatrix multiplication

SparseMatrix operator*(double a, const SparseMatrix &A) 
{
  SparseMatrix B(A);
  for (int k=0; k<A.nnz; ++k)  B.val[k] = a*A.val[k];
  return B;
}


// * overloaded: SparseMatrix-vector multiplication

Vector operator*(const SparseMatrix &A, const Vector &v) 
{
  Vector u(A.m);

  for (int i=1; i<=A.m; ++i) {
    for (int k=A.rowoff[i-1]; k<A.rowoff[i]; ++k)
      u(i) += A.val[k]*v(A.colind[k]);
  }

  return u;
}


// Get the diagonal of a sparse square matrix

DiagMatrix diagonal(const SparseMatrix &A) 
{
  if (A.m!=A.n) {
    std::cout << "Not a square matrix!\n";
    exit(-1);
  }

  DiagMatrix D(A.n);

  for (int i=1; i<=A.n; ++i) {
    for (int k=A.rowoff[i-1]; k<A.rowoff[i]; ++k) {
      if (i==A.colind[k]) {
        D.setEntry(i, A.val[k]);
        break;
      }
    }
  }

  return D;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// BlockDiagMatrix: Copy constructor

BlockDiagMatrix::BlockDiagMatrix(const BlockDiagMatrix &A) 
{
  n = A.n;
  nb = A.nb;

  db = new int[nb];
  bb = new int[nb];
  bv = new FullMatrix[nb];

  for (int ib=0; ib<nb; ++ib) {
    db[ib] = A.db[ib];
    bb[ib] = A.bb[ib];
    bv[ib] = A.bv[ib];
  }
}


// BlockDiagMatrix: Copy assignment

BlockDiagMatrix &BlockDiagMatrix::operator=(const BlockDiagMatrix &A) 
{
  if (this!=&A) {
    delete[] db, bb, bv;

    n = A.n;
    nb = A.nb;

    db = new int[nb];
    bb = new int[nb];
    bv = new FullMatrix[nb];

    for (int ib=0; ib<nb; ++ib) {
      db[ib] = A.db[ib];
      bb[ib] = A.bb[ib];
      bv[ib] = A.bv[ib];
    }
  }

  return *this;
}


// BlockDiagMatrix: resize to a zero matrix

void BlockDiagMatrix::resize(int numBlks, int *dimBlk) 
{
  int ib;
  delete[] db, bb, bv;

  nb = numBlks;
  db = new int[nb];
  bb = new int[nb];
  bv = new FullMatrix[nb];

  for (ib=0; ib<nb; ++ib)  db[ib] = dimBlk[ib];

  bb[0] = 0;
  for (ib=1; ib<nb; ++ib)  bb[ib] = bb[ib-1]+db[ib-1];

  n = 0;
  for (ib=0; ib<nb; ++ib)  n += db[ib];

  for (ib=0; ib<nb; ++ib)  bv[ib] = FullMatrix(db[ib],db[ib]);

  return;
}


// BlockDiagMatrix: save to a data file

void BlockDiagMatrix::save2file(char *filename) const 
{
  FILE *fp;
  int ii, jj, k;

  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit.");
    exit(-1);
  }

  fprintf(fp, "%4d  %6d\n", nb, n);

  for (k=0; k<nb; ++k) {
    fprintf(fp, "%4d   %2d\n", k+1, db[k]);
    for (ii=1; ii<=bv[k].rowSize(); ++ii) {
      for (jj=1; jj<=bv[k].colSize(); ++jj) {
        fprintf(fp, "  %+9.6f", bv[k](ii,jj));
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
  return;
}


// Overloading * for BlockDiagMatrix * Vector

Vector operator*(const BlockDiagMatrix &D, const Vector &v)
{
  // Assuming their sizes match
  int n = D.size();
  Vector u(n);
  Vector wu, wv;
  FullMatrix A;
  for (int ib=1; ib<=D.nb; ++ib) {
    A = D.block(ib);
    wv = v.getSubVector(D.beginPositionBlock(ib)+1, D.dimBlock(ib));
    wu = A * wv;
    u.setSubVector(D.beginPositionBlock(ib)+1, wu);
  }
  return u;
}


// Cholesky factorization for a block diagonal SPD matrix
// obtain a block diagonal of lower trigular matrices

BlockDiagMatrix CholeskyBlockDiagSPD(const BlockDiagMatrix &A) 
{
  int i, j, k, m, n;
  double sum1, sum2;

  BlockDiagMatrix B(A);
  FullMatrix C;

  for (m=1; m<=A.numBlocks(); ++m) {
    C = A.block(m);
    n = C.colSize();
    FullMatrix L(n,n);

    for (j=1; j<=n; ++j) {
      sum1 = 0.0;
      for (k=1; k<=(j-1); ++k)  sum1 += L(j,k)*L(j,k);
      L(j,j) = sqrt(C(j,j)-sum1);
      for (i=j+1; i<=n; ++i) {
        sum2 = 0.0;
        for (k=1; k<=(j-1); ++k)  sum2 += L(i,k)*L(j,k);
        L(i,j) = (C(i,j)-sum2)/L(j,j);
      }
    }

    B.setBlock(m,L);
  }

  return B;
}


// Inverting a block diagonal matrix by inverting each block using GEPP

BlockDiagMatrix inverseGEPP(const BlockDiagMatrix &A)
{
  BlockDiagMatrix B(A);
  int nb = B.nb;
  for (int i=0; i<nb; ++i)  B.bv[i] = inverseGEPP(B.bv[i]);
  return B;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// SparseBlockMatrix: Default empty constructor

SparseBlockMatrix::SparseBlockMatrix() 
{
  m = 0;  n = 0;  nnz = 0;
  nrb = 0;  ncb = 0;  nnzb = 0;
  drb = 0;  dcb = 0;
  brb = 0;  bcb = 0;
  bro = 0;  bci = 0;
  bv = 0;
}


// SparseBlockMatrix: Copy constructor

SparseBlockMatrix::SparseBlockMatrix(const SparseBlockMatrix &A)
{
  m = A.m;
  n = A.n;
  nnz = A.nnz;
  
  nrb = A.nrb;
  ncb = A.ncb;
  nnzb = A.nnzb;

  drb = new int[nrb];
  dcb = new int[ncb];

  brb = new int[nrb];
  bcb = new int[ncb];
  
  bro = new int[nrb+1];
  bci = new int[nnzb];
  
  bv = new FullMatrix[nnzb];

  for (int i=0; i<nrb; ++i)  drb[i] = A.drb[i];
  for (int j=0; j<ncb; ++j)  dcb[j] = A.dcb[j];

  for (int i=0; i<nrb; ++i)  brb[i] = A.brb[i];
  for (int j=0; j<ncb; ++j)  bcb[j] = A.bcb[j];

  for (int i=0; i<=nrb; ++i)  bro[i] = A.bro[i];
  for (int k=0; k<nnzb; ++k)  bci[k] = A.bci[k];
  
  for (int k=0; k<nnzb; ++k)  bv[k] = A.bv[k];
}


// SparseBlockMatrix: Copy assignment

SparseBlockMatrix &SparseBlockMatrix::operator=(const SparseBlockMatrix &A)
{
  if (this!=&A) {
    delete[] drb, dcb;
    delete[] brb, bcb;
    delete[] bro, bci;
    delete[] bv;

    m = A.m;
    n = A.n;
    nnz = A.nnz;
    
    nrb = A.nrb;
    ncb = A.ncb;
    nnzb = A.nnzb;
    
    drb = new int[nrb];
    dcb = new int[ncb];
    
    brb = new int[nrb];
    bcb = new int[ncb];
    
    bro = new int[nrb+1];
    bci = new int[nnzb];
    
    bv = new FullMatrix[nnzb];
    
    for (int i=0; i<nrb; ++i)  drb[i] = A.drb[i];
    for (int j=0; j<ncb; ++j)  dcb[j] = A.dcb[j];
    
    for (int i=0; i<nrb; ++i)  brb[i] = A.brb[i];
    for (int j=0; j<ncb; ++j)  bcb[j] = A.bcb[j];
    
    for (int i=0; i<=nrb; ++i)  bro[i] = A.bro[i];
    for (int k=0; k<nnzb; ++k)  bci[k] = A.bci[k];

    for (int k=0; k<nnzb; ++k)  bv[k] = A.bv[k];
  }
  
  return *this;
}


// SparseBlockMatrix:  Destructor

SparseBlockMatrix::~SparseBlockMatrix()
{
  delete[] drb, dcb;
  delete[] brb, bcb;
  delete[] bro, bci;
  delete[] bv;
}


// SparseBlockMatrix: The FullMatrix at block (ib,jb)
// Finished by RC 

FullMatrix SparseBlockMatrix::block(int ib, int jb) const 
{
  int bp;  // block position
  for (int k=bro[ib-1]; k<bro[ib]; ++k) {
    if (bci[k] == jb) {
      bp = k;
      break;
    } 
  }
  return bv[bp];
}


// SparseBlockMatrix:  The internal storage location for block (ib,jb)
// Finished by GH;  Revised by JKL
/*
int SparseBlockMatrix::location(int ib, int jb) const
{
  int loc = -1;
  for(int k=bro[ib-1]; k<bro[ib]; ++k)
  {
    if (bci[k]==jb)  {loc = k;  break;}
  }
  return loc;
}
*/


// SparseBlockMatrix: set block (ib,jb) to FullMatrix A

void SparseBlockMatrix::setBlock(int ib, int jb, const FullMatrix &A)
{
  int bp;  // block position

  for (int k=bro[ib-1]; k<bro[ib]; ++k) {
    if (bci[k] == jb) {
      bp = k;
      break;
    } 
  }

  if (bv[bp].rowSize()!=A.rowSize() || bv[bp].colSize()!=A.colSize())  
    bv[bp].resize(A.rowSize(), A.colSize());

  bv[bp] = A;
  return;
}


// SparseBlockMatrix: add FullMatrix A to block (ib,jb)

void SparseBlockMatrix::addBlock(int ib, int jb, const FullMatrix &A)
{
  int bp;  // block position

  for (int k=bro[ib-1]; k<bro[ib]; ++k){
    if (bci[k] == jb) {
      bp = k;
      break;
    } 
  }

  if (bv[bp].rowSize()!=A.rowSize() || bv[bp].colSize()!=A.colSize())  
    bv[bp].resize(A.rowSize(), A.colSize());

  bv[bp] = bv[bp] + A;
  return;
}


// SparseBlockMatrix:  zerout a specified row

void SparseBlockMatrix::zeroutRow(int ib, int iib) 
{
  FullMatrix A;

  A = block(ib,ib);
  A.zeroutRow(iib);
  setBlock(ib, ib, A);

  for (int k=bro[ib-1]; k<bro[ib]; ++k) {
    if (bci[k]==ib)  continue;
    A = bv[k];
    A.zeroutRow(iib);
    bv[k] = A;
  }

  return;
}


// SparseBlockMatrix:  zerout a specified column

void SparseBlockMatrix::zeroutColumn(int jb, int jjb)
{
  FullMatrix A;

  A = block(jb,jb);
  A.zeroutCol(jjb);
  setBlock(jb,jb,A);

  for (int i=1; i<=nrb; ++i) {
    if (i==jb)  continue;
    for (int k=bro[i-1]; k<bro[i]; ++k) {
      if (bci[k]==jb) {
        A = bv[k];
        A.zeroutCol(jjb);
        setBlock(i,jb,A);
      }
    }
  }

  return;
}


// SparseBlockMatrix:  unitify a specified row
// Usually used for enforcing Dirichlet boundary conditions

void SparseBlockMatrix::unitifyRow(int ib, int iib) 
{
  FullMatrix A;

  A = block(ib,ib);
  A.zeroutRow(iib);
  A(iib,iib) = 1; 
  setBlock(ib, ib, A);

  for (int k=bro[ib-1]; k<bro[ib]; ++k) {
    if (bci[k]==ib)  continue;
    A = bv[k];
    A.zeroutRow(iib);
    bv[k] = A;
  }

  return;
}


// SparseBlockMatrix:  unitify a specified column

void SparseBlockMatrix::unitifyColumn(int jb, int jjb)
{
  FullMatrix A;

  A = block(jb,jb);
  A.zeroutCol(jjb);
  A(jjb,jjb) = 1; 
  setBlock(jb,jb,A);

  for (int i=1; i<=nrb; ++i) {
    if (i==jb)  continue;
    for (int k=bro[i-1]; k<bro[i]; ++k) {
      if (bci[k]==jb) {
        A = bv[k];
        A.zeroutCol(jjb);
        setBlock(i,jb,A);
      }
    }
  }

  return;
}


// SparseBlockMatrix:
// resize a sparse block matrix by setting blocks to zero full matrices

void SparseBlockMatrix::resize(int numRowBands, int numColumnBands,
                               int *dimRowBand, int *dimColumnBand,
                               int *numBlocksInRowBand, int **columnBandIndex)
{
  int i, j, k;

  // Removing the existing data/info
  delete[] drb, dcb, brb, bcb, bro, bci, bv;
  m = 0;  n = 0;  nnz = 0;
  nrb = 0;  ncb = 0;  nnzb = 0;

  // Now new data/info
  nrb = numRowBands;
  ncb = numColumnBands;
  
  drb = new int[nrb];
  dcb = new int[ncb];
  
  for (i=0; i<nrb; ++i)  drb[i] = dimRowBand[i];
  for (j=0; j<ncb; ++j)  dcb[j] = dimColumnBand[j];

  brb = new int[nrb];  // Beginning position of row band
  bcb = new int[ncb];  // Beginning position of column band

  brb[0] = 1;  // Yes starts at 1
  bcb[0] = 1;  // Yes starts at 1
  
  for (i=1; i<nrb; ++i)  brb[i] = brb[i-1] + drb[i-1];
  for (j=1; j<ncb; ++j)  bcb[j] = bcb[j-1] + dcb[j-1];

  nnzb = 0;
  for (i=0; i<numRowBands; ++i)  nnzb += numBlocksInRowBand[i];

  bro = new int[nrb+1];
  bci = new int[nnzb];
  bv = new FullMatrix[nnzb];

  for (k=0; k<nnzb; ++k)  bci[k] = 0;

  k = 0;  // The "offset" always starts at 0
  for (i=0; i<nrb; ++i) {
    bro[i] = k;
    for (j=0; j<numBlocksInRowBand[i]; ++j) {
      bci[k] = columnBandIndex[i][j];
      k++;
    }
  }
  bro[nrb] = nnzb;
  
  m = 0;
  for (i=0; i<nrb; ++i)  m += drb[i];

  n = 0;
  for (j=0; j<ncb; ++j)  n += dcb[j];

  nnz = 0;
  for (i=0; i<nrb; ++i) {
    for (k=bro[i]; k<bro[i+1]; ++k) {
      if (bci[k]!=0) {
        nnz += drb[i] * dcb[bci[k]-1];
        bv[k].resize(drb[i], dcb[bci[k]-1]);  // Each block is a zero full mat.
      }
    }
  }

  return;
}


// SparseBlockMatrix: save to a data file

void SparseBlockMatrix::save2file(char *filename) const 
{
  FILE *fp;
  int i, ii, jj, k;

  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit.");
    exit(-1);
  }

  fprintf(fp, "%4d  %4d  %6d\n", m, n, nnz);
  fprintf(fp, "%4d  %4d  %6d\n", nrb, ncb, nnzb);

  for (i=0; i<nrb; ++i) {
    for (k=bro[i]; k<bro[i+1]; ++k) {
      fprintf(fp, "%4d   %4d   ", i+1, bci[k]);
      fprintf(fp, "%2d  %2d  ", brb[i], bcb[bci[k]-1]);
      fprintf(fp, "%2d  %2d\n", bv[k].rowSize(), bv[k].colSize());
      for (ii=1; ii<=bv[k].rowSize(); ++ii) {
        for (jj=1; jj<=bv[k].colSize(); ++jj) {
          fprintf(fp, "  %9.6f", bv[k](ii,jj));
        }
        fprintf(fp, "\n");
      }
    }
  }

  fclose(fp);
  return;
}


// JL20150603: TO BE CHECKED!!!
// SparseBlockMatrix: cout for a small size SparseBlockMatrix
std::ostream &operator<<(std::ostream &strm, const SparseBlockMatrix &A)
{
/*
 int i;

  cout << "nrb: " << A.nrb << "\n";
  cout << "ncb: " << A.ncb << "\n";
  cout << "nnzb: " << A.nnzb << "\n";

  cout << "drb: [ ";
  for (i=0; i<A.nrb; ++i)
    strm << A.drb[i] << " ";
  cout << "]" <<"\n";

  cout << "dcb: [ ";
  for (i=0; i<A.ncb; ++i)
    strm << A.dcb[i] << " ";
  cout << "]" << "\n";

  cout << "bro: [ ";
  for (i=0; i<=A.nrb; ++i)
    strm << A.bro[i] << " ";
  cout << "]" << "\n";

  cout << "bci: [ ";
  for (i=0; i<A.nnzb; ++i)
    strm << A.bci[i] << " ";
  cout << "]" << "\n";

  cout << "brb: [ ";
  for (i=0; i<=A.nrb; ++i)
    strm << A.brb[i] << " ";
  cout << "]" << "\n";

  cout << "bcb: [ ";
  for (i=0; i<=A.ncb; ++i)
    strm << A.bcb[i] << " ";
  cout << "]" << "\n";

  cout << "bv: " << "\n"; 
  for (i=0; i<A.nnzb; ++i)
    strm << A.bv[i];
  cout << "\n";

  cout << "m: "  << A.m << "\n";
  cout << "n: "  << A.n << "\n";
  cout << "nnz: "  << A.nnz << "\n";
*/
  return strm;
}


// * overloaded: scalar-SparseBlockMatrix multiplication

SparseBlockMatrix operator*(double a, const SparseBlockMatrix &A)
{
  SparseBlockMatrix B(A);
  for (int k=0; k<B.nnzb; ++k)  B.bv[k] = a * A.bv[k];
  return B;
}


// + overloaded: SparseBlockMatrix-SparseBlockMatrix addition
// Assuming they have the same pattern

SparseBlockMatrix operator+(const SparseBlockMatrix &A,
                            const SparseBlockMatrix &B)
{
  SparseBlockMatrix C(A);
  for (int k=0; k<C.nnzb; ++k)  C.bv[k] = A.bv[k] + B.bv[k];
  return C;
}


// * overloaded: SparseBlockMatrix-Vector multiplication
// Finished by RC;  Modified by JKL(20160925)

Vector operator*(const SparseBlockMatrix &A, const Vector &v)
{
  int i, j, k, l, uPos, vPos;

  Vector u(A.m);
  for (i=0; i<A.nrb; ++i) {
    for (k=1; k<=A.drb[i]; ++k) {
      uPos = A.brb[i] - 1 + k;
      u(uPos) = 0;
      for (j=A.bro[i]; j<A.bro[i+1]; ++j) {
        for (l=1; l<=A.dcb[A.bci[j]-1]; ++l) {
          vPos = A.bcb[A.bci[j]-1] - 1 + l;
          u(uPos) += A.bv[j](k,l) * v(vPos);
        }
      }
    }
  }
  /*
  std::cout << "i=" << i << "  ";
  std::cout << "A.brb[i]=" << A.brb[i] << "\n";
  std::cout << "k=" << k << "\n";
  std::cout << "uPos=" << uPos << "\n";
  std::cout << "l=" << l << "  " << std::flush;
  std::cout << "vPos=" << vPos << "\n" << std::flush;
  std::cout << "\n";
  */
  return u;
}


// Get the block diagonal matrix from a sparse block matrix
// assuming all the diagonal blocks are square
// not really a friend, SparseBlockMatrix public interfaces used

BlockDiagMatrix diagonal(const SparseBlockMatrix &A) 
{
  int ib;
  BlockDiagMatrix B;

  int numBlks = A.numRowBlocks();
  int *dimBlk = new int[numBlks];

  for (ib=1; ib<=numBlks; ++ib) 
    dimBlk[ib-1] = A.block(ib,ib).rowSize();

  B.resize(numBlks, dimBlk);

  for (ib=1; ib<=numBlks; ++ib)  
    B.setBlock(ib, A.block(ib,ib));

  return B;
}


// B=D*A: BlockDiagMatrix * SparseBlockMatrix
// Assuming D, A have compatible sparsity patterns
// Finished by GH;  Revised by JKL

SparseBlockMatrix operator*(const BlockDiagMatrix &D,
                            const SparseBlockMatrix &A)
{
  if (D.size() != A.m) {
    std::cout << "Matrix sizes mismatch!\n";
    exit(EXIT_FAILURE);
  }
  if (D.numberBlocks() != A.nrb) {
    std::cout << "Band counts incompatible!\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<A.nrb; ++i) {
    if (D.dimBlock(i+1) != A.drb[i]) {
      std::cout << "Band size mismatch!\n";
      exit(EXIT_FAILURE);
    }
  }
  
  // Sparsity pattern preserved
  SparseBlockMatrix B(A);

  // Looping through each row and multiply by the block across the row
  int k = 0;
  for (int i=0; i<B.nrb; ++i) {
    for (int j=B.bro[i]; j<B.bro[i+1]; ++j) {
      B.bv[k] = D.block(i+1) * A.bv[k];
      k++;
    }
  }

  return B;
}


// B=A*D: SparseBlockMatrix * BlockDiagMatrix
// Assuming A, D have compatible sparsity patterns
// Finished by GH;  Revised by JKL

SparseBlockMatrix operator*(const SparseBlockMatrix &A,
                            const BlockDiagMatrix &D)
{
  // std::cout << "A.n=" << A.n << "\n";
  // std::cout << "D.size()=" << D.size() << "\n";

  if (A.n != D.size()) {
    std::cout << "Matrix sizes mismatch!\n";
    exit(EXIT_FAILURE);
  }
  if ( A.ncb != D.numberBlocks()) {
    std::cout << "Band counts incompatible!\n";
    exit(EXIT_FAILURE);
  }
  for (int j=0; j<A.ncb; ++j) {
    if (A.dcb[j] != D.dimBlock(j+1)) {
      std::cout << "Band size mismatch!\n";
      exit(EXIT_FAILURE);
    }
  }
  
  // Sparsity pattern preserved
  SparseBlockMatrix B(A);
  
  // Looping through each row
  for (int i=1; i<=B.nrb; ++i) {  // The i-th row band of B
    for (int k=B.bro[i-1]; k<B.bro[i]; ++k) {
      int j = B.bci[k];  // (i,j)->k, the j-th column band of B
      if (j==i) {
        B.bv[k] = A.bv[k] * D.block(j);
      }
    }
  }

  return B;
}


// "SparseBlockMatrix * SparseBlockMatrix" nontrivial, to be implemented later
// Main difficulty: Efficiently setting up sparsity pattern of the product
/*
friend SparseBlockMatrix operator*(const SparseBlockMatrix &A,
                                   const SparseBlockMatrix &B);
*/


// JL20161109: NOT WORKING YET: SC lost symmetry
// BDSchur: Block Diagonal Schur complement (useful for WGFEMs or alike)
// Schur complement = A11 - A10 * A00^{-1} * A01
// Assuming compatible sparsity patterns:
// Furthermore A10 sparsity pattern = transpose of A01 sparsity pattern
//  input: A00, A01, A10, A11
// Output: SC, A00inv, A10hat

int BDSchur(SparseBlockMatrix &SC,
            BlockDiagMatrix &A00inv, SparseBlockMatrix &A10hat,
            const BlockDiagMatrix &A00, const SparseBlockMatrix &A01,
            const SparseBlockMatrix &A10, const SparseBlockMatrix &A11)
{
  A00inv = inverseGEPP(A00);
  A10hat = A10 * A00inv;

  // JL20160923: NOW using old-style loops; LATER ON using sets
  // std::set<int> IdxSet1, IdxSet2, IdxSet3;
  FullMatrix A, B, C;
  SC = A11;  // Schur complement sparsity pattern and initial values set to A11
  for (int i=1; i<=SC.nrb; ++i) {  // The i-th row band of SC
    for (int k=SC.bro[i-1]; k<SC.bro[i]; ++k) {
      int j = SC.bci[k];  // The j-th column band of SC
      SC.bv[k] = A11.bv[k];  // (i,j)->k
      // IdxSet1.clear();
      // IdxSet2.clear();
      // IdxSet3.clear();
      C = SC.bv[k];
      for (int k1=A10hat.bro[i-1]; k1<A10hat.bro[i]; ++k1) {
        int j1 = A10hat.bci[k1];  // The j1-th column band of A10hat
        A = A10hat.block(i,j1);
        int i1 = j1;  // The i1-th row band of A01
        B = A01.block(i1,j);
        // if (A.columnSize()==B.rowSize()) {
          C = C + (-1)*(A*B);
        // }
      }
      SC.bv[k] = C;
    }
  }

  return (0);  // If successful
}


// JL20161016: TO BE REVISED FOR EFFECIENCY: online?
// Assumption: All matrices and vectors are compatible
// BDSchur matrix-vector multiplication to be used in iterative solvers
/*
friend int BDSchurMatVec(Vector &y, const Vector &x,
                         BlockDiagMatrix &A00inv,
                         const SparseBlockMatrix &A01,
                         const SparseBlockMatrix &A10,
                         const SparseBlockMatrix &A11)
{
  y = A11*x - A10 * (A00inv * (A01*x));
  return (0);  // If successful
}
*/
