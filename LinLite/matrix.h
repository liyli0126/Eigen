// matrix.h
// Classes for (double) full, sparse, and sparse block matrices
// James Liu, Rachel Cali, Graham Harper, ColoState; 2007/01--2018/02

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include "vector.h"

  // Full small size matrices
  class FullMatrix {
  private:
    int m, n;
    double *p;
  public:
    FullMatrix() {m=0; n=0; p=0;}  // Default empty constructor
    FullMatrix(int numRows, int numCols);  // Constructing a zero matrix
    FullMatrix(const FullMatrix &A);  // Copy constructor
    FullMatrix &operator=(const FullMatrix &A);  // Copy assignment
    ~FullMatrix() {delete[] p; m=0; n=0; p=0;}  // Destructor
    int rowSize() const {return m;}
    int columnSize() const {return n;}
    int colSize() const {return n;}  // JL20160826: TO BE REPLACED BY columnSize()
    double &operator()(int i, int j) const;  // 1<=i<=m, 1<=j<=n
    void resize(int numRows, int numCols);  // resize to a zero matrix
    void zeroutRow(int rowPos);
    void zeroutColumn(int colPos);
    void zeroutCol(int colPos);  // JL20160826: TO BE REPLACED BY zeroutColumn()
    void save2file(char *filename) const;
  friend std::ostream &operator<<(std::ostream &strm, const FullMatrix &A);
  friend FullMatrix submatrix(const FullMatrix &A, int bgnRow, int bgnCol, 
                              int numRows, int numCols);
  friend FullMatrix transpose(const FullMatrix &A);
  friend FullMatrix inverseGEPP(const FullMatrix &A);
  friend FullMatrix operator+(const FullMatrix &A, const FullMatrix &B);
  friend FullMatrix operator*(const FullMatrix &A, const FullMatrix &B);
  friend FullMatrix tensorProduct(const FullMatrix &A, const FullMatrix &B);
  friend FullMatrix kron(const FullMatrix &A, const FullMatrix &B);
  friend FullMatrix operator*(double a, const FullMatrix &A);
  friend Vector rowVector(const FullMatrix &A, int i);
  friend Vector columnVector(const FullMatrix &A, int j);
  friend Vector operator*(const FullMatrix &A, const Vector &v);
  friend double colonProduct(const FullMatrix &A, const FullMatrix &B);
  };

  // Diagonal matrices
  // mainly for diagonal preconditioning when solving linear systems
  class DiagMatrix {
  private: 
    int n;
    double *p;
  public: 
    DiagMatrix() {n=0; p=0;}
    DiagMatrix(int order);
    DiagMatrix(const DiagMatrix &A);
    DiagMatrix &operator=(const DiagMatrix &A);
    ~DiagMatrix() {delete[] p;}
    int size() const {return n;}
    double entry(int i) const {return p[i-1];}
    void setEntry(int i, double a) {p[i-1]=a;}
    void addEntry(int i, double a) {p[i-1]+=a;}
    double det() const;
    void resize(int order);  // resize to the zero matrix
    void save2file(char *filename) const;
  };

  // Sparse matrices in the Compressed Row Storage (CRS) format
  class SparseMatrix {
  private: 
    int m, n;  // Matrix row, column sizes
    int nnz;  // Total number of nonzero entries
    int *rowoff;  // row off(set)>=0, length=m+1
    int *colind;  // col(umn) ind(ex)>=1, length=nnz
    double *val;  // val(ue), length=nnz
  public: 
    SparseMatrix() {m=0; n=0; nnz=0; rowoff=0; colind=0; val=0;}
    SparseMatrix(int numRows, int numCols, int numNonZeros);
    SparseMatrix(const SparseMatrix &A);  // Copy constructor
    SparseMatrix &operator=(const SparseMatrix &A);  // Copy assignment
    ~SparseMatrix() {delete[] rowoff,colind,val;}  // Destructor
    int rowSize() const {return m;}
    int columnSize() const {return n;}
    int numberNonzeros() const {return nnz;}
    double entry(int i, int j) const;
    void setEntry(int i, int j, double a);
    void addEntry(int i, int j, double a);
    void unitfyRow(int i);
    void zeroutRow(int i);
    void resize(int numRows, int numCols, int *numEntrsInRow, int **colIdx);
    void save2file(char *filename) const;
    int colSize() const {return n;}  // JL20160826: TO BE REPLACED BY columnSize()
    int numNonZeros() const {return nnz;}  // JL20160826: TO BE REPLACED BY numberNonzeros()
  friend SparseMatrix operator+(const SparseMatrix &A, const SparseMatrix &B);
  friend SparseMatrix operator*(double a, const SparseMatrix &A);
  friend Vector operator*(const SparseMatrix &A, const Vector &v);
  friend DiagMatrix diagonal(const SparseMatrix &A);
  };

  // Block diagonal (square) matrices
  class BlockDiagMatrix {
  private: 
    int n;  // size of the (square) matrix
    int nb;  // number of blocks
    int *db;  // dimension>=1 of each block
    int *bb;  // beginning position>=0 of each block
    FullMatrix *bv;  // each block is a square matrix
  public: 
    BlockDiagMatrix() {n=0; nb=0; db=0; bb=0; bv=0;}
    BlockDiagMatrix(const BlockDiagMatrix &A);
    BlockDiagMatrix &operator=(const BlockDiagMatrix &A);
    ~BlockDiagMatrix() {delete[] db,bb,bv;}
    int size() const {return n;}
    int numberBlocks() const {return nb;}
    int dimBlock(int ib) const {return db[ib-1];}
    int beginPositionBlock(int ib) const {return bb[ib-1];}
    FullMatrix block(int ib) const {return bv[ib-1];}
    void setBlock(int ib, const FullMatrix &A) {bv[ib-1]=A;}
    void addBlock(int ib, const FullMatrix &A) {bv[ib-1]=bv[ib-1]+A;}
    void resize(int numBlks, int *dimBlk);
    void save2file(char *filename) const;
    int numBlocks() const {return nb;}  // JL20160826: TO BE REPLACED BY numberBlocks()
    int bgnPosBlock(int ib) const {return bb[ib-1];}  // JL20160826: TO BE REPLACED BY beginPositionBlock()
  friend Vector operator*(const BlockDiagMatrix &D, const Vector &v);
  friend BlockDiagMatrix CholeskyBlockDiagSPD(const BlockDiagMatrix &A);
  friend BlockDiagMatrix inverseGEPP(const BlockDiagMatrix &A);
  };

  // Sparse block matrices in the Block Compressed Row Storage (BCRS) format
  // Block structure convenient for DGFEMs, WGFEMs, etc.
  // Assuming a tensor structure of the row & column bands
  class SparseBlockMatrix {
  private: 
    int m;  // number of rows of the matrix
    int n;  // number of columns of the matrix
    int nnz;  // total number of all entries in all nonzero blocks
    int nrb;  // number of row bands
    int ncb;  // number of column bands
    int nnzb;  // number of nonzero blocks
    int *drb;  // dim of row band (length=nrb)
    int *dcb;  // dim of column band (length=ncb)
    int *brb;  // (>=1) beginning position (in the whole matrix) of row band
    int *bcb;  // (>=1) beginning position (in the whole matrix) of column band
    int *bro;  // (>=0) (blockwise) row-band offset (length=nrb+1)
    int *bci;  // (>=1) (blockwise) column-band index (length=nnzb)
    FullMatrix *bv;  // length=nnzb, each block (value) is a full matrix
    // int location(int ib, int jb) const;  // >=-1: location for block (ib,jb)
 public: 
    SparseBlockMatrix();  // Default empty constructor
    SparseBlockMatrix(const SparseBlockMatrix &A);  // Copy constructor
    SparseBlockMatrix &operator=(const SparseBlockMatrix &A); // Copy assignment
    ~SparseBlockMatrix();
    int numberRowBands() const {return nrb;}  // number of row bands
    int numberColumnBands() const {return ncb;}  // number of column bands
    int numNonzeroBlocks() const {return nnzb;}  // number nonzero blocks
    int rowSize() const {return m;}
    int columnSize() const {return n;}
    int numberNonzeros() const {return nnz;}  // number of nonzero entries
    FullMatrix block(int ib, int jb) const;
    void setBlock(int ib, int jb, const FullMatrix &A);
    void addBlock(int ib, int jb, const FullMatrix &A);
    void zeroutRow(int ib, int iib);  // zero out a single row
    void zeroutColumn(int jb, int jjb);  // zero out a single column
    void unitifyRow(int ib, int iib);  // unitify a single row
    void unitifyColumn(int jb, int jjb);  // unitify a single column
    void resize(int numRowBands, int numColumnBands,
                int *dimRowBand, int *dimColumnBand,
                int *numBlocksInRowBand, int **columnBandIndex);
    void save2file(char *filename) const;
    int colSize() const {return n;}  // JL20160826: TO BE REPLACED BY columnSize()
    int numRowBlocks() const {return nrb;}  // JL20160826: TO BE REPLACED BY numberRowBands()
    int numColBlocks() const {return ncb;}  // JL20160826: TO BE REPLACED BY numberColumnBands()
    int numNonZeros() const {return nnz;}  // JL20160826: TO BE REPLACED BY numberNonzeros()
    void zeroutCol(int jb, int jjb);  // JL20160826: TO BE REPLACED BY zeroutColumn()
    void unitifyCol(int jb, int jjb);  // JL20160826: TO BE REPLACED BY unitifyColumn()
  friend std::ostream &operator<<(std::ostream &strm, const SparseBlockMatrix &A);
  friend SparseBlockMatrix operator*(double a, const SparseBlockMatrix &A);
  friend SparseBlockMatrix operator+(const SparseBlockMatrix &A,
                                     const SparseBlockMatrix &B);
  friend Vector operator*(const SparseBlockMatrix &A, const Vector &v);
  friend BlockDiagMatrix diagonal(const SparseBlockMatrix &A);
  friend SparseBlockMatrix operator*(const BlockDiagMatrix &D,
                                     const SparseBlockMatrix &A);
  friend SparseBlockMatrix operator*(const SparseBlockMatrix &A,
                                     const BlockDiagMatrix &D);
  //
  // "SparseBlockMatrix * SparseBlockMatrix" nontrivial, to be implemented later
  /*
  friend SparseBlockMatrix operator*(const SparseBlockMatrix &A,
                                     const SparseBlockMatrix &B);
  */
  //
  // JL20161109: NOT WORKING YET: SC lost symmetry
  // BDSchur: Block Diagonal Schur complement
  friend int BDSchur(SparseBlockMatrix &SC,
                     BlockDiagMatrix &A00inv, SparseBlockMatrix &A10hat,
                     const BlockDiagMatrix &A00, const SparseBlockMatrix &A01,
                     const SparseBlockMatrix &A10, const SparseBlockMatrix &A11);
  };

#endif  // MATRIX_H
