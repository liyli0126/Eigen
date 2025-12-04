// mat3.h
// A class for 3x3 double matrices
// James Liu, ColoState; 2007/01--2018/01

#ifndef Mat3_H
#define Mat3_H

#include <iostream>
#include "matrix.h"
#include "PtVec3d.h"

class Mat3 : public FullMatrix {
private:
  double entry[3][3];
public:
  Mat3();  // Default constructor
  //
  // A constructor based on aligning these vectors into a matrix
  // 'r' or 'c': row/column-wise alignment
  Mat3(char rc, const PtVec3d &v1, const PtVec3d &v2, const PtVec3d &v3);
  //
  // A constructor based on tensor product of two PtVec3d vectors
  // 'c' or 'r': the 1st vector is used as a column or row vector
  Mat3(char cr, const PtVec3d &v1, const PtVec3d &v2);
  //
  Mat3(const Mat3 &A);  // Copy constructor
  void reset();  // reset to a zero 3-by-3 matrix
  Mat3 &operator=(const Mat3 &A);  // Copy assignment
  double &operator()(int i, int j);  // element (i,j)
  PtVec3d row(int i) const;  // getting row i
  PtVec3d column(int j) const; // getting column j
  double det() const;
  Mat3 transpose() const;
  Mat3 inverse() const;
friend std::ostream &operator<<(std::ostream&strm, const Mat3&);
friend Mat3 operator*(double a, const Mat3 &A);
friend Mat3 operator*(const Mat3 &A, const Mat3 &B);
friend Mat3 operator+(const Mat3 &A, const Mat3 &B);
friend Mat3 operator-(const Mat3 &A, const Mat3 &B);
friend PtVec3d operator*(const Mat3 &A, const PtVec3d &v);
friend double colonProduct(const Mat3 &A, const Mat3 &B);
};

#endif  // Mat3_H
// mat3.h
