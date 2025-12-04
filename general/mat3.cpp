// mat3.cpp
// A class for 3x3 double matrices
// James Liu, ColoState; 2007/01--2018/01

#include <iostream>

#include "mat3.h"
#include "PtVec3d.h"


// Mat3: Default constructor
// all entries set to 0

Mat3::Mat3() 
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      entry[i][j] = 0.0;
}


// Mat3: A constructor based on aligning these vectors into a matrix
// 'r' or 'c': row/column-wise alignment
// 'r' or 'c': row/column-wise constructing from 3 vectors

Mat3::Mat3(char rc, const PtVec3d &v1, const PtVec3d &v2, const PtVec3d &v3)
{
  if (rc=='r') {
    entry[0][0] = v1.xCrd();
    entry[0][1] = v1.yCrd();
    entry[0][2] = v1.zCrd();

    entry[1][0] = v2.xCrd();
    entry[1][1] = v2.yCrd();
    entry[1][2] = v2.zCrd();

    entry[2][0] = v3.xCrd();
    entry[2][1] = v3.yCrd();
    entry[2][2] = v3.zCrd();
  }

  if (rc=='c') {
    entry[0][0] = v1.xCrd();
    entry[1][0] = v1.yCrd();
    entry[2][0] = v1.zCrd();

    entry[0][1] = v2.xCrd();
    entry[1][1] = v2.yCrd();
    entry[2][1] = v2.zCrd();

    entry[0][2] = v3.xCrd();
    entry[1][2] = v3.yCrd();
    entry[2][2] = v3.zCrd();
  }
}


// Mat3: A constructor based on tensor product of two PtVec3d vectors
// 'c' or 'r': the 1st vector is used as a column or row vector

Mat3::Mat3(char cr, const PtVec3d &v1, const PtVec3d &v2)
{
  double a[3], b[3];
  a[0] = v1.xCrd();  a[1] = v1.yCrd();  a[2] = v1.zCrd();
  b[0] = v2.xCrd();  b[1] = v2.yCrd();  b[2] = v2.zCrd();

  if (cr=='c') {
    for (int i=0; i<=2; ++i)
      for (int j=0; j<=2; ++j)
        entry[i][j] = a[i] * b[j];
  }
  if (cr=='r') {
    for (int i=0; i<=2; ++i)
      for (int j=0; j<=2; ++j)
        entry[i][j] = a[j] * b[i];
  }
}


// Mat3: Copy constructor

Mat3::Mat3(const Mat3 &A) 
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      entry[i][j] = A.entry[i][j];
}


// Mat3: Reset

void Mat3::reset()
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      entry[i][j] = 0.0;
  return;
}


// Mat3: Copy assignment

Mat3 &Mat3::operator=(const Mat3 &A) 
{
  if (this != &A) {
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        entry[i][j] = A.entry[i][j];
  }
  return *this;
}


// Mat3: Getting and setting element (i,j)

double &Mat3::operator()(int i, int j)
{
  return entry[i-1][j-1];
}


// Mat3: Getting row i

PtVec3d Mat3::row(int i) const
{
  return PtVec3d(entry[i-1][0], entry[i-1][1], entry[i-1][2]);
}


// Mat3: Getting column j

PtVec3d Mat3::column(int j) const
{
  return PtVec3d(entry[0][j-1], entry[1][j-1], entry[2][j-1]);
}


// Mat3: determinant

double Mat3::det() const 
{
  double x1 = entry[0][0];
  double y1 = entry[0][1];
  double z1 = entry[0][2];

  double x2 = entry[1][0];
  double y2 = entry[1][1];
  double z2 = entry[1][2];

  double x3 = entry[2][0];
  double y3 = entry[2][1];
  double z3 = entry[2][2];

  double d = x1*y2*z3 + x2*y3*z1 + x3*y1*z2;
  d = d - x3*y2*z1 - x1*y3*z2 - x2*y1*z3;

  return d;
}


// Mat3: transpose 

Mat3 Mat3::transpose() const 
{
  Mat3 B;

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      B.entry[i][j] = entry[j][i];

  return B;
}


// Mat3: inverse matrix

Mat3 Mat3::inverse() const 
{
  Mat3 B;

  B.entry[0][0] = entry[1][1]*entry[2][2]-entry[1][2]*entry[2][1];
  B.entry[0][1] = -(entry[1][0]*entry[2][2]-entry[1][2]*entry[2][0]);
  B.entry[0][2] = entry[1][0]*entry[2][1]-entry[1][1]*entry[2][0];

  B.entry[1][0] = -(entry[0][1]*entry[2][2]-entry[0][2]*entry[2][1]);
  B.entry[1][1] = entry[0][0]*entry[2][2]-entry[0][2]*entry[2][0];
  B.entry[1][2] = -(entry[0][0]*entry[2][1]-entry[0][1]*entry[2][0]);

  B.entry[2][0] = entry[0][1]*entry[1][2]-entry[0][2]*entry[1][1];
  B.entry[2][1] = -(entry[0][0]*entry[1][2]-entry[0][2]*entry[1][0]);
  B.entry[2][2] = entry[0][0]*entry[1][1]-entry[0][1]*entry[1][0];

  return (1.0/det())*B;
}


// Output the 3x3 matrix

std::ostream &operator<<(std::ostream&strm, const Mat3&A)
{
  std::cout << "\n";
  
  std::cout << "row1 = ";
  std::cout << A.entry[0][0] << ", ";
  std::cout << A.entry[0][1] << ", ";
  std::cout << A.entry[0][2] << "\n";
  
  std::cout << "row2 = ";
  std::cout << A.entry[1][0] << ", ";
  std::cout << A.entry[1][1] << ", ";
  std::cout << A.entry[1][2] << "\n";
  
  std::cout << "row3 = ";
  std::cout << A.entry[2][0] << ", ";
  std::cout << A.entry[2][1] << ", ";
  std::cout << A.entry[2][2] << "\n";
  
  return strm;
}


// Scalar-Mat3 multiplication

Mat3 operator*(double a, const Mat3 &A) 
{
  Mat3 B;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      B.entry[i][j] = a*A.entry[j][i];
  return B;
}


// Mat3-Mat3 multiplication

Mat3 operator*(const Mat3 &A, const Mat3 &B) 
{
  Mat3 C;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      C.entry[i][j] = A.entry[i][0]*B.entry[0][j]
                    + A.entry[i][1]*B.entry[1][j]
                    + A.entry[i][2]*B.entry[2][j];
  return C;
}


// Mat3-Mat3 addition

Mat3 operator+(const Mat3 &A, const Mat3 &B)
{
  Mat3 C;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      C.entry[i][j] = A.entry[i][j] + B.entry[i][j];
  return C;
}


// Mat3-Mat3 substraction

Mat3 operator-(const Mat3 &A, const Mat3 &B)
{
  Mat3 C;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      C.entry[i][j] = A.entry[i][j] - B.entry[i][j];
  return C;
}


// Mat3-OtVec3d (matrix-vector) multiplication

PtVec3d operator*(const Mat3 &A, const PtVec3d &v) 
{
  double x = v.xCrd();
  double y = v.yCrd();
  double z = v.zCrd();
  double s1 = A.entry[0][0]*x + A.entry[0][1]*y + A.entry[0][2]*z;
  double s2 = A.entry[1][0]*x + A.entry[1][1]*y + A.entry[1][2]*z;
  double s3 = A.entry[2][0]*x + A.entry[2][1]*y + A.entry[2][2]*z;
  return PtVec3d(s1, s2, s3);
}


// The colon product of two 3x3 matrices

double colonProduct(const Mat3 &A, const Mat3 &B)
{
  double dp = 0.0;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      dp += A.entry[i][j] * B.entry[i][j];
  return dp;
}

// mat3.cpp
