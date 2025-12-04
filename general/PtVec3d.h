// PtVec3d.h
// A class for points and vectors in 3-dim
// James Liu, ColoState; 2007/01--2018/01

#ifndef PtVec3d_H
#define PtVec3d_H

#include <cmath>
#include <iostream>

class PtVec3d {
private:
  double x, y, z;
public:
  // Default constructor:
  PtVec3d(double a=0, double b=0, double c=0) {x=a; y=b; z=c;}
  PtVec3d(const PtVec3d &P) {x=P.x; y=P.y; z=P.z;}  // Copy constructor
  PtVec3d &operator=(const PtVec3d &P);  // Copy assignment
  double xCrd() const {return x;}
  double yCrd() const {return y;}
  double zCrd() const {return z;}
  double l2norm() const {return sqrt(x*x+y*y+z*z);}
  PtVec3d operator-() const {return PtVec3d(-x,-y,-z);}
friend std::ostream &operator<<(std::ostream &strm, const PtVec3d &P);
friend double dist(const PtVec3d &P1, const PtVec3d &P2);
friend PtVec3d operator*(double a, const PtVec3d &v);
friend PtVec3d operator*(const PtVec3d &v, double a);
friend PtVec3d operator+(const PtVec3d &v1, const PtVec3d &v2);
friend PtVec3d operator-(const PtVec3d &v1, const PtVec3d &v2);
friend PtVec3d crossProduct(const PtVec3d &v1, const PtVec3d &v2);
friend double dotProduct(const PtVec3d &v1, const PtVec3d &v2);
// For tensor product of 2 PtVec3d vectors resulting in a 3x3 matrix,
// see mat3.h, mat3.cpp
};

#endif  // PtVec3d_H
// PtVec3d.h
