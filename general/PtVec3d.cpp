// PtVec3d.cpp
// A class for points and vectors in 3-dim
// James Liu, ColoState; 2007/01--2018/01

#include <cmath>
#include <iostream>
#include "PtVec3d.h"


// PtVec3d: Copy assignment

PtVec3d &PtVec3d::operator=(const PtVec3d &P) 
{
  if (this != &P) {
     x = P.x;  y = P.y;  z = P.z;
  }
  return *this;
}


// PtVec3d: output stream

std::ostream &operator<<(std::ostream &strm, const PtVec3d &P)
{
  std::cout << "(x,y,z)= " << P.x << ", " << P.y << ", " << P.z;
  return strm;
}


double dist(const PtVec3d &P1, const PtVec3d &P2)
{
  double x = P1.x - P2.x;
  double y = P1.y - P2.y;
  double z = P1.z - P2.z;
  return sqrt(x*x+y*y+z*z);
}


PtVec3d operator*(double a, const PtVec3d &v) 
{
  return PtVec3d(a*v.x, a*v.y, a*v.z);
}


PtVec3d operator*(const PtVec3d &v, double a)
{
  return PtVec3d(a*v.x, a*v.y, a*v.z);
}


PtVec3d operator+(const PtVec3d &v1, const PtVec3d &v2)
{
  return PtVec3d(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}


PtVec3d operator-(const PtVec3d &v1, const PtVec3d &v2) 
{
  return PtVec3d(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}


PtVec3d crossProduct(const PtVec3d &v1, const PtVec3d &v2)
{
  double x = v1.y * v2.z - v2.y * v1.z;
  double y = v1.z * v2.x - v2.z * v1.x;
  double z = v1.x * v2.y - v2.x * v1.y;
  return PtVec3d(x, y, z);
}


double dotProduct(const PtVec3d &v1, const PtVec3d &v2)
{
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}


/*
bool operator<=(const PtVec3d &P1, const PtVec3d &P2) 
{
  bool b = false;
  if (P1.x<=P2.x && P1.y<=P2.y && P1.z<=P2.z)  b = true;
  return b;
}
*/

// PtVec3d.cpp
