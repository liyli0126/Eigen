// EqnBC_Poisson3d_Ex12coscoscos.cpp
// James Liu, ColoState; 2014/07--2016/05

#include <cmath>
#include "mat3.h"
#include "PtVec3d.h"

#ifndef PI
#define PI 3.141592653589793
#endif


// Full diffusion/permeability tensor (3-by-3 matrix)

Mat3 fxnMatK(PtVec3d pt)
{
  Mat3 MatK;
  MatK(1,1) = 1;  MatK(2,2) = 1;  MatK(3,3) = 1;
  return MatK;
}


// The right-hand side function f in the PDE

double fxnf(PtVec3d pt)
{
  double x = pt.xCrd();
  double y = pt.yCrd();
  double z = pt.zCrd();
  double f = 3*PI*PI*cos(PI*x)*cos(PI*y)*cos(PI*z);
  return f;
}


// The exact pressure p (assuming known)

double fxnp(PtVec3d pt)
{
  double x = pt.xCrd();
  double y = pt.yCrd();
  double z = pt.zCrd();
  double p = cos(PI*x)*cos(PI*y)*cos(PI*z);
  return p;
}


// Dirichlet boundary conditions

double fxnpD(PtVec3d pt)
{
  double pD;
  pD = fxnp(pt);
  return pD;
}


// Neumann boundary conditions

double fxnuN(PtVec3d pt)
{
  double uN;
  // JL20141226: TO BE REVISED
  uN = 0.0;
  return uN;
}


// The known exact gradient

PtVec3d fxnpg(PtVec3d pt)
{
  double x = pt.xCrd();
  double y = pt.yCrd();
  double z = pt.zCrd();
  double g1 = -PI*sin(PI*x)*cos(PI*y)*cos(PI*z);
  double g2 = -PI*cos(PI*x)*sin(PI*y)*cos(PI*z);
  double g3 = -PI*cos(PI*x)*cos(PI*y)*sin(PI*z);
  return PtVec3d(g1,g2,g3);
}


// The known exact (Darcy) velocity

PtVec3d fxnu(PtVec3d pt)
{
  double x = pt.xCrd();
  double y = pt.yCrd();
  double z = pt.zCrd();
  double g1 = -PI*sin(PI*x)*cos(PI*y)*cos(PI*z);
  double g2 = -PI*cos(PI*x)*sin(PI*y)*cos(PI*z);
  double g3 = -PI*cos(PI*x)*cos(PI*y)*sin(PI*z);
  PtVec3d pg = PtVec3d(g1,g2,g3);
  Mat3 MatK = fxnMatK(pt);
  PtVec3d vel = (-1)*(MatK*pg);
  return vel;
}