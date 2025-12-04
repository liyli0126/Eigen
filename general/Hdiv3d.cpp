// JL20171102: TO BE FINISHED/REVISED
// Hdiv3d.cpp
// Broadly defined "Hdiv" bases 
// The Jacobian determinant value should be positive (tetra./hexa. orientation)
// James Liu, Graham Harper, ColoState; 2014/07--2017/11

#include <cmath>

#include "LinSys.h"
#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "mat3.h"
#include "PtVec3d.h"


////////////////////////////////////////////////////////////////////////////////
// For bricks
////////////////////////////////////////////////////////////////////////////////

// Brick: (3-dim) brick with coordinates (x0,y0,z0), (x1,y1,z1)
// RT[0] normalized basis functions:
//   [1;0;0], [0;1;0], [0;0;1], [X;0;0], [0;Y;0], [0;0;Z]
// where X=x-xc, Y=y-yc, Z=z-zc, (xc,yc,zc) element center


// Brick RT[0]: The Gram matrix and its inverse (both diagonal matrcies)

int Hdiv_BrickRT0_GramMat_NmlzBas(double x0, double y0, double z0,
  double x1, double y1, double z1, double GM[6], double GMI[6])
{
  double Deltax, Deltay, Deltaz;
  double Deltax2, Deltay2, Deltaz2;
  double vol, vol1;

  Deltax = x1 - x0;
  Deltay = y1 - y0;
  Deltaz = z1 - z0;

  Deltax2 = Deltax * Deltax;
  Deltay2 = Deltay * Deltay;
  Deltaz2 = Deltaz * Deltaz;
  
  vol = Deltax * Deltay * Deltaz;
  vol1 = 1.0/vol;
  
  // The diagonal Gram matrix (GM)
  GM[0] = vol;  GM[3] = vol*(Deltax2/12);
  GM[1] = vol;  GM[4] = vol*(Deltay2/12);
  GM[2] = vol;  GM[5] = vol*(Deltaz2/12);

  // The inverse of the diagonal matrix (GMI)
  GMI[0] = vol1;  GMI[3] = vol1*(12/Deltax2);
  GMI[1] = vol1;  GMI[4] = vol1*(12/Deltay2);
  GMI[2] = vol1;  GMI[5] = vol1*(12/Deltaz2);

  return(0);
}


// Brick RT[0]: Total normal fluxes on all 6 faces for all 6 nmlz.bas.fxns.
// row for faces: 0,1: x-left,right; 2,3: y-back,front; 4,5: z-bottom,fact
// column for basis functions: 3 std.unit.vectors, then those for X,Y,Z resp.

int Hdiv_BrickRT0_NmlFlux_NmlzBas(double x0, double y0, double z0,
  double x1, double y1, double z1, double NmlFlux[][6])
{
  double Deltax = x1 - x0;
  double Deltay = y1 - y0;
  double Deltaz = z1 - z0;
  double vol05 = 0.5*Deltax * Deltay * Deltaz;

  for (int i=0; i<6; ++i)
    for (int j=0; j<6; ++j)
      NmlFlux[i][j] = 0;

  NmlFlux[0][0] = -Deltay*Deltaz;
  NmlFlux[1][0] =  Deltay*Deltaz;
  
  NmlFlux[2][1] = -Deltaz*Deltax;
  NmlFlux[3][1] =  Deltaz*Deltax;

  NmlFlux[4][2] = -Deltax*Deltay;
  NmlFlux[5][2] =  Deltax*Deltay;

  NmlFlux[0][3] = vol05;
  NmlFlux[1][3] = vol05;

  NmlFlux[2][4] = vol05;
  NmlFlux[3][4] = vol05;

  NmlFlux[4][5] = vol05;
  NmlFlux[5][5] = vol05;

  return(0);  // if successful
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Brick RT[0]: The Gram matrix for normalized basis and
// full (3x3) diffusion/permeability matrix/tensor

FullMatrix Hdiv_BrickRT0_GramMatK_NmlzBas(const Brick &brick, Mat3 &MatK)
{
  int i, j;
  double Deltax, Deltay, Deltaz, vol;
  PtVec3d P0, P1;
  FullMatrix GramMatK(6,6);

  brick.getVertices(P0, P1);
  Deltax = P1.xCrd() - P0.xCrd();
  Deltay = P1.yCrd() - P0.yCrd();
  Deltaz = P1.zCrd() - P0.zCrd();
  vol = Deltax * Deltay * Deltaz;

  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GramMatK(i,j) = MatK(i,j) * vol;

  GramMatK(4,4) = MatK(1,1) * (Deltax*Deltax/12) * vol;
  GramMatK(5,5) = MatK(2,2) * (Deltay*Deltay/12) * vol;
  GramMatK(6,6) = MatK(3,3) * (Deltaz*Deltaz/12) * vol;

  return GramMatK;
}


// Brick RT[0]: Normal fluxes of the normalized basis functions
// 6x6: 6 rows for side faces, 6 columns for nmlz.bas.fxns.

FullMatrix Hdiv_BrickRT0_NmlFlux_NmlzBas(const Brick &brick)
{
  double Deltax, Deltay, Deltaz, vol05;
  PtVec3d P0, P1;
  FullMatrix NmlFluxRT0(6,6);

  brick.getVertices(P0, P1);
  Deltax = P1.xCrd() - P0.xCrd();
  Deltay = P1.yCrd() - P0.yCrd();
  Deltaz = P1.zCrd() - P0.zCrd();
  vol05 = 0.5 * Deltax * Deltay * Deltaz;

  NmlFluxRT0(1,1) = -Deltay*Deltaz;
  NmlFluxRT0(2,1) =  Deltay*Deltaz;

  NmlFluxRT0(3,2) = -Deltaz*Deltax;
  NmlFluxRT0(4,2) =  Deltaz*Deltax;

  NmlFluxRT0(5,3) = -Deltax*Deltay;
  NmlFluxRT0(6,3) =  Deltax*Deltay;

  NmlFluxRT0(1,4) = vol05;
  NmlFluxRT0(2,4) = vol05;

  NmlFluxRT0(3,5) = vol05;
  NmlFluxRT0(4,5) = vol05;

  NmlFluxRT0(5,6) = vol05;
  NmlFluxRT0(6,6) = vol05;

  return NmlFluxRT0;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Brick P1^3: Polynomial degree 1 vector 3-dim
// The Gram matrix for normalized basis functions (X,Y,Z)

int Hdiv_BrickP1v3_GramMat_NmlzBas(DiagMatrix &GM, const Brick &brick)
{
  int i;
  double Deltax, Deltay, Deltaz, x0, x1, y0, y1, z0, z1;
  double SX2, SY2, SZ2, vol;
  PtVec3d P0, P1;

  brick.getVertices(P0, P1);

  x0 = P0.xCrd();  y0 = P0.yCrd();  z0 = P0.zCrd();
  x1 = P1.xCrd();  y1 = P1.yCrd();  z1 = P1.zCrd();
  
  Deltax = x1 - x0;
  Deltay = y1 - y0;
  Deltaz = z1 - z0;

  vol = Deltax * Deltay * Deltaz;

  SX2 = (1.0/12) * (Deltax * Deltax) * vol;
  SY2 = (1.0/12) * (Deltay * Deltay) * vol;
  SZ2 = (1.0/12) * (Deltaz * Deltaz) * vol;

  GM.resize(12);

  GM.setEntry(1,vol);
  GM.setEntry(2,vol);
  GM.setEntry(3,vol);

  for (i=1; i<=3; ++i) {
    GM.setEntry(3*i+1, SX2);
    GM.setEntry(3*i+2, SY2);
    GM.setEntry(3*i+3, SZ2);
  }

  return(0);  // if successful
}


// Brick P1^3: Polynomial degree 1 vector 3-dim
// The Gram matrix for full 3-by-3 permeability tensor MatK
// and normalized basis functions (X,Y,Z)

int Hdiv_BrickP1v3_GramMatK_NmlzBas(FullMatrix &GMK, const Brick &brick,
                                    Mat3 &MatK)
{
  int i, j;
  double Deltax, Deltay, Deltaz, x0, x1, y0, y1, z0, z1;
  double SX2, SY2, SZ2, vol;
  PtVec3d P0, P1;
  FullMatrix MatA(3,3), MatB(9,9), MatKK(3,3);
  
  brick.getVertices(P0, P1);
  
  x0 = P0.xCrd();  y0 = P0.yCrd();  z0 = P0.zCrd();
  x1 = P1.xCrd();  y1 = P1.yCrd();  z1 = P1.zCrd();
  
  Deltax = x1 - x0;
  Deltay = y1 - y0;
  Deltaz = z1 - z0;
  
  vol = Deltax * Deltay * Deltaz;
  
  SX2 = (1.0/12) * (Deltax * Deltax) * vol;
  SY2 = (1.0/12) * (Deltay * Deltay) * vol;
  SZ2 = (1.0/12) * (Deltaz * Deltaz) * vol;
  
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      MatKK(i,j) = MatK(i,j);
  MatA(1,1) = SX2;  MatA(2,2) = SY2;  MatA(3,3) = SZ2;
  MatB = tensorProduct(MatKK, MatA);

  GMK.resize(12,12);
  //
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GMK(i,j) = MatK(i,j) * vol;
  //
  for (i=1; i<=9; ++i)
    for (j=1; j<=9; ++j)
      GMK(i+3,j+3) = MatB(i,j);

  return(0);  // if successful
}


////////////////////////////////////////////////////////////////////////////////
// For tetrahedra
////////////////////////////////////////////////////////////////////////////////

// Functions for tetrahedra

// Tetra: (3-dim) tetrahedron with coordinates (x[i],y[i],z[i]), i=0,1,2,3
// RT0 normalized basis functions:
//   [1;0;0], [0;1;0], [0;0;1], [X;Y;Z]
// where X=x-xc, Y=y-yc, Z=z-zc, (xc,yc,zc) element center


// The Gram matrix and its inverse for a tetrahedron
/*
int Hdiv_TetraRT0_GramMat_NmlzBas(double x[], double y[], double z[],
                                  double GM[4], double GMI[4])
{
  // double det0, det1, det2, det3;
  double a1, a2, a3, b1, b2, b3, c1, c2, c3, det;
  double vol, vol1, w, w1;
  double x10, x20, x30, x21, x31, x32;
  double y10, y20, y30, y21, y31, y32;
  double z10, z20, z30, z21, z31, z32;
  
  a1 = x[1] - x[0];  b1 = y[1] - y[0];  c1 = z[1] - z[0];
  a2 = x[2] - x[0];  b2 = y[2] - y[0];  c2 = z[2] - z[0];
  a3 = x[3] - x[0];  b3 = y[3] - y[0];  c3 = z[3] - z[0];
  det = a1*b2*c3 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 - a3*b2*c1;
  vol = det/6.0;  // 6=3!
  vol1 = 1.0/vol;
  
  // w = 1/(5*4*4) * ( X^T*A*X + Y^T*A*Y + Z^T*A*Z )
  // where A = [3,-1,-1,-1; -1,3,-1,-1; -1,-1,3,-1; -1,-1,-1,3]
  // and X = [x[0];x[1];x[2];x[3]], similarly for Y, Z
  
  x10 = x[1] - x[0];
  x20 = x[2] - x[0];  x21 = x[2] - x[1];
  x30 = x[3] - x[0];  x31 = x[3] - x[1];  x32 = x[3] - x[2];
  
  y10 = y[1] - y[0];
  y20 = y[2] - y[0];  y21 = y[2] - y[1];
  y30 = y[3] - y[0];  y31 = y[3] - y[1];  y32 = y[3] - y[2];
  
  z10 = z[1] - z[0];
  z20 = z[2] - z[0];  z21 = z[2] - z[1];
  z30 = z[3] - z[0];  z31 = z[3] - z[1];  z32 = z[3] - z[2];
  
  // w = \int\int\int_{Te} (X^2+Y^2+Z^2) dxdydz
  w = x10*x10 + x20*x20 + x30*x30 + x21*x21 + x31*x31 + x32*x32;
  w = y10*y10 + y20*y20 + y30*y30 + y21*y21 + y31*y31 + y32*y32 + w;
  w = z10*z10 + z20*z20 + z30*z30 + z21*z21 + z31*z31 + z32*z32 + w;
  w = 0.0125*vol*w;  // 0.0125 = 1/(5*4*4)
  w1 = 1.0/w;
  
  // The diagonal Gram matrix (GM)
  GM[0] = vol;  GM[1] = vol;  GM[2] = vol;  GM[3] = w;
  
  // The inverse of the diagonal matrix (GMI)
  GMI[0] = vol1;  GMI[1] = vol1;  GMI[2] = vol1;  GMI[3] = w1;
  
  return(0);  // if successful
}
*/


// GM: Order 4 diag.mat. for the normalized RT0 basis

int Hdiv_TetraRT0_GramMat_NmlzBas(DiagMatrix &GM, const Tetra &tetra,
                                  const GaussQuad &GQTe)
{
  int i, j;
  double vol, w;
  double x, y, z, xc, yc, zc, X, Y, Z;
  double SX2, SY2, SZ2, SXY, SXZ, SYZ;
  PtVec3d cntr, qp;

  //
  cntr = tetra.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  vol = tetra.volume();
  // vol1 = 1.0/vol;

  // Computing auxiliary integrals using the chosen Gaussian quadrature
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (i=0; i<GQTe.numberQuadraturePoints(); ++i) {
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQTe.numberVertices(); ++j) {
      qp = qp + GQTe.baryCoordinate(i,j) * tetra.vertex(j);
    }
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    SX2 += GQTe.weight(i) * vol * (X*X);
    SY2 += GQTe.weight(i) * vol * (Y*Y);
    SZ2 += GQTe.weight(i) * vol * (Z*Z);
    SXY += GQTe.weight(i) * vol * (X*Y);
    SXZ += GQTe.weight(i) * vol * (X*Z);
    SYZ += GQTe.weight(i) * vol * (Y*Z);
  }
  w = SX2 + SY2 + SZ2;
  // w1 = 1.0/w;

  // For the Gram matrix (GM)
  GM.resize(4);
  GM.setEntry(1, vol);
  GM.setEntry(2, vol);
  GM.setEntry(3, vol);
  GM.setEntry(4, w);
  
  // For the Gram matrix inverse (GMI)
  /*
  GMI.resize(4);
  GMI.setEntry(1, vol1);
  GMI.setEntry(2, vol1);
  GMI.setEntry(3, vol1);
  GMI.setEntry(4, w1);
  */

  return(0);  // if successful
}


// GMK: 4x4 "Gram" matrix with full permeability matrix MatK

int Hdiv_TetraRT0_GramMatK_NmlzBas(FullMatrix &GMK, const Tetra &tetra,
                                   Mat3 &MatK, const GaussQuad &GQTe)
{
  int i, j;
  double vol;
  double x, y, z, xc, yc, zc, X, Y, Z;
  double SX2, SY2, SZ2, SXY, SXZ, SYZ;
  PtVec3d cntr, qp;
  FullMatrix MatS(3,3);
  
  //
  cntr = tetra.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  vol = tetra.volume();
  
  // Computing auxiliary integrals using the chosen Gaussian quadrature
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (i=0; i<GQTe.numberQuadraturePoints(); ++i) {
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQTe.numberVertices(); ++j) {
      qp = qp + GQTe.baryCoordinate(i,j) * tetra.vertex(j);
    }
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    SX2 += GQTe.weight(i) * vol * (X*X);
    SY2 += GQTe.weight(i) * vol * (Y*Y);
    SZ2 += GQTe.weight(i) * vol * (Z*Z);
    SXY += GQTe.weight(i) * vol * (X*Y);
    SXZ += GQTe.weight(i) * vol * (X*Z);
    SYZ += GQTe.weight(i) * vol * (Y*Z);
  }
  
  //
  MatS(1,1) = SX2;  MatS(1,2) = SXY;  MatS(1,3) = SXZ;
  MatS(2,1) = SXY;  MatS(2,2) = SY2;  MatS(2,3) = SYZ;
  MatS(3,1) = SXZ;  MatS(3,2) = SYZ;  MatS(3,3) = SZ2;

  //
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GMK(i,j) = MatK(i,j) * vol;
  //
  GMK(4,4) = 0;
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GMK(4,4) += MatK(i,j) * MatS(i,j);

  return(0);  // if successful
}


// JL20150707: NEEDS DOUBLE-CHECK!!
// NmlFlux: 4x4 matrix, faces as rows, basis functions as columns

int Hdiv_TetraRT0_NmlFlux_NmlzBas(FullMatrix &NmlFlux, const Tetra &tetra)
{
  int i, j;
  double sa[4], sign[4];  // 4 side face areas and related signs (+1/-1)
  PtVec3d EmCntr, FcCntr[4], nml[4], w[4][4];
  Tri3d fc[4];
  NmlFlux.resize(4,4);
  
  // The tetrahedron center
  EmCntr = tetra.center();

  // The 4 faces, their areas, centers/midpoints, & unit outer normal vectors,
  // and nmlz.bas.fxn.values at these midpoints
  for (i=0; i<4; ++i) {
    fc[i] = tetra.face(i);
    sa[i] = fc[i].area();
    FcCntr[i] = fc[i].center();
    w[i][0] = PtVec3d(1,0,0);
    w[i][1] = PtVec3d(0,1,0);
    w[i][2] = PtVec3d(0,0,1);
    w[i][3] = FcCntr[i] - EmCntr;
    nml[i] = tetra.normal(i);
    sign[i] = 1;
    if (dotProduct(FcCntr[i]-EmCntr,nml[i])<0)  sign[i] = -1;
  }

  // Now the outer normal fluxes
  for (i=1; i<=4; ++i)
    for (j=1; j<=4; ++j)
      NmlFlux(i,j) = dotProduct(w[i-1][j-1],nml[i-1]) * sa[i-1];

  return(0);  // if successful
}


// KtransMat: 12x4 matrix, P1v3 as rows, RT0 as columns, not -K
// The K-trans. image of j-th RT0 bas.fxn. = lin. comb. of 12 P_1^3 bas.fxn.

int Hdiv_Tetra01_KtransMat(FullMatrix &KtransMat, Mat3 &MatK)
{
  int i, j;
  KtransMat.resize(12,4);

  //
  for (j=1; j<=3; ++j)
    for (i=1; i<=3; ++i)
      KtransMat(i,j) = MatK(i,j);
  //
  for (i=1; i<=3; ++i) {
    KtransMat(3*i+1,4) = MatK(i,1);
    KtransMat(3*i+2,4) = MatK(i,2);
    KtransMat(3*i+3,4) = MatK(i,3);
  }
  
  return(0);  // if successful
}


// ProjMat: 4X12 matrix, RT0 as rows, P1v3 as columns

int Hdiv_Tetra01_ProjP1v3RT0(FullMatrix &ProjMat, const Tetra &tetra,
                             const GaussQuad &GQTe)
{
  int i, j;
  double X, Y, Z, x, y, z, xc, yc, zc;
  double SX2, SY2, SZ2, SXY, SXZ, SYZ, vol;
  PtVec3d cntr, qp;  // w[0] not used
  DiagMatrix GM(4);
  ProjMat.resize(4,12);

  //
  cntr = tetra.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();

  //
  Hdiv_TetraRT0_GramMat_NmlzBas(GM, tetra, GQTe);
  vol = GM.entry(1);
  
  // Computing auxiliary integrals on the tetrahedral element
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (i=0; i<GQTe.numberQuadraturePoints(); ++i) {
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQTe.numberVertices(); ++j) {
     qp = qp + GQTe.baryCoordinate(i,j) * tetra.vertex(j);
    }
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    SX2 += GQTe.weight(i) * (X*X);
    SY2 += GQTe.weight(i) * (Y*Y);
    SZ2 += GQTe.weight(i) * (Z*Z);
    SXY += GQTe.weight(i) * (X*Y);
    SXZ += GQTe.weight(i) * (X*Z);
    SYZ += GQTe.weight(i) * (Y*Z);
  }
  SX2 *= vol;  SY2 *= vol;  SZ2 *= vol;
  SXY *= vol;  SXZ *= vol;  SYZ *= vol;

  // Computing the ProjMat(4,12) matrix
  ProjMat(1,1) = 1;
  ProjMat(2,2) = 1;
  ProjMat(3,3) = 1;
  //
  ProjMat(4,4) = SX2 / GM.entry(4);
  ProjMat(4,5) = SXY / GM.entry(4);
  ProjMat(4,6) = SXZ / GM.entry(4);
  //
  ProjMat(4,7) = SXY / GM.entry(4);
  ProjMat(4,8) = SY2 / GM.entry(4);
  ProjMat(4,9) = SYZ / GM.entry(4);
  //
  ProjMat(4,10) = SXZ / GM.entry(4);
  ProjMat(4,11) = SYZ / GM.entry(4);
  ProjMat(4,12) = SZ2 / GM.entry(4);

  return(0);  // if successful
}


// GM: Order 15 full matrix for the normalized RT1 basis

int Hdiv_TetraRT1_GramMat_NmlzBas(FullMatrix &GM, const Tetra &tetra,
                                  const GaussQuad &GQTe)
{
  int i, j;
  double vol, x, y, z, xc, yc, zc, X, Y, Z;
  double SX2, SY2, SZ2, SXY, SXZ, SYZ;
  double SX2Y, SX2Z, SXY2, SXZ2, SY2Z, SYZ2;
  double SX3, SY3, SZ3, SXYZ;
  PtVec3d cntr, qp;
  FullMatrix A(12,12), C(3,12), B(12,3), D(3,3);
  GM.resize(15,15);

  cntr = tetra.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
  vol = tetra.volume();
  
  // Computing auxiliary integrals
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  SX2Y = 0;  SXY2 = 0;  SX2Z = 0;  SXZ2 = 0;  SY2Z = 0;  SYZ2 = 0;
  SX3 = 0;  SY3 = 0;  SZ3 = 0;  SXYZ = 0;
  for (i=0; i<GQTe.numberQuadraturePoints(); ++i) {
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQTe.numberVertices(); ++j) {
      qp = qp + GQTe.baryCoordinate(i,j) * tetra.vertex(j);
    }
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;

    SX2 += GQTe.weight(i) * (X*X) * vol;
    SY2 += GQTe.weight(i) * (Y*Y) * vol;
    SZ2 += GQTe.weight(i) * (Z*Z) * vol;

    SXY += GQTe.weight(i) * (X*Y) * vol;
    SXZ += GQTe.weight(i) * (X*Z) * vol;
    SYZ += GQTe.weight(i) * (Y*Z) * vol;

    SX2Y += GQTe.weight(i) * (X*X*Y) * vol;
    SXY2 += GQTe.weight(i) * (X*Y*Y) * vol;

    SX2Z += GQTe.weight(i) * (X*X*Z) * vol;
    SXZ2 += GQTe.weight(i) * (X*Z*Z) * vol;

    SY2Z += GQTe.weight(i) * (Y*Y*Z) * vol;
    SYZ2 += GQTe.weight(i) * (Y*Z*Z) * vol;

    SX3 += GQTe.weight(i) * (X*X*X) * vol;
    SY3 += GQTe.weight(i) * (Y*Y*Y) * vol;
    SZ3 += GQTe.weight(i) * (Z*Z*Z) * vol;
    SXYZ += GQTe.weight(i) * (X*Y*Z) * vol;
  }

  // Assignments
  C(1,1) = SX2;  C(1,2) = SXY;  C(1,3) = SXZ;
  C(2,2) = SY2;  C(2,3) = SYZ;  C(3,3) = SZ2;
  C(2,1) = C(1,2);  C(3,1) = C(1,3);  C(3,2) = C(2,3);

  C(1,4) = SX3;   C(1,5) = SX2Y;  C(1,6) = SX2Z;
  C(2,5) = SXY2;  C(2,6) = SXYZ;  C(3,6) = SXZ2;
  C(2,4) = C(1,5);  C(3,4) = C(1,6);  C(3,5) = C(2,6);

  C(1,7) = SX2Y;  C(1,8) = SXY2;  C(1,9) = SXYZ;
  C(2,8) = SY3;   C(2,9) = SY2Z;  C(3,9) = SYZ2;
  C(2,10) = C(1,11);  C(3,10) = C(1,12);  C(3,11) = C(2,12);

  C(1,10) = SX2Z;  C(1,11) = SXYZ;  C(1,12) = SXZ2;
  C(2,11) = SY2Z;  C(2,12) = SYZ2;  C(3,12) = SZ3;
  C(2,10) = C(1,11);  C(3,10) = C(1,12);  C(3,11) = C(2,12);

  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (displacement gradient)

int Hdiv_TetraRT03_GramMat_NmlzBas(DiagMatrix &GM, const Tetra &tetra,
                                   const GaussQuad &GQTe)
{
  DiagMatrix A(4);
  Hdiv_TetraRT0_GramMat_NmlzBas(A, tetra, GQTe);
  
  GM.resize(12);
  for (int i=0; i<3; ++i) {
    GM.setEntry(4*i+1, A.entry(1));
    GM.setEntry(4*i+2, A.entry(2));
    GM.setEntry(4*i+3, A.entry(3));
    GM.setEntry(4*i+4, A.entry(4));
  }
  
  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (strain)

int Hdiv_TetraRT03_GramMatAvg_NmlzBas(FullMatrix &GMA, const Tetra &tetra,
                                      const GaussQuad &GQTe)
{
  PtVec3d cntr = tetra.center();
  double xc = cntr.xCrd();
  double yc = cntr.yCrd();
  double zc = cntr.zCrd();
  double vol = tetra.volume();

  // Computing auxiliary integrals
  double S1  = vol;
  double SX2 = 0,  SY2 = 0,  SZ2 = 0;
  double SXY = 0,  SXZ = 0,  SYZ = 0;
  for (int i=0; i<GQTe.numberQuadraturePoints(); ++i) {
    PtVec3d qp(0,0,0);
    for (int j=0; j<GQTe.numberVertices(); ++j) {
      qp = qp + GQTe.baryCoordinate(i,j) * tetra.vertex(j);
    }
    double X = qp.xCrd() - xc;
    double Y = qp.yCrd() - yc;
    double Z = qp.zCrd() - zc;
    SX2 += GQTe.weight(i) * (X*X) * vol;
    SY2 += GQTe.weight(i) * (Y*Y) * vol;
    SZ2 += GQTe.weight(i) * (Z*Z) * vol;
    SXY += GQTe.weight(i) * (X*Y) * vol;
    SXZ += GQTe.weight(i) * (X*Z) * vol;
    SYZ += GQTe.weight(i) * (Y*Z) * vol;
  }

  GMA.resize(12,12);
  GMA(1,1) = S1;  GMA(6,6) = S1;  GMA(11,11) = S1;
  GMA(2,2) = 0.5*S1;  GMA(3,3)   = 0.5*S1;
  GMA(5,5) = 0.5*S1;  GMA(7,7)   = 0.5*S1;
  GMA(9,9) = 0.5*S1;  GMA(10,10) = 0.5*S1;
  GMA(4,4)   = SX2 + 0.5*(SY2+SZ2);
  GMA(8,8)   = SY2 + 0.5*(SZ2+SX2);
  GMA(12,12) = SZ2 + 0.5*(SX2+SY2);
  GMA(2,5) = 0.5*S1;  GMA(3,9) = 0.5*S1;  GMA(7,10) = 0.5*S1;
  GMA(5,2) = 0.5*S1;  GMA(9,3) = 0.5*S1;  GMA(10,7) = 0.5*S1;
  GMA(4,8) = 0.5*SXY;  GMA(4,12) = 0.5*SXZ;  GMA(8,12) = 0.5*SYZ;
  GMA(8,4) = 0.5*SXY;  GMA(12,4) = 0.5*SXZ;  GMA(12,8) = 0.5*SYZ;

  return(0);  // If successful
}


////////////////////////////////////////////////////////////////////////////////
// For hexahedra
////////////////////////////////////////////////////////////////////////////////

// Functions for hexahedra

// Hexa AT0: The Gram matrix for (normalized+Piola) basis

int Hdiv_HexaAT0_BasFxnVal_NmlzPiolaBas(PtVec3d *val, const Hexa &hexa,
                                        const double xhat, const double yhat,
                                        const double zhat)
{
  val[0] = PtVec3d(0,0,0);  // Not used
  val[1] = PtVec3d(1,0,0);
  val[2] = PtVec3d(0,1,0);
  val[3] = PtVec3d(0,0,1);
  val[4] = hexa.trilinearmapping(xhat,yhat,zhat) - hexa.center();
  Mat3 PM = hexa.PiolaMatrix(xhat,yhat,zhat);
  val[5] = PM * PtVec3d(xhat,-yhat,0);
  val[6] = PM * PtVec3d(0,yhat,-zhat);
  return(0);  // If successful
}


int Hdiv_HexaAT0_GramMat_NmlzPiolaBas(FullMatrix &GM, const Hexa &hexa,
                                      const GaussQuad &GQH)
{
  PtVec3d *val = new PtVec3d[7];
  GM.resize(6,6);
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150425: REALLY NEEDED
    Hdiv_HexaAT0_BasFxnVal_NmlzPiolaBas(val,hexa,xhat,yhat,zhat);
    for (int i=1; i<=6; ++i) {
      for (int j=1; j<=6; ++j) {
        GM(i,j) = GM(i,j) + GQH.weight(k) * jac * dotProduct(val[i],val[j]);
      }
    }
  }
  return(0);  // If successful
}


// Hexa AT0: The Gram matrix for (normalized+Piola) basis and
// full (3x3) diffusion/permeability matrix/tensor

int Hdiv_HexaAT0_GramMatK_NmlzPiolaBas(FullMatrix &GM, FullMatrix &GMK,
                                       const Hexa &hexa, Mat3 MatK,
                                       const GaussQuad &GQH)
{
  PtVec3d *val = new PtVec3d[7];
  GM.resize(6,6);
  GMK.resize(6,6);
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150425: REALLY NEEDED
    Hdiv_HexaAT0_BasFxnVal_NmlzPiolaBas(val,hexa,xhat,yhat,zhat);
    for (int i=1; i<=6; ++i) {
      for (int j=1; j<=6; ++j) {
        GM(i,j)  = GM(i,j)  + GQH.weight(k) * jac * dotProduct(val[i],val[j]);
        GMK(i,j) = GMK(i,j) + GQH.weight(k) * jac * dotProduct(val[i],MatK*val[j]);
      }
    }
  }
  return(0);  // If successful
}


int Hdiv_HexaAT0_ProjCof_NmlzPiolaBas(FullMatrix &ProjCof,
                                      const Hexa &hexa, Mat3 MatK,
                                      const GaussQuad &GQH)
{
  Vector cof(6), RHS(6);
  FullMatrix GM(6,6), GMK(6,6);
  Hdiv_HexaAT0_GramMatK_NmlzPiolaBas(GM, GMK, hexa, MatK, GQH);
  ProjCof.resize(6,6);
  for (int j=1; j<=6; ++j) {
    for (int i=1; i<=6; ++i)  RHS(i) = GMK(i,j);
    cof = slvFullSpdSysCholesky(GM, RHS);
    for (int i=1; i<=6; ++i)  ProjCof(i,j) = cof(i);
  }
  return(0);  // If successful
}


// ???

int Hdiv_HexaAT0_NmlFlux_NmlzPiolaBas(FullMatrix &NmlFluxAT0,
                                      const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                      const GaussQuad &GQQ)
{
  NmlFluxAT0.resize(6,6);
  return(0);  // if successful
}


////////////////////////////////////////////////////////////////////////////////

// Hexa RT[0]: The Gram matrix for normalized basis

int Hdiv_HexaRT0_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                                 const GaussQuad &GQH)
{
  int j, k;
  double X, x, xc, xhat, Y, y, yc, yhat, Z, z, zc, zhat;
  double SX, SY, SZ, SX2, SY2, SZ2, SXY, SXZ, SYZ;
  double jac, vol;
  PtVec3d cntr, qp;
  
  // Assuming the hexahedron has been enriched
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  
  // Integration using the chosen Gaussian quadrature
  vol = 0;
  SX = 0;  SY = 0;  SZ = 0;
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (k=0; k<GQH.numberQuadraturePoints(); ++k) {
    xhat = GQH.CartesianCoordinate(k,0);
    yhat = GQH.CartesianCoordinate(k,1);
    zhat = GQH.CartesianCoordinate(k,2);
    // std::cout << "xhat,yaht,zhat= " << xhat << "  " << yhat << "  " << zhat << "\n";
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150425: REALLY NEEDED
    // std::cout << "jac=" << jac << "\n";
    qp = hexa.trilinearmapping(xhat,yhat,zhat);
    // An alternative for computing quadrature points
    /*
     qp = PtVec3d(0,0,0);
     for (j=0; j<GQH.numberVertices(); ++j) {
     qp = qp + GQH.baryCoordinate(k,j)*hexa.vertex(j);
     }
     */
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    // std::cout << "X,Y,Z= " << X << "  " << Y << "  " << Z << "\n";
    vol = vol + GQH.weight(k) * jac;
    SX = SX + GQH.weight(k) * jac * X;
    SY = SY + GQH.weight(k) * jac * Y;
    SZ = SZ + GQH.weight(k) * jac * Z;
    SX2 = SX2 + GQH.weight(k) * jac * (X*X);
    SY2 = SY2 + GQH.weight(k) * jac * (Y*Y);
    SZ2 = SZ2 + GQH.weight(k) * jac * (Z*Z);
    SXY = SXY + GQH.weight(k) * jac * (X*Y);
    SXZ = SXZ + GQH.weight(k) * jac * (X*Z);
    SYZ = SYZ + GQH.weight(k) * jac * (Y*Z);
  }
  // std::cout << "SX2=" << SX2 << "\n";
  
  // Assignment for the Gram matrix
  GM.resize(6,6);
  GM(1,1) = vol;  GM(2,2) = vol;  GM(3,3) = vol;
  GM(4,4) = SX2;  GM(5,5) = SY2;  GM(6,6) = SZ2;
  GM(1,4) = SX;  GM(2,5) = SY;  GM(3,6) = SZ;
  GM(4,1) = GM(1,4);  GM(5,2) = GM(2,5);  GM(6,3) = GM(3,6);
  // std::cout << GM << "\n";
  
  return(0);  // if successful
}


// Hexa RT[0]: The Gram matrix for normalized basis and
// full (3x3) diffusion/permeability matrix/tensor

int Hdiv_HexaRT0_GramMatK_NmlzBas(FullMatrix &GM, FullMatrix &GMK,
                                  const Hexa &hexa, Mat3 MatK,
                                  const GaussQuad &GQH)
{
  int i, j, k;
  double X, x, xc, xhat, Y, y, yc, yhat, Z, z, zc, zhat;
  double SX, SY, SZ, SX2, SY2, SZ2, SXY, SXZ, SYZ;
  double jac, vol;
  PtVec3d cntr, qp;
  
  // Assuming the hexahedron has been enriched
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  
  // Integration using the chosen Gaussian quadrature
  vol = 0;
  SX = 0;  SY = 0;  SZ = 0;
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (k=0; k<GQH.numberQuadraturePoints(); ++k) {
    xhat = GQH.CartesianCoordinate(k,0);
    yhat = GQH.CartesianCoordinate(k,1);
    zhat = GQH.CartesianCoordinate(k,2);
    // std::cout << "xhat,yaht,zhat= " << xhat << "  " << yhat << "  " << zhat << "\n";
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150425: REALLY NEEDED
    // std::cout << "jac=" << jac << "\n";
    qp = hexa.trilinearmapping(xhat,yhat,zhat);
    // An alternative for computing quadrature points
    /*
     qp = PtVec3d(0,0,0);
     for (j=0; j<GQH.numberVertices(); ++j) {
     qp = qp + GQH.baryCoordinate(k,j)*hexa.vertex(j);
     }
     */
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    // std::cout << "X,Y,Z= " << X << "  " << Y << "  " << Z << "\n";
    vol = vol + GQH.weight(k) * jac;
    SX = SX + GQH.weight(k) * jac * X;
    SY = SY + GQH.weight(k) * jac * Y;
    SZ = SZ + GQH.weight(k) * jac * Z;
    SX2 = SX2 + GQH.weight(k) * jac * (X*X);
    SY2 = SY2 + GQH.weight(k) * jac * (Y*Y);
    SZ2 = SZ2 + GQH.weight(k) * jac * (Z*Z);
    SXY = SXY + GQH.weight(k) * jac * (X*Y);
    SXZ = SXZ + GQH.weight(k) * jac * (X*Z);
    SYZ = SYZ + GQH.weight(k) * jac * (Y*Z);
  }
  // std::cout << "SX2=" << SX2 << "\n";
  
  // Assignment for the Gram matrix
  GM.resize(6,6);
  GM(1,1) = vol;  GM(2,2) = vol;  GM(3,3) = vol;
  GM(4,4) = SX2;  GM(5,5) = SY2;  GM(6,6) = SZ2;
  GM(1,4) = SX;  GM(2,5) = SY;  GM(3,6) = SZ;
  GM(4,1) = GM(1,4);  GM(5,2) = GM(2,5);  GM(6,3) = GM(3,6);
  // std::cout << GM << "\n";
  
  // std::cout << MatK << "\n";
  
  // Assignment for Gram matrix with the diffusion/permeability matrix MatK
  GMK.resize(6,6);
  //
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GMK(i,j) = MatK(i,j) * vol;
  //
  GMK(1,4) = MatK(1,1)*SX;  GMK(1,5) = MatK(1,2)*SY;  GMK(1,6) = MatK(1,3)*SZ;
  GMK(2,4) = MatK(2,1)*SX;  GMK(2,5) = MatK(2,2)*SY;  GMK(2,6) = MatK(2,3)*SZ;
  GMK(3,4) = MatK(3,1)*SX;  GMK(3,5) = MatK(3,2)*SY;  GMK(3,6) = MatK(3,3)*SZ;
  for (i=4; i<=6; ++i)
    for (j=1; j<=3; ++j)
      GMK(i,j) = GMK(j,i);
  //
  GMK(4,4) = MatK(1,1)*SX2;  GMK(4,5) = MatK(1,2)*SXY;  GMK(4,6) = MatK(1,3)*SXZ;
  GMK(5,5) = MatK(2,2)*SY2;  GMK(5,6) = MatK(2,3)*SYZ;  GMK(6,6) = MatK(3,3)*SZ2;
  GMK(5,4) = GMK(4,5);  GMK(6,4) = GMK(4,6);  GMK(6,5) = GMK(5,6);
  // std::cout << GMK << "\n";
  
  return(0);  // if successful
}


// ???

int Hdiv_HexaRT0_NmlFlux_NmlzBas(FullMatrix &NmlFluxRT0,
                                 const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                 const GaussQuad &GQQ)
{
  int i, j, k;
  double dp, jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat;
  PtVec3d cntr, nml, qp, w[7];  // w[0] not used
  Quadri3d quadri;

  //
  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
  
  //
  NmlFluxRT0.resize(6,6);
  w[1] = PtVec3d(1,0,0);
  w[2] = PtVec3d(0,1,0);
  w[3] = PtVec3d(0,0,1);
  for (i=1; i<=6; ++i) {  // the i-th face
    quadri = SQF[i-1];
    for (j=1; j<=6; ++j)  NmlFluxRT0(i,j) = 0;
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      qp = quadri.bilinearmapping(xhat, yhat);
      nml = quadri.normal(xhat, yhat);
      x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
      X = x - xc;  Y = y - yc;  Z = z - zc;
      w[4] = PtVec3d(X, 0, 0);
      w[5] = PtVec3d(0, Y, 0);
      w[6] = PtVec3d(0, 0, Z);
      for (j=1; j<=6; ++j) {
        dp = - dotProduct(w[j], nml);  // Notice the negative sign!
        NmlFluxRT0(i,j) += GQQ.weight(k) * jac * dp;
      }
    }
    for (j=1; j<=6; ++j)
      NmlFluxRT0(i,j) = NmlFluxRT0(i,j) * sign[i-1];
  }

  return(0);  // if successful
}


// Hdiv: Hexa0 & Hexa1: K-transformation matrix from RT[0] to P_1^3

int Hdiv_Hexa01_KtransMat(FullMatrix &KtransMat, Mat3 &MatK)
{
  int i, j;

  KtransMat.resize(12,6);
  
  // Computing the KtransMat(12,6) that already contains the NEGATIVE sign
  // the -K transformation image of the j-th RT[0] bas.fxn.
  // as a linear comb. of 12 P_1^3 bas.fxns.
  for (j=1; j<=3; ++j)
    for (i=1; i<=3; ++i)
      KtransMat(i,j) = -MatK(i,j);
  //
  KtransMat(4,4) = -MatK(1,1);  KtransMat(5,5) = -MatK(1,2);  KtransMat(6,6) = -MatK(1,3);
  KtransMat(7,4) = -MatK(2,1);  KtransMat(8,5) = -MatK(2,2);  KtransMat(9,6) = -MatK(2,3);
  KtransMat(10,4)= -MatK(3,1);  KtransMat(11,5)= -MatK(3,2);  KtransMat(12,6)= -MatK(3,3);

  return(0);  // if successful
}


// Hdiv: Hexa0 & Hexa1: Projecting (12) P_1^3 basis functions into RT[0]

int Hdiv_Hexa01_ProjP1v3RT0(FullMatrix &ProjMat, const Hexa &hexa,
                            const GaussQuad &GQH)
{
  int i, j, k;
  double jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat, zhat;
  double SX, SY, SZ, SX2, SY2, SZ2, SXY, SXZ, SYZ;
  PtVec3d cntr, qp;  // w[0] not used
  Vector cof(6), RHS(6);
  FullMatrix GM(6,6), GMK(6,6);

  ProjMat.resize(6,12);

  //
  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();

  //
  Hdiv_HexaRT0_GramMat_NmlzBas(GM, hexa, GQH);
  
  // Computing auxiliary integrals on the hexahedral element
  SX = 0;  SY = 0;  SZ = 0;
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    // Alternative:
    /*
     qp = PtVec3d(0,0,0);
     for (j=0; j<GQH.numberVertices(); ++j) {
     qp = qp + GQH.baryCoordinate(i,j) * hexa.vertex(j);
     }
     */
    qp = hexa.trilinearmapping(xhat, yhat, zhat);
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    // S += GQH.weight(i) * jac;
    SX += GQH.weight(i) * jac * X;
    SY += GQH.weight(i) * jac * Y;
    SZ += GQH.weight(i) * jac * Z;
    SX2 += GQH.weight(i) * jac * (X*X);
    SY2 += GQH.weight(i) * jac * (Y*Y);
    SZ2 += GQH.weight(i) * jac * (Z*Z);
    SXY += GQH.weight(i) * jac * (X*Y);
    SXZ += GQH.weight(i) * jac * (X*Z);
    SYZ += GQH.weight(i) * jac * (Y*Z);
  }
  
  // Computing the ProjMat(6,12) matrix
  ProjMat(1,1) = 1;  ProjMat(2,2) = 1;  ProjMat(3,3) = 1;
  ProjMat(4,4) = 1;  ProjMat(5,8) = 1;  ProjMat(6,12) = 1;
  //
  // For the rest, 6-by-6 linear systems need to be solved
  //
  RHS.resize(6);
  RHS(2) = SX;  RHS(5) = SXY;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,5) = cof(i);
  //
  RHS.resize(6);
  RHS(3) = SX;  RHS(6) = SXZ;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,6) = cof(i);
  //
  RHS.resize(6);
  RHS(1) = SY;  RHS(4) = SXY;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,7) = cof(i);
  //
  RHS.resize(6);
  RHS(3) = SY;  RHS(6) = SYZ;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,9) = cof(i);
  //
  RHS.resize(6);
  RHS(1) = SZ;  RHS(4) = SXZ;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,10) = cof(i);
  //
  RHS.resize(6);
  RHS(2) = SZ;  RHS(5) = SYZ;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=6; ++i)  ProjMat(i,11) = cof(i);
  
  return(0);  // if successful
}


// Hexa BDM[1]: The Gram matrix for normalized basis functions and
// full (3x3) diffusion/permeability matrix/tensor

int Hdiv_HexaBDM1_GramMatK_NmlzBas(FullMatrix &GM, FullMatrix &GMK,
                                   const Hexa &hexa, Mat3 MatK,
                                   const GaussQuad &GQH)
{
  int i, j;
  double IntgrlX, IntgrlY, IntgrlZ;
  double IntgrlX2, IntgrlY2, IntgrlZ2;
  double IntgrlYZ, IntgrlZX, IntgrlXY;
  double jac, vol;
  double X, Y, Z;
  double x, y, z;
  double xc, yc, zc;
  double xhat, yhat, zhat;
  PtVec3d cntr, qp;
  
  // Assuming the hexahedron has been enriched
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  
  // Computing integrals using Gaussian quadrature
  // vol = 0;
  IntgrlX = 0;  IntgrlX2 = 0;  IntgrlYZ = 0;
  IntgrlY = 0;  IntgrlY2 = 0;  IntgrlZX = 0;
  IntgrlZ = 0;  IntgrlZ2 = 0;  IntgrlXY = 0;
  for (i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQH.numberVertices(); ++j) {
      qp = qp + GQH.baryCoordinate(i,j)*hexa.vertex(j);
    }
    // qp = hexa.trilinearmapping(xhat,yhat,zhat);
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    // vol = vol + jac*GQH.weight(i);
    IntgrlX = IntgrlX + X*jac*GQH.weight(i);
    IntgrlY = IntgrlY + Y*jac*GQH.weight(i);
    IntgrlZ = IntgrlZ + Z*jac*GQH.weight(i);
    IntgrlX2 = IntgrlX2 + (X*X)*jac*GQH.weight(i);
    IntgrlY2 = IntgrlY2 + (Y*Y)*jac*GQH.weight(i);
    IntgrlZ2 = IntgrlZ2 + (Z*Z)*jac*GQH.weight(i);
    IntgrlYZ = IntgrlYZ + (Y*Z)*jac*GQH.weight(i);
    IntgrlZX = IntgrlZX + (Z*X)*jac*GQH.weight(i);
    IntgrlXY = IntgrlXY + (X*Y)*jac*GQH.weight(i);
  }
  
  // Assignment for Gram matrix
  // The same as that for RT[0]
  GM.resize(12,12);
  GM(1,1) = vol;  GM(2,2) = vol;  GM(3,3) = vol;
  GM(4,4) = IntgrlX2;  GM(5,5) = IntgrlY2;  GM(6,6) = IntgrlZ2;
  GM(1,4) = IntgrlX;  GM(2,5) = IntgrlY;  GM(3,6) = IntgrlZ;
  GM(4,1) = GM(1,4);  GM(5,2) = GM(2,5);  GM(6,3) = GM(3,6);
  // New stuff for BDM[1] in the hierarchy
  GM(1,7)  = IntgrlY;  GM(1,8)  = IntgrlZ;
  GM(2,9)  = IntgrlZ;  GM(2,10) = IntgrlX;
  GM(3,11) = IntgrlX;  GM(3,12) = IntgrlY;
  //
  GM(4,7)  = IntgrlXY;  GM(4,8)  = IntgrlZX;
  GM(5,9)  = IntgrlYZ;  GM(5,10) = IntgrlZX;
  GM(6,11) = IntgrlZX;  GM(6,12) = IntgrlYZ;
  //
  GM(7,1) = GM(1,7);  GM(8,1)  = GM(1,8);
  GM(9,2) = GM(2,9);  GM(10,2) = GM(2,10);
  GM(11,3) = GM(3,11);  GM(12,3) = GM(3,12);
  //
  GM(7,4)  = GM(4,7);   GM(8,4)  = GM(4,8);
  GM(9,5)  = GM(5,9);   GM(10,5) = GM(5,10);
  GM(11,6) = GM(6,11);  GM(12,6) = GM(6,12);
  //
  GM(7,7) = IntgrlY2;  GM(7,8) = IntgrlYZ;
  GM(8,7) = GM(7,8);  GM(8,8) = IntgrlZ2;
  //
  GM(9,9) = IntgrlZ2;  GM(9,10) = IntgrlZX;
  GM(10,9) = GM(9,10);  GM(10,10) = IntgrlX2;
  //
  GM(11,11) = IntgrlX2;  GM(11,12) = IntgrlXY;
  GM(12,11) = GM(11,12);  GM(12,12) = IntgrlY2;
  
  // Assignment for Gram matrix with permeability matrix PermK
  GMK.resize(12,12);
  
  return(0);  // if successful
}


int Hdiv_Hexa1_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                               const GaussQuad &GQH)
{
  int i, j;
  double jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat, zhat;
  double S, SX, SY, SZ, SX2, SY2, SZ2, SXY, SXZ, SYZ;
  PtVec3d cntr, qp;
  FullMatrix MatI(3,3);
  FullMatrix MatA, MatB;  // Their sizes are to be determined later

  // Assuming the hexahedron has been enriched
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();

  // Computing integrals using the chosen Gaussian quadrature
  S = 0;
  SX = 0;  SY = 0;  SZ = 0;
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  // std::cout << "GQH.numberQuadraturePoints()= "
  //   << GQH.numberQuadraturePoints() << "\n";
  for (i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    // Alternative:
    /*
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQH.numberVertices(); ++j) {
      qp = qp + GQH.baryCoordinate(i,j) * hexa.vertex(j);
    }
    */
    qp = hexa.trilinearmapping(xhat, yhat, zhat);
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    S += GQH.weight(i) * jac;
    SX += GQH.weight(i) * jac * X;
    SY += GQH.weight(i) * jac * Y;
    SZ += GQH.weight(i) * jac * Z;
    SX2 += GQH.weight(i) * jac * (X*X);
    SY2 += GQH.weight(i) * jac * (Y*Y);
    SZ2 += GQH.weight(i) * jac * (Z*Z);
    SXY += GQH.weight(i) * jac * (X*Y);
    SXZ += GQH.weight(i) * jac * (X*Z);
    SYZ += GQH.weight(i) * jac * (Y*Z);
  }
  
  // JL20150520: TO BE REVISED: mat3 as a subclass of FullMatrix
  MatI(1,1) = 1;  MatI(2,2) = 1;  MatI(3,3) = 1;
  
  // Assignments for the Gram matrix
  GM.resize(12,12);
  //
  GM(1,1) = S;  GM(2,2) = S;  GM(3,3) = S;
  //
  MatA.resize(1,3);
  MatA(1,1) = SX;  MatA(1,2) = SY;  MatA(1,3) = SZ;
  MatB = tensorProduct(MatI, MatA);
  for (i=1; i<=3; ++i) {
    for (j=1; j<=9; ++j) {
      GM(i,j+3) = MatB(i,j);
      GM(j+3,i) = MatB(i,j);
    }
  }
  //
  MatA.resize(3,3);
  MatA(1,1) = SX2;  MatA(1,2) = SXY;  MatA(1,3) = SXZ;
  MatA(2,1) = SXY;  MatA(2,2) = SY2;  MatA(2,3) = SYZ;
  MatA(3,1) = SXZ;  MatA(3,2) = SYZ;  MatA(3,3) = SZ2;
  MatB = tensorProduct(MatI, MatA);
  for (i=1; i<=9; ++i)
    for (j=1; j<=9; ++j)
      GM(i+3,j+3) = MatB(i,j);
  
  return(0);  // if successful
}


// ???

int Hdiv_Hexa1_GramMatK_NmlzBas(FullMatrix &GMK, const Hexa &hexa,
                                Mat3 &MatK, const GaussQuad &GQH)
{
  int i, j;
  double jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat, zhat;
  double S, SX, SY, SZ, SX2, SY2, SZ2, SXY, SXZ, SYZ;
  PtVec3d cntr, qp;
  FullMatrix MatKK(3,3);
  FullMatrix MatA, MatB;  // Their sizes are to be determined later

  // Assuming the hexahedron has been enriched
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();

  // Computing integrals using the chosen Gaussian quadrature
  S = 0;
  SX = 0;  SY = 0;  SZ = 0;
  SX2 = 0;  SY2 = 0;  SZ2 = 0;
  SXY = 0;  SXZ = 0;  SYZ = 0;
  for (int i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    // Alternative:
    /*
    qp = PtVec3d(0,0,0);
    for (j=0; j<GQH.numberVertices(); ++j) {
      qp = qp + GQH.baryCoordinate(i,j) * hexa.vertex(j);
    }
    */
    qp = hexa.trilinearmapping(xhat, yhat, zhat);
    x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
    X = x - xc;  Y = y - yc;  Z = z - zc;
    S = S + GQH.weight(i) * jac;
    SX = SX + GQH.weight(i) * jac * X;
    SY = SY + GQH.weight(i) * jac * Y;
    SZ = SZ + GQH.weight(i) * jac * Z;
    SX2 = SX2 + GQH.weight(i) * jac * (X*X);
    SY2 = SY2 + GQH.weight(i) * jac * (Y*Y);
    SZ2 = SZ2 + GQH.weight(i) * jac * (Z*Z);
    SXY = SXY + GQH.weight(i) * jac * (X*Y);
    SXZ = SXZ + GQH.weight(i) * jac * (X*Z);
    SYZ = SYZ + GQH.weight(i) * jac * (Y*Z);
  }
  
  // JL20150520: TO BE REVISED: mat3 as a subclass of FullMatrix
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      MatKK(i,j) = MatK(i,j);
  
  // Assignments for the Gram matrix with permeability matrix MatK
  GMK.resize(12,12);
  //
  for (i=1; i<=3; ++i)
    for (j=1; j<=3; ++j)
      GMK(i,j) = MatKK(i,j) * S;
  //
  MatA.resize(1,3);
  MatA(1,1) = SX;  MatA(1,2) = SY;  MatA(1,3) = SZ;
  MatB = tensorProduct(MatKK, MatA);
  for (i=1; i<=3; ++i) {
    for (j=1; j<=9; ++j) {
      GMK(i,j+3) = MatB(i,j);
      GMK(j+3,i) = MatB(i,j);
    }
  }
  //
  MatA.resize(3,3);
  MatA(1,1) = SX2;  MatA(1,2) = SXY;  MatA(1,3) = SXZ;
  MatA(2,1) = SXY;  MatA(2,2) = SY2;  MatA(2,3) = SYZ;
  MatA(3,1) = SXZ;  MatA(3,2) = SYZ;  MatA(3,3) = SZ2;
  MatB = tensorProduct(MatKK, MatA);
  for (i=1; i<=9; ++i)
    for (j=1; j<=9; ++j)
      GMK(i+3,j+3) = MatB(i,j);
  
  return(0);  // if successful
}


// Matrix-version, to be used for elasticity (displacement gradient)
// Actually GM should be a block diagonal matrix
// with each block having size 6x6

int Hdiv_HexaAT03_GramMat_NmlzPiolaBas(FullMatrix &GM, const Hexa &hexa,
                                       const GaussQuad &GQH)
{
  FullMatrix A(6,6);
  Hdiv_HexaAT0_GramMat_NmlzPiolaBas(A, hexa, GQH);
  
  GM.resize(18,18);
  for (int i=1; i<=6; ++i) {
    for (int j=1; j<=6; ++j) {
      GM(i,   j   ) = A(i,j);
      GM(i+6, j+6 ) = A(i,j);
      GM(i+12,j+12) = A(i,j);
    }
  }
  
  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (strain)

int Hdiv_HexaAT03_GramMatAvg_NmlzPiolaBas(FullMatrix &GMA, const Hexa &hexa,
                                          const GaussQuad &GQH)
{
  // Assuming the hexahedron has been enriched
  PtVec3d cntr = hexa.center();
  double xc = cntr.xCrd();
  double yc = cntr.yCrd();
  double zc = cntr.zCrd();

  // Computing integrals using the chosen Gaussian quadrature
  double S1 = 0,  SX = 0,  SY = 0,  SZ = 0;
  double SX2 = 0,  SY2 = 0, SZ2 = 0;
  double SXY = 0,  SXZ = 0, SYZ = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    PtVec3d qp = hexa.trilinearmapping(xhat, yhat, zhat);
    double x = qp.xCrd();
    double y = qp.yCrd();
    double z = qp.zCrd();
    double X = x - xc;
    double Y = y - yc;
    double Z = z - zc;
    S1 = S1 + GQH.weight(k) * jac;
    SX = SX + GQH.weight(k) * jac * X;
    SY = SY + GQH.weight(k) * jac * Y;
    SZ = SZ + GQH.weight(k) * jac * Z;
    SX2 = SX2 + GQH.weight(k) * jac * (X*X);
    SY2 = SY2 + GQH.weight(k) * jac * (Y*Y);
    SZ2 = SZ2 + GQH.weight(k) * jac * (Z*Z);
    SXY = SXY + GQH.weight(k) * jac * (X*Y);
    SXZ = SXZ + GQH.weight(k) * jac * (X*Z);
    SYZ = SYZ + GQH.weight(k) * jac * (Y*Z);
  }
  
  GMA.resize(18,18);
  // Main diagonal block #1
  GMA(1,1) = S1;   GMA(2,2) = 0.5*S1;   GMA(3,3) = 0.5*S1;
  GMA(4,4) = SX2;  GMA(5,5) = 0.5*SY2;  GMA(6,6) = 0.5*SZ2;
  GMA(1,4) = SX;   GMA(2,5) = 0.5*SY;   GMA(3,6) = 0.5*SZ;
  GMA(4,1) = SX;   GMA(5,2) = 0.5*SY;   GMA(6,3) = 0.5*SZ;
  // Main diagonal block #2
  GMA(7,7)   = 0.5*S1;   GMA(8,8)   = S1;   GMA(9,9)   = 0.5*S1;
  GMA(10,10) = 0.5*SX2;  GMA(11,11) = SY2;  GMA(12,12) = 0.5*SZ2;
  GMA(7,10)  = 0.5*SX;   GMA(8,11)  = SY;   GMA(9,12)  = 0.5*SZ;
  GMA(10,7)  = 0.5*SX;   GMA(11,8)  = SY;   GMA(12,9)  = 0.5*SZ;
  // Main diagonal block #3
  GMA(13,13) = 0.5*S1;   GMA(14,14) = 0.5*S1;   GMA(15,15) = S1;
  GMA(16,16) = 0.5*SX2;  GMA(17,17) = 0.5*SY2;  GMA(18,18) = SZ2;
  GMA(13,16) = 0.5*SX;   GMA(14,17) = 0.5*SY;   GMA(15,18) = SZ;
  GMA(16,13) = 0.5*SX;   GMA(17,14) = 0.5*SY;   GMA(18,15) = SZ;
  // Off-diagonal block #(1,2), #(2,1)
  GMA(2,7) = 0.5*S1;  GMA(2,10) = 0.5*SX;  GMA(5,7) = 0.5*SY;  GMA(5,10) = 0.5*SXY;
  GMA(7,2) = 0.5*S1;  GMA(10,2) = 0.5*SX;  GMA(7,5) = 0.5*SY;  GMA(10,5) = 0.5*SXY;
  // Off-diagonal block #(1,3), #(3,1)
  GMA(3,13) = 0.5*S1;  GMA(3,16) = 0.5*SX;  GMA(6,13) = 0.5*SZ;  GMA(6,16) = 0.5*SXZ;
  GMA(13,3) = 0.5*S1;  GMA(16,3) = 0.5*SX;  GMA(13,6) = 0.5*SZ;  GMA(16,6) = 0.5*SXZ;
  // Off-diagonal block #(2,3), #(3,2)
  GMA(9,14) = 0.5*S1;  GMA(9,17) = 0.5*SY;  GMA(12,14) = 0.5*SZ;  GMA(12,17) = 0.5*SYZ;
  GMA(14,9) = 0.5*S1;  GMA(17,9) = 0.5*SY;  GMA(14,12) = 0.5*SZ;  GMA(17,12) = 0.5*SYZ;
  
  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (displacement gradient)
// Actually GM should be a block diagonal matrix 
// with each block having size 6x6 

int Hdiv_HexaRT03_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                                  const GaussQuad &GQH) 
{
  FullMatrix A(6,6);
  Hdiv_HexaRT0_GramMat_NmlzBas(A, hexa, GQH);

  GM.resize(18,18);
  for (int i=1; i<=6; ++i) {
	  for (int j=1; j<=6; ++j) {
	    GM(i,   j   ) = A(i,j);
	    GM(i+6, j+6 ) = A(i,j);
	    GM(i+12,j+12) = A(i,j);
	  }
  }

  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (strain)

int Hdiv_HexaRT03_GramMatAvg_NmlzBas(FullMatrix &GMA, const Hexa &hexa,
                                     const GaussQuad &GQH)
{
  // Assuming the hexahedron has been enriched
  PtVec3d cntr = hexa.center();
  double xc = cntr.xCrd();
  double yc = cntr.yCrd();
  double zc = cntr.zCrd();

  // Computing integrals using the chosen Gaussian quadrature
  double S1 = 0,  SX = 0,  SY = 0,  SZ = 0;
  double SX2 = 0,  SY2 = 0, SZ2 = 0;
  double SXY = 0,  SXZ = 0, SYZ = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    /*
     qp = PtVec3d(0,0,0);
     for (j=0; j<GQH.numberVertices(); ++j) {
     qp = qp + GQH.baryCoordinate(i,j) * hexa.vertex(j);
     }
     */
    PtVec3d qp = hexa.trilinearmapping(xhat, yhat, zhat);
    double x = qp.xCrd();
    double y = qp.yCrd();
    double z = qp.zCrd();
    double X = x - xc;
    double Y = y - yc;
    double Z = z - zc;
    S1 = S1 + GQH.weight(k) * jac;
    SX = SX + GQH.weight(k) * jac * X;
    SY = SY + GQH.weight(k) * jac * Y;
    SZ = SZ + GQH.weight(k) * jac * Z;
    SX2 = SX2 + GQH.weight(k) * jac * (X*X);
    SY2 = SY2 + GQH.weight(k) * jac * (Y*Y);
    SZ2 = SZ2 + GQH.weight(k) * jac * (Z*Z);
    SXY = SXY + GQH.weight(k) * jac * (X*Y);
    SXZ = SXZ + GQH.weight(k) * jac * (X*Z);
    SYZ = SYZ + GQH.weight(k) * jac * (Y*Z);
  }
  
  GMA.resize(18,18);
  // Main diagonal block #1
  GMA(1,1) = S1;   GMA(2,2) = 0.5*S1;   GMA(3,3) = 0.5*S1;
  GMA(4,4) = SX2;  GMA(5,5) = 0.5*SY2;  GMA(6,6) = 0.5*SZ2;
  GMA(1,4) = SX;   GMA(2,5) = 0.5*SY;   GMA(3,6) = 0.5*SZ;
  GMA(4,1) = SX;   GMA(5,2) = 0.5*SY;   GMA(6,3) = 0.5*SZ;
  // Main diagonal block #2
  GMA(7,7)   = 0.5*S1;   GMA(8,8)   = S1;   GMA(9,9)   = 0.5*S1;
  GMA(10,10) = 0.5*SX2;  GMA(11,11) = SY2;  GMA(12,12) = 0.5*SZ2;
  GMA(7,10)  = 0.5*SX;   GMA(8,11)  = SY;   GMA(9,12)  = 0.5*SZ;
  GMA(10,7)  = 0.5*SX;   GMA(11,8)  = SY;   GMA(12,9)  = 0.5*SZ;
  // Main diagonal block #3
  GMA(13,13) = 0.5*S1;   GMA(14,14) = 0.5*S1;   GMA(15,15) = S1;
  GMA(16,16) = 0.5*SX2;  GMA(17,17) = 0.5*SY2;  GMA(18,18) = SZ2;
  GMA(13,16) = 0.5*SX;   GMA(14,17) = 0.5*SY;   GMA(15,18) = SZ;
  GMA(16,13) = 0.5*SX;   GMA(17,14) = 0.5*SY;   GMA(18,15) = SZ;
  // Off-diagonal block #(1,2), #(2,1)
  GMA(2,7) = 0.5*S1;  GMA(2,10) = 0.5*SX;  GMA(5,7) = 0.5*SY;  GMA(5,10) = 0.5*SXY;
  GMA(7,2) = 0.5*S1;  GMA(10,2) = 0.5*SX;  GMA(7,5) = 0.5*SY;  GMA(10,5) = 0.5*SXY;
  // Off-diagonal block #(1,3), #(3,1)
  GMA(3,13) = 0.5*S1;  GMA(3,16) = 0.5*SX;  GMA(6,13) = 0.5*SZ;  GMA(6,16) = 0.5*SXZ;
  GMA(13,3) = 0.5*S1;  GMA(16,3) = 0.5*SX;  GMA(13,6) = 0.5*SZ;  GMA(16,6) = 0.5*SXZ;
  // Off-diagonal block #(2,3), #(3,2)
  GMA(9,14) = 0.5*S1;  GMA(9,17) = 0.5*SY;  GMA(12,14) = 0.5*SZ;  GMA(12,17) = 0.5*SYZ;
  GMA(14,9) = 0.5*S1;  GMA(17,9) = 0.5*SY;  GMA(14,12) = 0.5*SZ;  GMA(17,12) = 0.5*SYZ;
  
  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (displacement gradient)

int Hdiv_HexaP033_GramMat_NmlzBas(DiagMatrix &GM, const Hexa &hexa,
                                  const GaussQuad &GQH)
{
  // Computing hexa. volume using the chosen Gaussian quadrature
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    vol += GQH.weight(k) * jac;
  }

  GM.resize(9);
  for (int i=1; i<=9; ++i)  GM.setEntry(i,vol);
  
  return(0);  // If successful
}


// Matrix-version, to be used for elasticity (strain)

int Hdiv_HexaP033_GramMatAvg_NmlzBas(FullMatrix &GMA, const Hexa &hexa,
                                     const GaussQuad &GQH)
{
  // Computing hexa. volume using the chosen Gaussian quadrature
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // It is safer to use fabs
    vol += GQH.weight(k) * jac;
  }
  double hvol = 0.5*vol;

  GMA.resize(9,9);
  GMA(1,1) = vol;  GMA(5,5) = vol;  GMA(9,9) = vol;
  GMA(2,2) = hvol;  GMA(3,3) = hvol;
  GMA(4,4) = hvol;  GMA(6,6) = hvol;
  GMA(7,7) = hvol;  GMA(8,8) = hvol;
  GMA(2,4) = hvol;  GMA(3,7) = hvol;  GMA(6,8) = hvol;
  GMA(4,2) = hvol;  GMA(7,3) = hvol;  GMA(8,6) = hvol;

  return(0);  // If successful
}

// Hdiv3d.cpp
