// WG3d_Hexa.cpp
// This file will be expanded as we develop more weak Galerkin finite elements 
// James Liu, ColoState; 2014/07--2017/11

// Weak Galerkin 3-dim finite elements on hexahedra:
// For Darcy equation
//   (P0,P0;AT0):                   Good for any convex hexahedron with flat faces
//   (Q0,Q0;RT[0]):                 Good if not too much distortion
//   (Q0,Q1;P1^3):                  Unsuccessful so far
//   (P1,P0;P0^3) with stabilizer:  Theoretically good 
// For elasticity equation 
//   (P0^3,P0^3;AT0^3,P0):          Good for any convex hexahedron with flat faces
//   (Q0^3,Q0^3;RT[0]^3,Q0):        Good if not too much distortion

#include <iostream>

#include "LinSys.h"
#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "mat3.h"
#include "PtVec3d.h"


////////////////////////////////////////////////////////////////////////////////
// For WG(P0,P0;AT0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P0,P0;AT0): Element mass matrix (1-by-1, or just a number)

FullMatrix WG3d_HexaP0P0AT0_EltMassMat(const Hexa &hexa, const GaussQuad &GQH)
{
  double jac, vol, xhat, yhat, zhat;
  FullMatrix EMM(1,1);
  vol = 0;
  for (int i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol = vol + GQH.weight(i)*jac;
  }
  EMM(1,1) = vol;
  return EMM;
}


// WG3d: Hexa(P0,P0;AT0): Coeffs. of disc.wk.grad. in AT0 nmlz.Piola.bas. (7-by-6)

FullMatrix WG3d_HexaP0P0AT0_CofDiscWkGrad_NmlzPiolaBas(const Hexa &hexa,
                                                       const GaussQuad &GQH,
                                                       const GaussQuad &GQQ)
{
  int sign[6];
  Quadri3d SQF[6];
  Vector cof(6), RHS(6);
  FullMatrix GM(6,6), GMK(6,6), MatCDWG(7,6), MatRHS(7,6);

  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol = vol + GQH.weight(k) * jac;
  }

  // Center
  PtVec3d EmCntr = hexa.center();

  // S? quadrilateral faces and their signes
  SQF[0] = Quadri3d(hexa.vertex(0),hexa.vertex(3),hexa.vertex(7),hexa.vertex(4));
  SQF[1] = Quadri3d(hexa.vertex(1),hexa.vertex(2),hexa.vertex(6),hexa.vertex(5));
  SQF[2] = Quadri3d(hexa.vertex(0),hexa.vertex(4),hexa.vertex(5),hexa.vertex(1));
  SQF[3] = Quadri3d(hexa.vertex(3),hexa.vertex(7),hexa.vertex(6),hexa.vertex(2));
  SQF[4] = Quadri3d(hexa.vertex(0),hexa.vertex(1),hexa.vertex(2),hexa.vertex(3));
  SQF[5] = Quadri3d(hexa.vertex(4),hexa.vertex(5),hexa.vertex(6),hexa.vertex(7));
  sign[0] = -1;  sign[1] = 1;
  sign[2] = -1;  sign[3] = 1;
  sign[4] = -1;  sign[5] = 1;
  
  // For the interior WG basis functions
  MatRHS(1,4) = -3*vol;

  // For the face WG functions
  for (int j=1; j<=6; ++j) {
    PtVec3d FcCntr = SQF[j-1].center();
    PtVec3d w = FcCntr - EmCntr;
    PtVec3d nml = sign[j-1] * SQF[j-1].normal(0.5,0.5);
    double ar = SQF[j-1].area();  // Assuming it is flat
    MatRHS(j+1,1) = nml.xCrd() * ar;
    MatRHS(j+1,2) = nml.yCrd() * ar;
    MatRHS(j+1,3) = nml.zCrd() * ar;
    MatRHS(j+1,4) = dotProduct(w,nml) * ar;
  }
  MatRHS(3,5) =  1;  // Face #2
  MatRHS(5,5) = -1;  // Face #4
  MatRHS(5,6) =  1;  // Face #4
  MatRHS(7,6) = -1;  // Face #6
  
  Hdiv_HexaAT0_GramMat_NmlzPiolaBas(GM, hexa, GQH);
  
  for (int i=1; i<=7; ++i) {
    for (int j=1; j<=6; ++j)  RHS(j) = MatRHS(i,j);
    cof = slvFullSpdSysCholesky(GM, RHS);
    for (int j=1; j<=6; ++j)  MatCDWG(i,j) = cof(j);
  }

  return MatCDWG;
}


// WG3d: Hexa(P0,P0;AT0): Element grad-grad matrix (7-by-7)
// with diffusion/permeability (3x3) matrix MatK

FullMatrix WG3d_HexaP0P0AT0_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ)
{
  Vector cofi(6), cofj(6);
  FullMatrix EGGMK(7,7), GM(6,6), GMK(6,6), MatCDWG(7,6);
  
  Hdiv_HexaAT0_GramMatK_NmlzPiolaBas(GM, GMK, hexa, MatK, GQH);
  MatCDWG = WG3d_HexaP0P0AT0_CofDiscWkGrad_NmlzPiolaBas(hexa, GQH, GQQ);
  
  for (int i=1; i<=7; ++i) {
    for (int k=1; k<=6; ++k)  cofi(k) = MatCDWG(i,k);
    for (int j=1; j<=7; ++j) {
      for (int k=1; k<=6; ++k)  cofj(k) = MatCDWG(j,k);
      EGGMK(i,j) = dotProd(cofi, GMK*cofj);
    }
  }

  return EGGMK;
}


// WG3d: Hexa(P0,P0;AT0):
// 6 faces, 7 WG bas.fxns., 6 AT0 nmlz.bas.fxns., No projection for Kw
int WG3d_HexaP0P0AT0_NmlFluxK_NmlzPiolaBas_NoProj(FullMatrix &NmlFluxWG,
                                                  const Hexa &hexa, Mat3 &MatK,
                                                  Quadri3d SQF[6], int sign[6],
                                                  const GaussQuad &GQH,
                                                  const GaussQuad &GQQ)
{
  int i, j, k;
  double dp, jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat;
  Quadri3d quadri;
  PtVec3d cntr, nml, qp, w[7];  // w[0] not used
  FullMatrix CDWG(7,6), TransCDWG(6,7), NmlFluxRT0(6,6);
  
  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
  
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
      w[4] = PtVec3d(X,0,0);
      w[5] = PtVec3d(0,Y,0);
      w[6] = PtVec3d(0,0,Z);
      for (j=1; j<=6; ++j) {
        dp = - dotProduct(MatK*w[j], nml);  // Notice the negative sign!
        NmlFluxRT0(i,j) += GQQ.weight(k) * jac * dp;
      }
    }
    for (j=1; j<=6; ++j)
      NmlFluxRT0(i,j) = NmlFluxRT0(i,j) * sign[i-1];
  }
  
  // std::cout << NmlFluxRT0 << "\n";
  
  CDWG = WG3d_HexaP0P0AT0_CofDiscWkGrad_NmlzPiolaBas(hexa, GQH, GQQ);
  TransCDWG = transpose(CDWG);
  
  NmlFluxWG.resize(6,7);
  NmlFluxWG = NmlFluxRT0 * TransCDWG;
  
  return(0);  // if successful
}


////////////////////////////////////////////////////////////////////////////////
// For WG(Q0,Q0;RT[0])Hexa
////////////////////////////////////////////////////////////////////////////////


// WG3d: Hexa(Q0,Q0;RT[0]): Element mass matrix (1-by-1)
// Just a number: the volume of the hexahedron

FullMatrix WG3d_HexaQ0Q0RT0_EltMassMat(const Hexa &hexa, const GaussQuad &GQH)
{
  double jac, vol, xhat, yhat, zhat;
  FullMatrix EMM(1,1);
  vol = 0;
  for (int i=0; i<GQH.numberQuadraturePoints(); ++i) {
    xhat = GQH.CartesianCoordinate(i,0);
    yhat = GQH.CartesianCoordinate(i,1);
    zhat = GQH.CartesianCoordinate(i,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol = vol + GQH.weight(i)*jac;
  }
  EMM(1,1) = vol;
  return EMM;
}


// WG3d: Hexa(Q0,Q0;RT[0]): Coeff. of disc.wk.grad. in RT[0] nmlz.bas. (7-by-6)

FullMatrix WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(const Hexa &hexa,
                                                  const GaussQuad &GQH,
                                                  const GaussQuad &GQQ)
{
  int i, j, k;
  double jac, vol;
  double X, Y, Z;
  double x, y, z;
  double xc, yc, zc;
  double xhat, yhat, zhat;
  PtVec3d cntr, nml, qp, w[6];
  Quadri3d quadri;
  Mat3 MatK;
  FullMatrix GM(6,6), GMK(6,6), MatCDWG(7,6);
  Vector cof(6), RHS(6);

  //
  cntr = hexa.center();
  xc = cntr.xCrd();
  yc = cntr.yCrd();
  zc = cntr.zCrd();
  
  // JL20150522: THIS IS CORRECT BUT NEEDS TO BE REVISED
  MatK(1,1) = 1;
  MatK(2,2) = 1;
  MatK(3,3) = 1;

  //
  Hdiv_HexaRT0_GramMatK_NmlzBas(GM, GMK, hexa, MatK, GQH);
  /*
  std::cout << GM << "\n";
  std::cout << GMK << "\n";
  */

  // Computing the hexahedron volume
  vol = 0;
  for (k=0; k<GQH.numberQuadraturePoints(); ++k) {
    xhat = GQH.CartesianCoordinate(k,0);
    yhat = GQH.CartesianCoordinate(k,1);
    zhat = GQH.CartesianCoordinate(k,2);
    jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol = vol + GQH.weight(k) * jac;
  }

  // JL20150425: IT SEEMS NOT NEEDED.
  // Computing the areas of the 6 side quadrilaterals of the hexahedron
  /*
  for (i=0; i<6; ++i) {
    quadri = hexa.face(i);
    sa[i] = 0;
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      sa[i] = sa[i] + GQQ.weight(k) * jac;
    }
  }
  */

  // For the interior function
  RHS(1) = 0;  RHS(4) = -vol;
  RHS(2) = 0;  RHS(5) = -vol;
  RHS(3) = 0;  RHS(6) = -vol;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (j=1; j<=6; ++j)  MatCDWG(1,j) = cof(j);
  
  // For the boundary/face functions
  w[0] = PtVec3d(1,0,0);
  w[1] = PtVec3d(0,1,0);
  w[2] = PtVec3d(0,0,1);
  for (i=0; i<6; ++i) {  // i-th face
    // std::cout << "face(i)=" << i << "\n";
    quadri = hexa.face(i);
    /*
    for (j=0; j<4; ++j) {
      std::cout << quadri.vertex(j) << "\n";
    }
    std::cout << "\n";
    */
    for (j=1; j<=6; ++j)  RHS(j) = 0;  // j-th normalized basis of RT[0]
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      qp = quadri.bilinearmapping(xhat,yhat);
      nml = quadri.normal(xhat, yhat);
      // std::cout << "nml=" << nml << "\n";
      x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
      X = x - xc;  Y = y - yc;  Z = z - zc;
      w[3] = PtVec3d(X,0,0);
      w[4] = PtVec3d(0,Y,0);
      w[5] = PtVec3d(0,0,Z);
      for (j=1; j<=6; ++j) {
        RHS(j) = RHS(j) + GQQ.weight(k) * jac * dotProduct(w[j-1],nml);
      }
    }
    cof = slvFullSpdSysCholesky(GM, RHS);
    // std::cout << RHS << "\n";
    for (j=1; j<=6; ++j)  MatCDWG(i+2,j) = cof(j);
  }

  return MatCDWG;
}


// WG3d: Hexa(Q0,Q0;RT[0]): Element grad-grad matrix (7-by-7)
// with diffusion/permeability (3x3) matrix MatK

FullMatrix WG3d_HexaQ0Q0RT0_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ)
{
  int i, j, k;
  FullMatrix GM(6,6), GMK(6,6);
  FullMatrix MatCDWG(7,6);
  FullMatrix EGGMK(7,7);
  Vector cofi(6), cofj(6);

  // std::cout << MatK << "\n";
  Hdiv_HexaRT0_GramMatK_NmlzBas(GM, GMK, hexa, MatK, GQH);
  // std::cout << GM << "\n";
  // std::cout << GMK << "\n";
  MatCDWG = WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(hexa, GQH, GQQ);
  // std::cout << MatCDWG << "\n";

  for (i=1; i<=7; ++i) {
    for (k=1; k<=6; ++k)  cofi(k) = MatCDWG(i,k);
    for (j=1; j<=7; ++j) {
      for (k=1; k<=6; ++k)  cofj(k) = MatCDWG(j,k);
      EGGMK(i,j) = dotProd(cofi, GMK*cofj);
    }
  }

  return EGGMK;
}


// WG3d: Hexa(Q0,Q0;RT[0]):
// 6 faces, 7 WG bas.fxns., 6 RT[0] nmlz.bas.fxns., Yes projection for Kw
/*
int WG3d_HexaQ0Q0RT0_NmlFluxK_NmlzBas_Proj(FullMatrix &NmlFluxWG,
                                           const Hexa &hexa, Mat3 &MatK,
                                           Quadri3d SQF[6], int sign[6],
                                           const GaussQuad &GQH,
                                           const GaussQuad &GQQ)
{
  FullMatrix CDWG(7,6), TransCDWG(6,7);
  FullMatrix NmlFluxRT0(6,6), KtransMat(12,6), ProjMat(6,12);

  Hdiv_HexaRT0_NmlFlux_NmlzBas(NmlFluxRT0, hexa, SQF, sign, GQQ);
  Hdiv_Hexa01_KtransMat(KtransMat, MatK);
  Hdiv_Hexa01_ProjP1v3RT0(ProjMat, hexa, GQH);

  CDWG = WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(hexa, GQH, GQQ);
  TransCDWG = transpose(CDWG);

  NmlFluxWG = (NmlFluxRT0 * (ProjMat*KtransMat)) * TransCDWG;

  return(0);  // if successful
}
*/


// WG3d: Hexa(Q0,Q0;RT[0]):
// 6 faces, 7 WG bas.fxns., 6 RT[0] nmlz.bas.fxns., No projection for Kw

int WG3d_HexaQ0Q0RT0_NmlFluxK_NmlzBas_NoProj(FullMatrix &NmlFluxWG,
                                             const Hexa &hexa, Mat3 &MatK,
                                             Quadri3d SQF[6], int sign[6],
                                             const GaussQuad &GQH,
                                             const GaussQuad &GQQ)
{
  int i, j, k;
  double dp, jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat;
  Quadri3d quadri;
  PtVec3d cntr, nml, qp, w[7];  // w[0] not used
  FullMatrix CDWG(7,6), TransCDWG(6,7), NmlFluxRT0(6,6);

  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();

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
      w[4] = PtVec3d(X,0,0);
      w[5] = PtVec3d(0,Y,0);
      w[6] = PtVec3d(0,0,Z);
      for (j=1; j<=6; ++j) {
        dp = - dotProduct(MatK*w[j], nml);  // Notice the negative sign!
        NmlFluxRT0(i,j) += GQQ.weight(k) * jac * dp;
      }
    }
    for (j=1; j<=6; ++j)
      NmlFluxRT0(i,j) = NmlFluxRT0(i,j) * sign[i-1];
  }

  // std::cout << NmlFluxRT0 << "\n";
  
  CDWG = WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(hexa, GQH, GQQ);
  TransCDWG = transpose(CDWG);

  NmlFluxWG.resize(6,7);
  NmlFluxWG = NmlFluxRT0 * TransCDWG;

  // std::cout << NmlFluxRT0 << "\n";
  
  return(0);  // if successful
}


////////////////////////////////////////////////////////////////////////////////
// For WG(Q0,Q1;P1^3)Hexa
////////////////////////////////////////////////////////////////////////////////


// JL20150601: TO BE REVISED FOR EFFECIENCY
// WG3d: Hexa1(Q0,Q1;P1^3): Coeff.disc.wk.grad. in nmlz.bas.: matrix 12-by-25
// nmlz.bas.fxns as rows;  WG bas.fxns. as columns
/*
int WG3d_Hexa1_CofDiscWkGrad_NmlzBas(FullMatrix &CDWG, const Hexa &hexa,
                                     Quadri3d SQF[6], int sign[6],
                                     const GaussQuad &GQH, const GaussQuad &GQQ)
{
  int i, j, k, m;
  double jac, vol, x, y, z, xc, yc, zc, xhat, yhat, X, Y, Z;
  PtVec3d cntr, nml, qp, w[13];  // w[0] not used
  Quadri3d quadri;
  Vector cof(12), RHS(12), RHS1(12), RHS2(12), RHS3(12), RHS4(12);
  FullMatrix GM(12,12);
  CDWG.resize(12,25);

  //
  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();

  // Getting the Gram matrix of the normalized basis for P1^3
  Hdiv_Hexa1_GramMat_NmlzBas(GM, hexa, GQH);
  vol = GM(1,1);

  // For WG bas.fxn.(Q0): 1 constant function in tetrahedron interior
  // by default: RHS(1,2,3, 5,6, 7,9, 10,11) = 0
  RHS(4) = -vol;  RHS(8) = -vol;  RHS(12) = -vol;
  cof = slvFullSpdSysCholesky(GM, RHS);
  for (i=1; i<=12; ++i)  CDWG(i,1) = cof(i);

  // For WG bas.fxn.(Q1): 4*6=24 bilinear functions on 6 quadrilateral faces
  // w: 12 normalized (vector) basis functions for P1^3 on a hexahedron
  w[1] = PtVec3d(1,0,0);
  w[2] = PtVec3d(0,1,0);
  w[3] = PtVec3d(0,0,1);
  //
  for (m=0; m<6; ++m) {  // the m-th face, 4 bas.fxns. on each face
    quadri = SQF[m];

    // For 4 WG scalar basis functions on each face:
    //   (1-xhat)*(1-yhat), xhat*(1-yhat), xhat*yhat, (1-xhat)*yhat
    RHS1.resize(12);
    RHS2.resize(12);
    RHS3.resize(12);
    RHS4.resize(12);
    //
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      qp = quadri.bilinearmapping(xhat, yhat);
      nml = quadri.normal(xhat, yhat);
      x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
      X = x - xc;  Y = y - yc;  Z = z - zc;
      w[4] = PtVec3d(X,0,0);  w[5] = PtVec3d(Y,0,0);  w[6] = PtVec3d(Z,0,0);
      w[7] = PtVec3d(0,X,0);  w[8] = PtVec3d(0,Y,0);  w[9] = PtVec3d(0,Z,0);
      w[10]= PtVec3d(0,0,X);  w[11]= PtVec3d(0,0,Y);  w[12]= PtVec3d(0,0,Z);
      for (i=1; i<=12; ++i) {
        RHS1(i) += GQQ.weight(k) * jac * dotProduct(w[i],nml) * (1-xhat) * (1-yhat);
        RHS2(i) += GQQ.weight(k) * jac * dotProduct(w[i],nml) *     xhat * (1-yhat);
        RHS3(i) += GQQ.weight(k) * jac * dotProduct(w[i],nml) *     xhat * yhat;
        RHS4(i) += GQQ.weight(k) * jac * dotProduct(w[i],nml) * (1-xhat) * yhat;
      }
    }
    //
    RHS1 = sign[m] * RHS1;
    cof = slvFullSpdSysCholesky(GM, RHS1);
    for (i=1; i<=12; ++i)  CDWG(i,4*m+2) = cof(i);
    //
    RHS2 = sign[m] * RHS2;
    cof = slvFullSpdSysCholesky(GM, RHS2);
    for (i=1; i<=12; ++i)  CDWG(i,4*m+3) = cof(i);
    //
    RHS3 = sign[m] * RHS3;
    cof = slvFullSpdSysCholesky(GM, RHS3);
    for (i=1; i<=12; ++i)  CDWG(i,4*m+4) = cof(i);
    //
    RHS4 = sign[m] * RHS4;
    cof = slvFullSpdSysCholesky(GM, RHS4);
    for (i=1; i<=12; ++i)  CDWG(i,4*m+5) = cof(i);
  }
  
  return(0);  // if successful
}
*/


// WG3d: Hexa1(Q0,Q1;P1^3): Element grad-grad matrix w/ permeability MatK: 25*25
// (25=1+4*6) WG bas.fxns.: 1 const.fxn. interior, 4 Q1 bas.fxns. each face
// "Hdiv"(P1^3) normalized basis functions are used
/*
FullMatrix WG3d_Hexa1_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                      Quadri3d SQF[6], int sign[6],
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ)
{
  int i, j, k, l;
  FullMatrix CDWG(12,25), EGGMK(25,25), GMK(12,12), B;

  Hdiv_Hexa1_GramMatK_NmlzBas(GMK, hexa, MatK, GQH);
  // GMK.save2file("GMK.dat");

  WG3d_Hexa1_CofDiscWkGrad_NmlzBas(CDWG, hexa, SQF, sign, GQH, GQQ);
  // CDWG.save2file("CDWG.dat");
  // B = transpose(CDWG);
  // B.save2file("CDWG.dat");
  // std::cout << CDWG << "\n\n";

  // JL20150619: NEEDS DOUBLE-CHECK!!!
  for (i=1; i<=25; ++i) {
    for (j=1; j<=25; ++j) {
      EGGMK(i,j) = 0;
      for (k=1; k<=12; ++k)
        for (l=1; l<=12; ++l)
          EGGMK(i,j) += CDWG(k,i) * GMK(k,l) * CDWG(l,j);
    }
  }

  return EGGMK;
}
*/


// WG3d: Hexa1(Q0,Q1;P1^3): 6*25
// 6 faces, 25 WG bas.fxns., 12 P1^3 nmlz.bas.fxns., No projection needed for Kw
/*
FullMatrix WG3d_Hexa1_NmlFluxK_NmlzBas(const Hexa &hexa, Mat3 &MatK,
                                       Quadri3d SQF[6], int sign[6],
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ)
{
  int i, j, k;
  double dp, jac, X, Y, Z, x, y, z, xc, yc, zc, xhat, yhat;
  PtVec3d cntr, nml, qp, w[13];  // w[0] not used
  Quadri3d quadri;
  FullMatrix CDWG(12,25), NmlFluxP1v3(6,12), NmlFluxWG(6,25);

  cntr = hexa.center();
  xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
  
  w[1] = PtVec3d(1,0,0);
  w[2] = PtVec3d(0,1,0);
  w[3] = PtVec3d(0,0,1);

  for (i=1; i<=6; ++i) {  // the i-th face
    quadri = SQF[i-1];
    for (j=1; j<=12; ++j)  NmlFluxP1v3(i,j) = 0;
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      qp = quadri.bilinearmapping(xhat, yhat);
      nml = quadri.normal(xhat, yhat);
      x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
      X = x - xc;  Y = y - yc;  Z = z - zc;
      w[4] = PtVec3d(X,0,0);  w[5] = PtVec3d(Y,0,0);  w[6] = PtVec3d(Z,0,0);
      w[7] = PtVec3d(0,X,0);  w[8] = PtVec3d(0,Y,0);  w[9] = PtVec3d(0,Z,0);
      w[10]= PtVec3d(0,0,X);  w[11]= PtVec3d(0,0,Y);  w[12]= PtVec3d(0,0,Z);
      for (j=1; j<=12; ++j) {
        dp = - dotProduct(MatK*w[j], nml);  // Notice the negative sign!
        NmlFluxP1v3(i,j) += GQQ.weight(k) * jac * dp;
      }
    }
    for (j=1; j<=12; ++j)
      NmlFluxP1v3(i,j) = NmlFluxP1v3(i,j) * sign[i-1];
  }
  // std::cout << NmlFluxP1v3 << "\n";
  
  WG3d_Hexa1_CofDiscWkGrad_NmlzBas(CDWG, hexa, SQF, sign, GQH, GQQ);

  // NmlFluxWG.resize(6,25);
  NmlFluxWG = NmlFluxP1v3 * CDWG;

  return NmlFluxWG;
}
*/


// JL20160827: NEEDS DOUBLE-CHECK
////////////////////////////////////////////////////////////////////////////////
// For WG(P1,P0;P0^3)Hexa with stabilizer
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P1,P0;P0^3): Coeff. of disc.wk.grad. in P0^3 nat.bas. (10-by-3)

FullMatrix WG3d_HexaP1P0P03_CofDiscWkGrad_NatBas(const Hexa &hexa,
                                                 Quadri3d SQF[6], int sign[6],
                                                 const GaussQuad &GQH,
                                                 const GaussQuad &GQQ)
{
  FullMatrix CDWGB(10,3);

  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  double vol1 = 1/vol;

  // NOTE: For 1<=i<=4: coeff.disc.wk.grad.=0
  for (int i=5; i<=10; ++i) {
    Quadri3d quadri = SQF[i-5];
    PtVec3d vtmp(0,0,0);
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d nml = sign[i-5] * quadri.normal(xhat,yhat);
      vtmp = vtmp + (GQQ.weight(k)*jac) * nml;
    }
    vtmp = vol1 * vtmp;
    CDWGB(i,1) = vtmp.xCrd();
    CDWGB(i,2) = vtmp.yCrd();
    CDWGB(i,3) = vtmp.zCrd();
  }

  return CDWGB;
}


// WG3d: Hexa(P1,P0;P0^3): Element grad-grad matrix (10-by-10)

FullMatrix WG3d_HexaP1P0P03_EltGradGradMatK(const Hexa &hexa, Mat3 &MatK,
                                            Quadri3d SQF[6], int sign[6],
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ)
{
  FullMatrix EGGMK(10,10);
  FullMatrix CDWGB(10,3);
  FullMatrix GMK(3,3);
  
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }

  for (int i=1; i<=3; ++i)
    for (int j=1; j<=3; ++j)
      GMK(i,j) = vol * MatK(i,j);

  CDWGB = WG3d_HexaP1P0P03_CofDiscWkGrad_NatBas(hexa, SQF, sign, GQH, GQQ);

  // JL20160827: TO BE REVISED FOR EFFECIENCY
  // NOTE: i,j start at 5, because coeff.disc.wk.grad.=0 for 1<=i<=4
  for (int i=5; i<=10; ++i) {
    for (int j=5; j<=10; ++j) {
      EGGMK(i,j) = 0;
      for (int k=1; k<=3; ++k) {
        for (int l=1; l<=3; ++l) {
          EGGMK(i,j) += GMK(k,l) * CDWGB(i,k) * CDWGB(j,l);
        }
      }
    }
  }

  return EGGMK;
}


// JL20160827: TO BE REVISED FOR BETTER ACCURACY (better than midpoint values)
// WG3d: Hexa(P1,P0;P0^3): Elementwise stabilizer matrix (10-by-10)

FullMatrix WG3d_HexaP1P0P03_EltStabMat(const Hexa &hexa,
                                       Quadri3d SQF[6], int sign[6],
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ)
{
  FullMatrix ESM(10,10);
  double ar[6];  // ar for quadri.surf. area
  double Xm[6], Ym[6], Zm[6];  // midpoint coordinates for X,Y,Z on 6 faces
  
  PtVec3d cntr = hexa.center();
  double xc = cntr.xCrd();
  double yc = cntr.yCrd();
  double zc = cntr.zCrd();
  
  for (int i=1; i<=6; ++i) {
    PtVec3d midpt = SQF[i-1].center();
    Xm[i-1] = midpt.xCrd() - xc;
    Ym[i-1] = midpt.yCrd() - yc;
    Zm[i-1] = midpt.zCrd() - zc;
  }

  for (int i=1; i<=6; ++i) {
    Quadri3d quadri = SQF[i-1];
    ar[i-1] = 0;
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      ar[i-1] += GQQ.weight(k) * jac;
    }
  }

  // For element interior basis functions themselves
  for (int j=1; j<=6; ++j) {
    ESM(1,1) += ar[j-1];
    ESM(1,2) += ar[j-1]*Xm[j-1];
    ESM(1,3) += ar[j-1]*Ym[j-1];
    ESM(1,4) += ar[j-1]*Zm[j-1];
    ESM(2,2) += ar[j-1]*Xm[j-1]*Xm[j-1];
    ESM(2,3) += ar[j-1]*Ym[j-1]*Xm[j-1];
    ESM(2,4) += ar[j-1]*Zm[j-1]*Xm[j-1];
    ESM(3,3) += ar[j-1]*Ym[j-1]*Ym[j-1];
    ESM(3,4) += ar[j-1]*Zm[j-1]*Ym[j-1];
    ESM(4,4) += ar[j-1]*Zm[j-1]*Zm[j-1];
  }
  for (int i=2; i<=4; ++i)
    for (int j=1; j<=i-1; ++j)
      ESM(i,j) = ESM(j,i);

  // For interaction of interior bas.fxns. and face bas.fxns.
  for (int j=1; j<=6; ++j) {
    ESM(1,4+j) = -ar[j-1];          ESM(4+j,1) = ESM(1,4+j);
    ESM(2,4+j) = -ar[j-1]*Xm[j-1];  ESM(4+j,2) = ESM(2,4+j);
    ESM(3,4+j) = -ar[j-1]*Ym[j-1];  ESM(4+j,3) = ESM(3,4+j);
    ESM(4,4+j) = -ar[j-1]*Zm[j-1];  ESM(4+j,4) = ESM(4,4+j);
  }

  // For the face basis functions themselves
  for (int i=1; i<=6; ++i)  ESM(4+i,4+i) = ar[i-1];

  return ESM;
}


////////////////////////////////////////////////////////////////////////////////
// For WG(P0^3,P0^3;AT0^3,P0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Coeff. of disc.wk.grad. in AT0^3 nmlz.Piola.bas. for 21 WG bas.fxns
// 21-by-18 matrix

FullMatrix WG3d_HexaP03P03AT03P0_DiscWkGradBasFxn_NmlzPiolaBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                               const GaussQuad &GQH, const GaussQuad &GQQ)
{
  FullMatrix GM(18,18);
  Hdiv_HexaAT03_GramMat_NmlzPiolaBas(GM, hexa, GQH);
  
  // Hexa. volume needed
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  PtVec3d cntr = hexa.center();
  
  // Auxiliary
  FullMatrix AllRHS(21,18);
  AllRHS(1,4) = -vol;  AllRHS(2,10) = -vol;  AllRHS(3,16) = -vol;
  AllRHS(1,5) = -vol;  AllRHS(2,11) = -vol;  AllRHS(3,17) = -vol;
  AllRHS(1,6) = -vol;  AllRHS(2,12) = -vol;  AllRHS(3,18) = -vol;
  for (int i=1; i<=6; ++i) {
    Quadri3d quadri = SQF[i-1];
    
    // JL20170505: Need side surface area !
    double ar = 0;  // For face area
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      ar += GQQ.weight(k) * jac;
    }
    ar = ar * sign[i-1];
    
    // Face center
    PtVec3d FcCntr = quadri.center();
    double Xm = FcCntr.xCrd() - cntr.xCrd();
    double Ym = FcCntr.yCrd() - cntr.yCrd();
    double Zm = FcCntr.zCrd() - cntr.zCrd();
    
    // JL20170505: TO BE CHANGED TO averaged normal !!
    PtVec3d nml = quadri.normal(0.5, 0.5);  // normal at center
    
    // All right-hand-side
    AllRHS(3*i+1, 1) = ar * nml.xCrd();
    AllRHS(3*i+1, 2) = ar * nml.yCrd();
    AllRHS(3*i+1, 3) = ar * nml.zCrd();
    AllRHS(3*i+1, 4) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+1, 5) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+1, 6) = (Zm * ar) * nml.zCrd();
    //
    AllRHS(3*i+2, 7) = ar * nml.xCrd();
    AllRHS(3*i+2, 8) = ar * nml.yCrd();
    AllRHS(3*i+2, 9) = ar * nml.zCrd();
    AllRHS(3*i+2,10) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+2,11) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+2,12) = (Zm * ar) * nml.zCrd();
    //
    AllRHS(3*i+3,13) = ar * nml.xCrd();
    AllRHS(3*i+3,14) = ar * nml.yCrd();
    AllRHS(3*i+3,15) = ar * nml.zCrd();
    AllRHS(3*i+3,16) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+3,17) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+3,18) = (Zm * ar) * nml.zCrd();
  }
  
  // Now solving size-18 lin.sys. for coeffs.
  FullMatrix CDWGB(21,18);
  Vector cof(18);
  Vector RHS(18);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=18; ++j)  RHS(j) = AllRHS(i,j);
    cof = slvFullSpdSysCholesky(GM, RHS);
    for (int j=1; j<=18; ++j)  CDWGB(i,j) = cof(j);
  }
  
  return CDWGB;
}


// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Coeff. of disc.wk.div. in P0 nmlz.bas. for 21 WG bas.fxns
// 21-dim vector

Vector WG3d_HexaP03P03AT03P0_DiscWkDivBasFxn_NmlzBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                     const GaussQuad &GQH, const GaussQuad &GQQ)
{
  Vector CDWDB(21);

  // Need hexa. volume
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  double vol1 = 1/vol;

  // NOTE: For the interior 3 bas.fxns, disc.wk.div.=0
  for (int i=0; i<6; ++i) {
    Quadri3d quadri = SQF[i];
    PtVec3d vtmp(0,0,0);
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d nml = quadri.normal(xhat,yhat);
      vtmp = vtmp + (GQQ.weight(k)*jac) * nml;
    }
    vtmp = (sign[i] * vol1) * vtmp;

    CDWDB(3+3*i+1) = vtmp.xCrd();
    CDWDB(3+3*i+2) = vtmp.yCrd();
    CDWDB(3+3*i+3) = vtmp.zCrd();
  }

  return CDWDB;
}


// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Elementwise discrete weak strain-strain matrix

FullMatrix WG3d_HexaP03P03AT03P0_EltStrnStrnMat_NmlzPiolaBas(const Hexa &hexa,
                                                             Quadri3d SQF[6], int sign[6],
                                                             const GaussQuad &GQH, const GaussQuad &GQQ)
{
  FullMatrix CDWGB(21,18);
  FullMatrix GMA(18,18);
  CDWGB = WG3d_HexaP03P03AT03P0_DiscWkGradBasFxn_NmlzPiolaBas(hexa, SQF, sign, GQH, GQQ);
  Hdiv_HexaAT03_GramMatAvg_NmlzPiolaBas(GMA, hexa, GQH);
  FullMatrix ESSM(21,21);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=21; ++j) {
      ESSM(i,j) = 0;
      for (int k=1; k<=18; ++k) {
        for (int l=1; l<=18; ++l) {
          ESSM(i,j) = ESSM(i,j) + CDWGB(i,k) * GMA(k,l) * CDWGB(j,l);
        }
      }
    }
  }
  return ESSM;
}


// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Elementwise discrete weak div-div matrix

FullMatrix WG3d_HexaP03P03AT03P0_EltDivDivMat_NmlzBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                      const GaussQuad &GQH, const GaussQuad &GQQ)
{
  Vector CDWDB(21);
  CDWDB = WG3d_HexaP03P03AT03P0_DiscWkDivBasFxn_NmlzBas(hexa, SQF, sign, GQH, GQQ);
  // Need hexa. volume
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  FullMatrix EDDM(21,21);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=21; ++j) {
      EDDM(i,j) = CDWDB(i) * CDWDB(j) * vol;
    }
  }
  return EDDM;
}


////////////////////////////////////////////////////////////////////////////////
// For WG(Q0^3,Q0^3;RT[0]^3,Q0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Coeff. of disc.wk.grad. in RT[0]^3 nmlz.bas. for 21 WG bas.fxns 
// 21-by-18 matrix

FullMatrix WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkGradBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  FullMatrix GM(18,18);
  Hdiv_HexaRT03_GramMat_NmlzBas(GM, hexa, GQH);

  // Hexa. volume needed
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  PtVec3d cntr = hexa.center();
 
  // Auxiliary
  FullMatrix AllRHS(21,18);
  AllRHS(1,4) = -vol;  AllRHS(2,10) = -vol;  AllRHS(3,16) = -vol;
  AllRHS(1,5) = -vol;  AllRHS(2,11) = -vol;  AllRHS(3,17) = -vol;
  AllRHS(1,6) = -vol;  AllRHS(2,12) = -vol;  AllRHS(3,18) = -vol;
  for (int i=1; i<=6; ++i) {
    Quadri3d quadri = SQF[i-1];

    // JL20170505: Need side surface area !
    double ar = 0;  // For face area
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      ar += GQQ.weight(k) * jac;
    }
    ar = ar * sign[i-1];

    // Face center
    PtVec3d FcCntr = quadri.center();
    double Xm = FcCntr.xCrd() - cntr.xCrd();
    double Ym = FcCntr.yCrd() - cntr.yCrd();
    double Zm = FcCntr.zCrd() - cntr.zCrd();

    // JL20170505: TO BE CHANGED TO averaged normal !!
    PtVec3d nml = quadri.normal(0.5, 0.5);  // normal at center

    // All right-hand-side
    AllRHS(3*i+1, 1) = ar * nml.xCrd();
    AllRHS(3*i+1, 2) = ar * nml.yCrd();
    AllRHS(3*i+1, 3) = ar * nml.zCrd();
    AllRHS(3*i+1, 4) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+1, 5) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+1, 6) = (Zm * ar) * nml.zCrd();
    //
    AllRHS(3*i+2, 7) = ar * nml.xCrd();
    AllRHS(3*i+2, 8) = ar * nml.yCrd();
    AllRHS(3*i+2, 9) = ar * nml.zCrd();
    AllRHS(3*i+2,10) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+2,11) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+2,12) = (Zm * ar) * nml.zCrd();
    //
    AllRHS(3*i+3,13) = ar * nml.xCrd();
    AllRHS(3*i+3,14) = ar * nml.yCrd();
    AllRHS(3*i+3,15) = ar * nml.zCrd();
    AllRHS(3*i+3,16) = (Xm * ar) * nml.xCrd();
    AllRHS(3*i+3,17) = (Ym * ar) * nml.yCrd();
    AllRHS(3*i+3,18) = (Zm * ar) * nml.zCrd();
  }

  // Now solving size-18 lin.sys. for coeffs.
  FullMatrix CDWGB(21,18);
  Vector cof(18);
  Vector RHS(18);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=18; ++j)  RHS(j) = AllRHS(i,j);
    cof = slvFullSpdSysCholesky(GM, RHS);
    for (int j=1; j<=18; ++j)  CDWGB(i,j) = cof(j);
  }

  return CDWGB;
}


// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Coeff. of disc.wk.div. in Q0 nmlz.bas. for 21 WG bas.fxns
// 21-dim vector

Vector WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkDivBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  Vector CDWDB(21);

  // Need hexa. volume
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  double vol1 = 1/vol;

  // NOTE: For the interior 3 bas.fxns, disc.wk.div.=0
  for (int i=0; i<6; ++i) {
    Quadri3d quadri = SQF[i];
    PtVec3d vtmp(0,0,0);
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d nml = quadri.normal(xhat,yhat);
      vtmp = vtmp + (GQQ.weight(k)*jac) * nml;
    }
    vtmp = (sign[i] * vol1) * vtmp;
    CDWDB(3+3*i+1) = vtmp.xCrd();
    CDWDB(3+3*i+2) = vtmp.yCrd();
    CDWDB(3+3*i+3) = vtmp.zCrd();
  }
  
  return CDWDB;
}


// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Elementwise discrete weak starin-strain matrix

FullMatrix WG3d_HexaQ03Q03RT03Q0_NmlBas_EltStrnStrnMat(const Hexa &hexa,
                                                       Quadri3d SQF[6], int sign[6],
                                                       const GaussQuad &GQH,
                                                       const GaussQuad &GQQ)
{
  FullMatrix CDWGB(21,18);
  CDWGB = WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkGradBasFxn(hexa, SQF, sign, GQH, GQQ);

  FullMatrix GMA(18,18);
  Hdiv_HexaRT03_GramMatAvg_NmlzBas(GMA, hexa, GQH);

  FullMatrix ESSM(21,21);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=21; ++j) {
      ESSM(i,j) = 0;
      for (int k=1; k<=18; ++k) {
        for (int l=1; l<=18; ++l) {
          ESSM(i,j) = ESSM(i,j) + CDWGB(i,k) * GMA(k,l) * CDWGB(j,l);
        }
      }
    }
  }

  return ESSM;
}


// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Elementwise discrete weak div-div matrix

FullMatrix WG3d_HexaQ03Q03RT03Q0_NmlBas_EltDivDivMat(const Hexa &hexa,
                                                     Quadri3d SQF[6], int sign[6],
                                                     const GaussQuad &GQH,
                                                     const GaussQuad &GQQ)
{
  Vector CDWDB(21);
  CDWDB = WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkDivBasFxn(hexa, SQF, sign, GQH, GQQ);

  // Need hexa. volume
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  
  FullMatrix EDDM(21,21);
  for (int i=1; i<=21; ++i) {
    for (int j=1; j<=21; ++j) {
      EDDM(i,j) = CDWDB(i) * CDWDB(j) * vol;
    }
  }

  return EDDM;
}


////////////////////////////////////////////////////////////////////////////////
// For WG(P1^3,Prm;P0^{3x3},P0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Coeff. of disc.wk.grad. in P0^{3x3} nmlz.bas. for 48 WG bas.fxns
// 48-by-9 matrix

FullMatrix WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkGradBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  // Gram mat. for P0^3 is a scalar matrix with entry = volume;
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  double vol1 = 1/vol;

  FullMatrix CDWGB(48,9);
  FullMatrix FaceRHS(6,9);
  for (int i=1; i<=6; ++i) {
    FaceRHS.resize(6,9);
    Quadri3d quadri = SQF[i-1];
    PtVec3d cntr = quadri.center();
    double xm = cntr.xCrd();
    double ym = cntr.yCrd();
    double zm = cntr.zCrd();
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d nml = quadri.normal(xhat,yhat);
      PtVec3d qp = quadri.bilinearmapping(xhat,yhat);
      double xtilde = qp.xCrd() - xm;
      double ytilde = qp.yCrd() - ym;
      double ztilde = qp.zCrd() - zm;
      FaceRHS(1,1) = FaceRHS(1,1) + GQQ.weight(k) * jac * nml.xCrd();
      FaceRHS(1,2) = FaceRHS(1,2) + GQQ.weight(k) * jac * nml.yCrd();
      FaceRHS(1,3) = FaceRHS(1,3) + GQQ.weight(k) * jac * nml.zCrd();
      //
      FaceRHS(2,4) = FaceRHS(2,4) + GQQ.weight(k) * jac * nml.xCrd();
      FaceRHS(2,5) = FaceRHS(2,5) + GQQ.weight(k) * jac * nml.yCrd();
      FaceRHS(2,6) = FaceRHS(2,6) + GQQ.weight(k) * jac * nml.zCrd();
      //
      FaceRHS(3,7) = FaceRHS(3,7) + GQQ.weight(k) * jac * nml.xCrd();
      FaceRHS(3,8) = FaceRHS(3,8) + GQQ.weight(k) * jac * nml.yCrd();
      FaceRHS(3,9) = FaceRHS(3,9) + GQQ.weight(k) * jac * nml.zCrd();
      //
      FaceRHS(4,4) = FaceRHS(4,4) + GQQ.weight(k) * jac * nml.xCrd() *   ztilde ;
      FaceRHS(4,5) = FaceRHS(4,5) + GQQ.weight(k) * jac * nml.yCrd() *   ztilde ;
      FaceRHS(4,6) = FaceRHS(4,6) + GQQ.weight(k) * jac * nml.zCrd() *   ztilde ;
      FaceRHS(4,7) = FaceRHS(4,7) + GQQ.weight(k) * jac * nml.xCrd() * (-ytilde);
      FaceRHS(4,8) = FaceRHS(4,8) + GQQ.weight(k) * jac * nml.yCrd() * (-ytilde);
      FaceRHS(4,9) = FaceRHS(4,9) + GQQ.weight(k) * jac * nml.zCrd() * (-ytilde);
      //
      FaceRHS(5,1) = FaceRHS(5,1) + GQQ.weight(k) * jac * nml.xCrd() * (-ztilde);
      FaceRHS(5,2) = FaceRHS(5,2) + GQQ.weight(k) * jac * nml.yCrd() * (-ztilde);
      FaceRHS(5,3) = FaceRHS(5,3) + GQQ.weight(k) * jac * nml.zCrd() * (-ztilde);
      FaceRHS(5,7) = FaceRHS(5,7) + GQQ.weight(k) * jac * nml.xCrd() *   xtilde ;
      FaceRHS(5,8) = FaceRHS(5,8) + GQQ.weight(k) * jac * nml.yCrd() *   xtilde ;
      FaceRHS(5,9) = FaceRHS(5,9) + GQQ.weight(k) * jac * nml.zCrd() *   xtilde ;
      //
      FaceRHS(6,1) = FaceRHS(6,1) + GQQ.weight(k) * jac * nml.xCrd() *   ytilde ;
      FaceRHS(6,2) = FaceRHS(6,2) + GQQ.weight(k) * jac * nml.yCrd() *   ytilde ;
      FaceRHS(6,3) = FaceRHS(6,3) + GQQ.weight(k) * jac * nml.zCrd() *   ytilde ;
      FaceRHS(6,4) = FaceRHS(6,4) + GQQ.weight(k) * jac * nml.xCrd() * (-xtilde);
      FaceRHS(6,5) = FaceRHS(6,5) + GQQ.weight(k) * jac * nml.yCrd() * (-xtilde);
      FaceRHS(6,6) = FaceRHS(6,6) + GQQ.weight(k) * jac * nml.zCrd() * (-xtilde);
    }
    FaceRHS = sign[i-1] * FaceRHS;
    for (int ii=1; ii<=6; ++ii) {
      for (int j=1; j<=9; ++j) {
        CDWGB(12+6*(i-1)+ii,j) = FaceRHS(ii,j)*vol1;
      }
    }
  }
  
  return CDWGB;
}


// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Coeff. of disc.wk.div. in P0 nmlz.bas. for 48 WG bas.fxns
// 48-dim vector

Vector WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkDivBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  // Gram mat. for P0^3 is a scalar matrix with entry = volume;
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  double vol1 = 1/vol;

  Vector CDWDB(48);
  Vector RHS(6);
  for (int i=1; i<=6; ++i) {
    RHS.resize(6);
    Quadri3d quadri = SQF[i-1];
    PtVec3d cntr = quadri.center();
    double xm = cntr.xCrd();
    double ym = cntr.yCrd();
    double zm = cntr.zCrd();
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d nml = quadri.normal(xhat,yhat);
      double nml1 = nml.xCrd();
      double nml2 = nml.yCrd();
      double nml3 = nml.zCrd();
      PtVec3d qp = quadri.bilinearmapping(xhat,yhat);
      double xtilde = qp.xCrd() - xm;
      double ytilde = qp.yCrd() - ym;
      double ztilde = qp.zCrd() - zm;
      RHS(1) = RHS(1) + GQQ.weight(k) * jac * nml.xCrd();
      RHS(2) = RHS(2) + GQQ.weight(k) * jac * nml.yCrd();
      RHS(3) = RHS(3) + GQQ.weight(k) * jac * nml.zCrd();
      RHS(4) = RHS(4) + GQQ.weight(k) * jac * (nml2 * ztilde - nml3 * ytilde);
      RHS(5) = RHS(5) + GQQ.weight(k) * jac * (nml3 * xtilde - nml1 * ztilde);
      RHS(6) = RHS(6) + GQQ.weight(k) * jac * (nml1 * ytilde - nml2 * xtilde);
    }
    RHS = sign[i-1] * RHS;
    for (int ii=1; ii<=6; ++ii)  CDWDB(12+6*(i-1)+ii) = RHS(ii)*vol1;
  }

  return CDWDB;
}


// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise discrete weak strain-strain matrix
// 48-by-48 matrix

FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltStrnStrnMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  FullMatrix CDWGB(48,9);
  CDWGB = WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkGradBasFxn(hexa,SQF,sign,GQH,GQQ);
  
  FullMatrix GMA(9,9);
  Hdiv_HexaP033_GramMatAvg_NmlzBas(GMA, hexa, GQH);

  FullMatrix ESSM(48,48);
  for (int i=1; i<=48; ++i) {
    for (int j=1; j<=48; ++j) {
      ESSM(i,j) = 0;
      for (int k=1; k<=9; ++k) {
        for (int l=1; l<=9; ++l) {
          ESSM(i,j) = ESSM(i,j) + CDWGB(i,k) * GMA(k,l) * CDWGB(j,l);
        }
      }
    }
  }

  return ESSM;
}


// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise discrete weak div-div matrix
// 48-by-48 matrix

FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltDivDivMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  Vector CDWDB(48);
  CDWDB = WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkDivBasFxn(hexa,SQF,sign,GQH,GQQ);
  
  // Need hexa. volume
  // Computing the hexahedron volume
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol += GQH.weight(k) * jac;
  }
  
  FullMatrix EDDM(48,48);
  for (int i=1; i<=48; ++i) {
    for (int j=1; j<=48; ++j) {
      EDDM(i,j) = CDWDB(i) * CDWDB(j) * vol;
    }
  }

  return EDDM;
}


// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise stabilization matrix
// 48-by-48 matrix

FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltStabMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ)
{
  int i0, j0;
  Vector cof(6), SingleRHS(6);
  FullMatrix ESM(48,48), PrmGM(6,6), ProjCof(12,6), RHS(12,6), SubMat;

  // std::cout << "In WG3d_HexaP13PrmP033P0_NmlzBas_EltStabMat ..." << "\n";
  
  // Hexahedral center and diameter
  PtVec3d cntr = hexa.center();
  double xc = cntr.xCrd();
  double yc = cntr.yCrd();
  double zc = cntr.zCrd();
  double diam = hexa.diameter();

  for (int j=1; j<=6; ++j) {
    Quadri3d quadri = SQF[j-1];

    // Computing quadrilateral face area and center
    double ar = 0;  // For face area
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      ar += GQQ.weight(k) * jac;
    }
    ar = fabs(ar);
    // std::cout << "ar=" << ar << "\n";
    
    PtVec3d FcCntr = quadri.center();
    double xm = FcCntr.xCrd();
    double ym = FcCntr.yCrd();
    double zm = FcCntr.zCrd();
    
    // Computing auxiliary quantities for Gram matrix of face bas.fxns.
    double S1=0,  Sx=0,  Sy=0,  Sz=0;
    double Sx2=0,  Sy2=0,  Sz2=0,  Sxy=0,  Sxz=0,  Syz=0;
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d qp = quadri.bilinearmapping(xhat,yhat);
      double xtilde = qp.xCrd() - xm;
      double ytilde = qp.yCrd() - ym;
      double ztilde = qp.zCrd() - zm;
      S1 += GQQ.weight(k)*jac;
      Sx += GQQ.weight(k)*jac * xtilde;
      Sy += GQQ.weight(k)*jac * ytilde;
      Sz += GQQ.weight(k)*jac * ztilde;
      Sx2 += GQQ.weight(k)*jac * (xtilde*xtilde);
      Sy2 += GQQ.weight(k)*jac * (ytilde*ytilde);
      Sz2 += GQQ.weight(k)*jac * (ztilde*ztilde);
      Sxy += GQQ.weight(k)*jac * (xtilde*ytilde);
      Sxz += GQQ.weight(k)*jac * (xtilde*ztilde);
      Syz += GQQ.weight(k)*jac * (ytilde*ztilde);
    }
    // Setting up Gram matrix of face bas.fxns.
    PrmGM.resize(6,6);
    //
    PrmGM(1,1) = S1;  PrmGM(2,2) = S1;  PrmGM(3,3) = S1;
    PrmGM(4,4) = Sy2 + Sz2;
    PrmGM(5,5) = Sx2 + Sz2;
    PrmGM(6,6) = Sx2 + Sy2;
    PrmGM(4,5) = -Sxy;  PrmGM(4,6) = -Sxz;  PrmGM(5,6) = -Syz;
    PrmGM(5,4) = -Sxy;  PrmGM(6,4) = -Sxz;  PrmGM(6,5) = -Syz;
    PrmGM(1,5) = -Sz;  PrmGM(1,6) =  Sy;  PrmGM(2,6) = -Sx;
    PrmGM(5,1) = -Sz;  PrmGM(6,1) =  Sy;  PrmGM(6,2) = -Sx;
    PrmGM(2,4) =  Sz;  PrmGM(3,4) = -Sy;  PrmGM(3,5) =  Sx;
    PrmGM(4,2) =  Sz;  PrmGM(4,3) = -Sy;  PrmGM(5,3) =  Sx;
    // std::cout << PrmGM << "\n";
    //
    // Computing auxiliary quantities for FaceRHS for projection
    // Note S1, Sx, Sy, Sz already calculated
    double SX=0,  SXx=0,  SXy=0,  SXz=0;
    double SY=0,  SYx=0,  SYy=0,  SYz=0;
    double SZ=0,  SZx=0,  SZy=0,  SZz=0;
    for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      double xhat = GQQ.CartesianCoordinate(k,0);
      double yhat = GQQ.CartesianCoordinate(k,1);
      double jac = quadri.JacobianDeterminant(xhat, yhat);
      jac = fabs(jac);  // JL20150425: REALLY NEEDED?
      PtVec3d qp = quadri.bilinearmapping(xhat,yhat);
      double xtilde = qp.xCrd() - xm;
      double ytilde = qp.yCrd() - ym;
      double ztilde = qp.zCrd() - zm;
      double X = qp.xCrd() - xc;
      double Y = qp.yCrd() - yc;
      double Z = qp.zCrd() - zc;
      SX += GQQ.weight(k)*jac * X;
      SY += GQQ.weight(k)*jac * Y;
      SZ += GQQ.weight(k)*jac * Z;
      SXx += GQQ.weight(k)*jac * (X*xtilde);
      SYx += GQQ.weight(k)*jac * (Y*xtilde);
      SZx += GQQ.weight(k)*jac * (Z*xtilde);
      SXy += GQQ.weight(k)*jac * (X*ytilde);
      SYy += GQQ.weight(k)*jac * (Y*ytilde);
      SZy += GQQ.weight(k)*jac * (Z*ytilde);
      SXz += GQQ.weight(k)*jac * (X*ztilde);
      SYz += GQQ.weight(k)*jac * (Y*ztilde);
      SZz += GQQ.weight(k)*jac * (Z*ztilde);
    }
    // Setting up FaceRHS
    RHS.resize(12,6);
    //
    RHS(1,1) = S1;  RHS(1,5) = -Sz;   RHS(1,6) = Sy;
    RHS(2,1) = SX;  RHS(2,5) = -SXz;  RHS(2,6) = SXy;
    RHS(3,1) = SY;  RHS(3,5) = -SYz;  RHS(3,6) = SYy;
    RHS(4,1) = SZ;  RHS(4,5) = -SZz;  RHS(4,6) = SZy;
    //
    RHS(5,2) = S1;  RHS(5,4) = Sz;   RHS(5,6) = -Sx;
    RHS(6,2) = SX;  RHS(6,4) = SXz;  RHS(6,6) = -SXx;
    RHS(7,2) = SY;  RHS(7,4) = SYz;  RHS(7,6) = -SYx;
    RHS(8,2) = SZ;  RHS(8,4) = SZz;  RHS(8,6) = -SZx;
    //
    RHS( 9,3) = S1;  RHS( 9,4) = -Sy;   RHS( 9,5) = Sx;
    RHS(10,3) = SX;  RHS(10,4) = -SXy;  RHS(10,5) = SXx;
    RHS(11,3) = SY;  RHS(11,4) = -SYy;  RHS(11,5) = SYx;
    RHS(12,3) = SZ;  RHS(12,4) = -SZy;  RHS(12,5) = SZx;
    //
    // RHS = sign[j-1]*RHS;
    // Solving for projection coefficients
    for (int i=1; i<=12; ++i) {
      for (int k=1; k<=6; ++k)  SingleRHS(k) = RHS(i,k);
      cof = slvFullSpdSysCholesky(PrmGM, SingleRHS);
      for (int k=1; k<=6; ++k)  ProjCof(i,k) = cof(k);
    }
    // std::cout << "ProjCof=" << ProjCof << "\n";
    // For interior-interior interaction
    SubMat = (ar/diam) * (ProjCof * (PrmGM*transpose(ProjCof)));
    // std::cout << "SubMat=" << SubMat << "\n";
    
    for (int ii=1; ii<=12; ++ii) {
      for (int jj=1; jj<=12; ++jj) {
        ESM(ii,jj) = ESM(ii,jj) + SubMat(ii,jj);
      }
    }
    // std::cout << "ESM=" << ESM << "\n";
    // For interior-face-interior interaction
    SubMat = -(ar/diam) * ProjCof*PrmGM;
    i0 = 0;
    j0 = 12 + 6*(j-1);
    for (int ii=1; ii<=12; ++ii) {
      for (int jj=1; jj<=6; ++jj) {
        ESM(i0+ii,j0+jj) = SubMat(ii,jj);
        ESM(j0+jj,i0+ii) = SubMat(ii,jj);
      }
    }
    // For face itself
    SubMat = (ar/diam) * PrmGM;
    i0 = 12 + 6*(j-1);
    j0 = 12 + 6*(j-1);
    for (int ii=1; ii<=6; ++ii) {
      for (int jj=1; jj<=6; ++jj) {
        ESM(i0+ii,j0+jj) = SubMat(ii,jj);
      }
    }
  }

  return ESM;
}

// WG3d_Hexa.cpp
