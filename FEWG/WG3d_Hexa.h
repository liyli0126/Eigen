// WG3d_Hexa.h
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

#ifndef WG3D_HEXA_H
#define WG3D_HEXA_H

#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "mat3.h"


////////////////////////////////////////////////////////////////////////////////
// For WG(P0,P0;AT0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P0,P0;AT0): Element mass matrix (1-by-1, or just a number)
FullMatrix WG3d_HexaP0P0AT0_EltMassMat(const Hexa &hexa, const GaussQuad &GQH);

// WG3d: Hexa(P0,P0;AT0): Coeffs. of disc.wk.grad. in AT0 nmlz.bas. (7-by-6)
FullMatrix WG3d_HexaP0P0AT0_CofDiscWkGrad_NmlzPiolaBas(const Hexa &hexa,
                                                       const GaussQuad &GQH,
                                                       const GaussQuad &GQQ);

// WG3d: Hexa(P0,P0;AT0): Element grad-grad matrix (7-by-7)
// with diffusion/permeability (3x3) matrix MatK
FullMatrix WG3d_HexaP0P0AT0_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ);

// WG3d: Hexa(P0,P0;AT0):
// 6 faces, 7 WG bas.fxns., 6 AT0 nmlz.bas.fxns., No projection for Kw
int WG3d_HexaP0P0AT0_NmlFluxK_NmlzPiolaBas_NoProj(FullMatrix &NmlFluxWG,
                                                  const Hexa &hexa, Mat3 &MatK,
                                                  Quadri3d SQF[6], int sign[6],
                                                  const GaussQuad &GQH,
                                                  const GaussQuad &GQQ);


////////////////////////////////////////////////////////////////////////////////
// For WG(Q0,Q0;RT[0])Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(Q0,Q0;RT[0]): Element mass matrix (1-by-1, or just a number)
FullMatrix WG3d_HexaQ0Q0RT0_EltMassMat(const Hexa &hexa, const GaussQuad &GQH);

// WG3d: Hexa(Q0,Q0;RT[0]): Coeffs. of disc.wk.grad. in RT[0] nmlz.bas. (7-by-6)
FullMatrix WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(const Hexa &hexa,
                                                  const GaussQuad &GQH,
                                                  const GaussQuad &GQQ);

// WG3d: Hexa(Q0,Q0;RT[0]): Element grad-grad matrix (7-by-7)
// with diffusion/permeability (3x3) matrix MatK
FullMatrix WG3d_HexaQ0Q0RT0_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ);

// WG3d: Hexa(Q0,Q0;RT[0]):
// 6 faces, 7 WG bas.fxns., 6 RT[0] nmlz.bas.fxns., Yes projection for Kw
/*
int WG3d_HexaQ0Q0RT0_NmlFluxK_NmlzBas_Proj(FullMatrix &NmlFluxWG,
                                           const Hexa &hexa, Mat3 &MatK,
                                           Quadri3d SQF[6], int sign[6],
                                           const GaussQuad &GQH,
                                           const GaussQuad &GQQ);
*/

// WG3d: Hexa(Q0,Q0;RT[0]):
// 6 faces, 7 WG bas.fxns., 6 RT[0] nmlz.bas.fxns., No projection for Kw
int WG3d_HexaQ0Q0RT0_NmlFluxK_NmlzBas_NoProj(FullMatrix &NmlFluxWG,
                                             const Hexa &hexa, Mat3 &MatK,
                                             Quadri3d SQF[6], int sign[6],
                                             const GaussQuad &GQH,
                                             const GaussQuad &GQQ);


/*
////////////////////////////////////////////////////////////////////////////////
// For WG(Q0,Q1;P1^3)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa1(Q0,Q1;P1^3): Coeff.disc.wk.grad. in nmlz.bas.: matrix 12-by-25
// 12 nmlz.bas.fxns as rows;  25 WG bas.fxns. as columns
FullMatrix WG3d_Hexa1_CofDiscWkGrad_NmlzBas(const Hexa &hexa,
                                            Quadri3d SQF[6], int sign[6],
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ);

// WG3d: Hexa1(Q0,Q1;P1^3): Element grad-grad matrix w/ permeability MatK: 25*25
FullMatrix WG3d_Hexa1_EltGradGradMatK(const Hexa &hexa, Mat3 MatK,
                                      Quadri3d SQF[6], int sign[6],
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ);

// WG3d: Hexa1(Q0,Q1;P1^3): 6*25
FullMatrix WG3d_Hexa1_NmlFluxK_NmlzBas(const Hexa &hexa, Mat3 &MatK,
                                       Quadri3d SQF[6], int sign[6],
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);
*/


////////////////////////////////////////////////////////////////////////////////
// For WG(P1,P0;P0^3)Hexa with stabilizer
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P1,P0;P0^3): Coeff. of disc.wk.grad. in P0^3 nat.bas. (10-by-3)
FullMatrix WG3d_HexaP1P0P03_CofDiscWkGrad_NatBas(const Hexa &hexa,
                                                 Quadri3d SQF[6], int sign[6],
                                                 const GaussQuad &GQH,
                                                 const GaussQuad &GQQ);

// WG3d: Hexa(P1,P0;P0^3): Element grad-grad matrix (10-by-10)
FullMatrix WG3d_HexaP1P0P03_EltGradGradMatK(const Hexa &hexa, Mat3 &MatK,
                                            Quadri3d SQF[6], int sign[6],
                                            const GaussQuad &GQH,
                                            const GaussQuad &GQQ);

// WG3d: Hexa(P1,P0;P0^3): Elementwise stabilizer matrix (10-by-10)
FullMatrix WG3d_HexaP1P0P03_EltStabMat(const Hexa &hexa,
                                       Quadri3d SQF[6], int sign[6],
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);


////////////////////////////////////////////////////////////////////////////////
// For WG(P0^3,P0^3;AT0^3,P0)Hexa
////////////////////////////////////////////////////////////////////////////////
/*
// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Coeff. of disc.wk.grad. in AT0^3 nmlz.bas. for 21 WG bas.fxns
// 21-by-18 matrix
FullMatrix WG3d_HexaP03P03AT03P0_DiscWkGradBasFxn_NmlzPiolaBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                                  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Coeff. of disc.wk.div. in P0 nmlz.bas. for 21 WG bas.fxns
// 21-dim vector
Vector WG3d_HexaP03P03AT03P0_DiscWkDivBasFxn_NmlzBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                     const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Elementwise discrete weak strain-strain matrix
FullMatrix WG3d_HexaP03P03AT03P0_EltStrnStrnMat_NmlzPiolaBas(const Hexa &hexa,
                                                             Quadri3d SQF[6], int sign[6],
                                                             const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P0^3,P0^3;AT0^3,P0):
// Elementwise discrete weak div-div matrix
FullMatrix WG3d_HexaP03P03AT03P0_EltDivDivMat_NmlzBas(const Hexa &hexa, Quadri3d SQF[6], int sign[6],
                                                      const GaussQuad &GQH, const GaussQuad &GQQ);
*/

////////////////////////////////////////////////////////////////////////////////
// For WG(Q0^3,Q0^3;RT[0]^3,Q0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Coeff. of disc.wk.grad. in RT[0]^3 nmlz.bas. for 21 WG bas.fxns 
// 21-by-18 matrix 
FullMatrix WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkGradBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Coeff. of disc.wk.div. in Q0 nmlz.bas. for 21 WG bas.fxns
// 21-dim vector
Vector WG3d_HexaQ03Q03RT03Q0_CofNmlzBas_DiscWkDivBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Elementwise discrete weak strain-strain matrix
FullMatrix WG3d_HexaQ03Q03RT03Q0_NmlBas_EltStrnStrnMat(const Hexa &hexa,
                                                       Quadri3d SQF[6], int sign[6],
                                                       const GaussQuad &GQH,
                                                       const GaussQuad &GQQ);

// WG3d: Hexa(Q0^3,Q0^3;RT[0]^3,Q0):
// Elementwise discrete weak div-div matrix
FullMatrix WG3d_HexaQ03Q03RT03Q0_NmlBas_EltDivDivMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);


////////////////////////////////////////////////////////////////////////////////
// For WG(P1^3,Prm;P0^{3x3},P0)Hexa
////////////////////////////////////////////////////////////////////////////////

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Coeff. of disc.wk.grad. in P0^{3x3} nmlz.bas. for 48 WG bas.fxns
// 48-by-9 matrix
FullMatrix WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkGradBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Coeff. of disc.wk.div. in P0 nmlz.bas. for 48 WG bas.fxns
// 48-dim vector
Vector WG3d_HexaP13PrmP033P0_CofNmlzBas_DiscWkDivBasFxn(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise discrete weak strain-strain matrix
// 48-by-48 matrix
FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltStrnStrnMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise discrete weak div-div matrix
// 48-by-48 matrix
FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltDivDivMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

// WG3d: Hexa(P1^3,Prm;P0^{3x3},P0):
// Elementwise stabilization matrix
// 48-by-48 matrix
FullMatrix WG3d_HexaP13PrmP033P0_NmlzBas_EltStabMat(
  const Hexa &hexa, Quadri3d SQF[6], int sign[6],
  const GaussQuad &GQH, const GaussQuad &GQQ);

#endif  // WG3D_HEXA_H
// WG3d_Hexa.h
