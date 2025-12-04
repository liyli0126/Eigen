// Darcy3d_WG_HexaMesh.h
// James Liu, Graham Harper, ColoState; 2014/07--2018/02

// Solving Darcy in 3d on a hexahedral mesh:
// -- Using WG(P0,P0;AT0) elements:     Good for convex hexahedron with flat faces
// -- Using WG(Q0,Q0;RT[0]) elements:   Good if not too much distortion
// -- Using WG(P1,P0;P0^3)+stabilizer:  Theoretically good
// -- Using WG(Q0,Q1;P1^3) elements:    Unsuccessful

#include <vector>

#include "matrix.h"
#include "vector.h"

#include "GaussQuad.h"
#include "HexaMesh.h"
#include "mat3.h"


////////////////////////////////////////////////////////////////////////////////
// For solving Darcy eqn. in 3d by the lowest order WG(Q0,Q0;RT[0])Hexa elements
// This is a new type of finite elements being investigated
////////////////////////////////////////////////////////////////////////////////

// The single-matrix approach: 4 functions with suffix "1"

// The Schur-complement approach: 4 functions

int Darcy3d_WG_HexaP0P0AT0_AsmGlbMats(BlockDiagMatrix &GlbMatEE,
                                      SparseBlockMatrix &GlbMatEC,
                                      SparseBlockMatrix &GlbMatCE,
                                      SparseBlockMatrix &GlbMatCC,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ);

int Darcy3d_WG_HexaP0P0AT0_AsmBndryConds(Vector &GlbVecD, Vector &GlbVecN,
                                         double (*fxnpD)(PtVec3d),
                                         double (*fxnuN)(PtVec3d),
                                         const HexaMesh &mesh,
                                         const GaussQuad &GQQ);

int Darcy3d_WG_HexaP0P0AT0_AsmSource(Vector &GlbRhsE, double (*fxnf)(PtVec3d),
                                     const HexaMesh &mesh, const GaussQuad &GQH);

int Darcy3d_WG_HexaP0P0AT0_ModiLinSys(SparseBlockMatrix &GlbMatEF,
                                      SparseBlockMatrix &GlbMatFE,
                                      SparseBlockMatrix &GlbMatFF,
                                      Vector &GlbVecE, Vector &GlbVecF,
                                      Vector &GlbVecD, Vector &GlbVecN,
                                      const HexaMesh &mesh);

// Common functions

int Darcy3d_WG_HexaP0P0AT0_Err(double &L2ErrPres, double &L2ErrVel,
                               double &L2ErrFlux, const HexaMesh &mesh,
                               const Vector &NumerPresEm,
                               const FullMatrix &NumerVelCofAT0,
                               double (*fxnp)(PtVec3d),
                               PtVec3d (*fxnu)(PtVec3d));

int Darcy3d_WG_HexaP0P0AT0_PresVelFlux(Vector &PresEm, FullMatrix &CofAT0Vel,
                                       FullMatrix &NmlFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Vector &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);

int Darcy3d_WG_HexaP0P0AT0_ProjPres(Vector &ProjPresEm, Vector &ProjPresFc,
                                    double (*fxnp)(PtVec3d),
                                    double (*fxnpD)(PtVec3d),
                                    const HexaMesh &mesh,
                                    const GaussQuad &GQH,
                                    const GaussQuad &GQQ);

int Darcy3d_WG_HexaP0P0AT0_VeriNmlFluxCnty(const HexaMesh &mesh,
                                           const FullMatrix &NmlFlux);


////////////////////////////////////////////////////////////////////////////////
// For solving Darcy eqn. in 3d by the lowest order WG(Q0,Q0;RT[0])Hexa elements
// This is a new type of finite elements being investigated
////////////////////////////////////////////////////////////////////////////////

// The single-matrix approach: 4 functions with suffix "1"

int Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1(Vector &GlbVecDirichlet,
                                          Vector &GlbVecNeumann,
                                          double (*fxnpD)(PtVec3d),
                                          double (*fxnuN)(PtVec3d),
                                          const HexaMesh &mesh,
                                          const GaussQuad &GQQ);

/*
int Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds(Vector &GlbVecD, Vector &GlbVecN,
                                         double (*fxnpD)(PtVec3d),
                                         double (*fxnuN)(PtVec3d),
                                         const Vector &ABC, const HexaMesh &mesh,
                                         const GaussQuad &GQQ);
*/

int Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1(SparseMatrix &GlbMat,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      const GaussQuad &GQH, const GaussQuad &GQQ);

int Darcy3d_WG_HexaQ0Q0RT0_AsmSource1(Vector &GlbVecSource,
                                      double (*fxnf)(PtVec3d),
                                      const HexaMesh &mesh,
                                      const GaussQuad &GQH);

int Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1(SparseMatrix &GlbMat, Vector &GlbRHS,
                                       Vector &GlbVecSource,
                                       Vector &GlbVecDirichlet,
                                       Vector &GlbVecNeumann,
                                       const HexaMesh &mesh);

// The Schur-complement approach: 4 functions

int Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMats(BlockDiagMatrix &GlbMatEE,
                                      SparseBlockMatrix &GlbMatEC,
                                      SparseBlockMatrix &GlbMatCE,
                                      SparseBlockMatrix &GlbMatCC,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ);

int Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds(Vector &GlbVecD, Vector &GlbVecN,
                                         double (*fxnpD)(PtVec3d),
                                         double (*fxnuN)(PtVec3d),
                                         const HexaMesh &mesh,
                                         const GaussQuad &GQQ);

int Darcy3d_WG_HexaQ0Q0RT0_AsmSource(Vector &GlbRhsE, double (*fxnf)(PtVec3d),
                                     const HexaMesh &mesh, const GaussQuad &GQH);

int Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys(SparseBlockMatrix &GlbMatEF,
                                      SparseBlockMatrix &GlbMatFE,
                                      SparseBlockMatrix &GlbMatFF,
                                      Vector &GlbVecE, Vector &GlbVecF,
                                      Vector &GlbVecD, Vector &GlbVecN,
                                      const HexaMesh &mesh);

// Common functions 

int Darcy3d_WG_HexaQ0Q0RT0_Err(double &L2ErrPres, double &L2ErrVel,
                               double &L2ErrFlux, const HexaMesh &mesh,
                               const Vector &NumerPresEm,
                               const FullMatrix &NumerVelCofRT0,
                               double (*fxnp)(PtVec3d),
                               PtVec3d (*fxnu)(PtVec3d));

int Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(Vector &PresEm, FullMatrix &CofRT0Vel,
                                       FullMatrix &NmlFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Vector &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);

int Darcy3d_WG_HexaQ0Q0RT0_ProjPres(Vector &ProjPresEm, Vector &ProjPresFc,
                                    double (*fxnp)(PtVec3d),
                                    double (*fxnpD)(PtVec3d),
                                    const HexaMesh &mesh,
                                    const GaussQuad &GQH,
                                    const GaussQuad &GQQ);

int Darcy3d_WG_HexaQ0Q0RT0_VeriNmlFluxCnty(const HexaMesh &mesh,
                                           const FullMatrix &NmlFlux);

////////////////////////////////////////////////////////////////////////////////
// For solving Darcy in 3d by WG(P1,P0;P0^3)Hexa elements with a stabilizer
////////////////////////////////////////////////////////////////////////////////

// The single-matrix approach: 4 functions with suffix "1"

int Darcy3d_WG_HexaP1P0P03_AsmBndryConds1(Vector &GlbVecDirichlet,
                                          Vector &GlbVecNeumann,
                                          double (*fxnpD)(PtVec3d),
                                          double (*fxnuN)(PtVec3d),
                                          const HexaMesh &mesh,
                                          const GaussQuad &GQQ);

int Darcy3d_WG_HexaP1P0P03_AsmGlbMat1(SparseMatrix &GlbMat,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      double rho,
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ);

int Darcy3d_WG_HexaP1P0P03_AsmSource1(Vector &GlbVecSource,
                                      double (*fxnf)(PtVec3d),
                                      const HexaMesh &mesh,
                                      const GaussQuad &GQH);

int Darcy3d_WG_HexaP1P0P03_ModiLinSys1(SparseMatrix &GlbMat, Vector &GlbRHS,
                                       const Vector &GlbVecSource,
                                       const Vector &GlbVecDirichlet,
                                       const Vector &GlbVecNeumann,
                                       const HexaMesh &mesh);

// The Schur-complement approach: 4 functions

int Darcy3d_WG_HexaP1P0P03_AsmBndryConds(Vector &GlbVecD, Vector &GlbVecN,
                                         double (*fxnpD)(PtVec3d),
                                         double (*fxnuN)(PtVec3d),
                                         const HexaMesh &mesh,
                                         const GaussQuad &GQQ);

int Darcy3d_WG_HexaP1P0P03_AsmGlbMats(BlockDiagMatrix &GlbMatEE,
                                      SparseBlockMatrix &GlbMatEF,
                                      SparseBlockMatrix &GlbMatFE,
                                      SparseBlockMatrix &GlbMatFF,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      double rho,
                                      const GaussQuad &GQH,
                                      const GaussQuad &GQQ);

int Darcy3d_WG_HexaP1P0P03_AsmSource(Vector &GlbVecE,
                                     double (*fxnf)(PtVec3d),
                                     const HexaMesh &mesh,
                                     const GaussQuad &GQH);

int Darcy3d_WG_HexaP1P0P03_ModiLinSys(SparseBlockMatrix &GlbMatEF,
                                      SparseBlockMatrix &GlbMatFE,
                                      SparseBlockMatrix &GlbMatFF,
                                      Vector &GlbVecE, Vector &GlbVecF,
                                      Vector &GlbVecD, Vector &GlbVecN,
                                      const HexaMesh &mesh);

// Common functions

int Darcy3d_WG_HexaP1P0P03_ProjPresVel(Vector &CofProjPresEm,
                                       Vector &CofProjPresFc,
                                       double (*fxnp)(PtVec3d),
                                       PtVec3d (*fxnu)(PtVec3d),
                                       const HexaMesh &mesh,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);

// Darcy3d_WG_HexaMesh.h
