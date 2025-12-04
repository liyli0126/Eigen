// JL20171102: TO BE FINISHED/REVISED By Graham
// Hdiv3d.h
// Broadly defined "Hdiv" bases 
// James Liu, Graham Harper, ColoState; 2014/07--2017/11

#include "matrix.h"
#include "cell3d.h"
#include "GaussQuad.h"
#include "mat3.h"
#include "PtVec3d.h"


////////////////////////////////////////////////////////////////////////////////
// For bricks
////////////////////////////////////////////////////////////////////////////////

int Hdiv_BrickRT0_GramMat_NmlzBas(double x0, double y0, double z0,
  double x1, double y1, double z1, double GM[6], double GMI[6]);

int Hdiv_BrickRT0_NmlFlux_NmlzBas(double x0, double y0, double z0,
  double x1, double y1, double z1, double NmlFlux[][6]);

FullMatrix Hdiv_BrickRT0_GramMatK_NmlzBas(const Brick &brick, Mat3 &MatK);

FullMatrix Hdiv_BrickRT0_NmlFlux_NmlzBas(const Brick &brick);

int Hdiv_BrickP1v3_GramMat_NmlzBas(DiagMatrix &GM, const Brick &brick);

int Hdiv_BrickP1v3_GramMatK_NmlzBas(FullMatrix &GMK, const Brick &brick,
                                    Mat3 &MatK);


////////////////////////////////////////////////////////////////////////////////
// For tetrahedra
////////////////////////////////////////////////////////////////////////////////

int Hdiv_TetraRT0_GramMat_NmlzBas(double x[], double y[], double z[],
                                  double GM[4], double GMI[4]);

// GM: Order 4 diag.mat. for the normalized RT0 basis
int Hdiv_TetraRT0_GramMat_NmlzBas(DiagMatrix &GM, const Tetra &tetra,
                                  const GaussQuad &GQTe);

// GMK: 4x4 "Gram" matrix with full permeability matrix MatK
int Hdiv_TetraRT0_GramMatK_NmlzBas(FullMatrix &GMK, const Tetra &tetra,
                                   Mat3 &MatK, const GaussQuad &GQTe);

// NmlFlux: 4x4 matrix, faces as rows, basis functions as columns
/*
int Hdiv_TetraRT0_NmlFlux_NmlzBas(FullMatrix &NmlFlux, const Tetra &tetra,
                                  Tri3d STF[4], int sign[4],
                                  const GaussQuad &GQT);
*/
int Hdiv_TetraRT0_NmlFlux_NmlzBas(FullMatrix &NmlFlux, const Tetra &tetra);

// KtransMat: 12x4 matrix, P1v3 as rows, RT0 as columns, not -K
int Hdiv_Tetra01_KtransMat(FullMatrix &KtransMat, Mat3 &MatK);

// ProjMat: 4X12 matrix, RT0 as rows, P1v3 as columns
int Hdiv_Tetra01_ProjP1v3RT0(FullMatrix &ProjMat, const Tetra &tetra,
                             const GaussQuad &GQH);

// GM: Order 15 full matrix for the normalized RT1 basis
int Hdiv_TetraRT1_GramMat_NmlzBas(FullMatrix &GM, const Tetra &tetra,
                                  const GaussQuad &GQTe);

// Matrix-version, to be used for elasticity (displacement gradient)
int Hdiv_TetraRT03_GramMat_NmlzBas(DiagMatrix &GM, const Tetra &tetra,
                                   const GaussQuad &GQTe);

// Matrix-version, to be used for elasticity (strain)
int Hdiv_TetraRT03_GramMatAvg_NmlzBas(FullMatrix &GMA, const Tetra &tetra,
                                      const GaussQuad &GQTe);


////////////////////////////////////////////////////////////////////////////////
// For hexahedra
////////////////////////////////////////////////////////////////////////////////

int Hdiv_HexaAT0_BasFxnVal_NmlzPiolaBas(PtVec3d *val, const Hexa &hexa,
                                        const double xhat, const double yhat,
                                        const double zhat);

int Hdiv_HexaAT0_GramMat_NmlzPiolaBas(FullMatrix &GM, const Hexa &hexa,
                                      const GaussQuad &GQH);

int Hdiv_HexaAT0_GramMatK_NmlzPiolaBas(FullMatrix &GM, FullMatrix &GMK,
                                       const Hexa &hexa, Mat3 MatK,
                                       const GaussQuad &GQH);

int Hdiv_HexaAT0_ProjCof_NmlzPiolaBas(FullMatrix &ProjCof,
                                      const Hexa &hexa, Mat3 MatK,
                                      const GaussQuad &GQH);

int Hdiv_HexaAT0_NmlFlux_NmlzPiolaBas(FullMatrix &NmlFluxAT0, const Hexa &hexa,
                                      Quadri3d SQF[6], int sign[6],
                                      const GaussQuad &GQQ);

int Hdiv_HexaRT0_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                                 const GaussQuad &GQH);

int Hdiv_HexaRT0_GramMatK_NmlzBas(FullMatrix &GM, FullMatrix &GMK,
                                  const Hexa &hexa, Mat3 MatK,
                                  const GaussQuad &GQH);

int Hdiv_HexaRT0_NmlFlux_NmlzBas(FullMatrix &NmFluxRT0, const Hexa &hexa,
                                 Quadri3d SQF[6], int sign[6],
                                 const GaussQuad &GQQ);

int Hdiv_Hexa01_KtransMat(FullMatrix &KtransMat, Mat3 &MatK);

int Hdiv_Hexa01_ProjP1v3RT0(FullMatrix &ProjMat, const Hexa &hexa,
                            const GaussQuad &GQH);

int Hdiv_HexaBDM1_GramMatK_NmlzBas(FullMatrix &GM, FullMatrix &GMK,
                                   const Hexa &hexa, Mat3 MatK,
                                   const GaussQuad &GQH);

int Hdiv_Hexa1_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                               const GaussQuad &GQH);

int Hdiv_Hexa1_GramMatK_NmlzBas(FullMatrix &GMK, const Hexa &hexa,
                                Mat3 &MatK, const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (displacement gradient)
int Hdiv_HexaAT03_GramMat_NmlzPiolaBas(FullMatrix &GM, const Hexa &hexa,
                                       const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (strain)
int Hdiv_HexaAT03_GramMatAvg_NmlzPiolaBas(FullMatrix &GMA, const Hexa &hexa,
                                          const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (displacement gradient)
int Hdiv_HexaRT03_GramMat_NmlzBas(FullMatrix &GM, const Hexa &hexa,
                                  const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (strain)
int Hdiv_HexaRT03_GramMatAvg_NmlzBas(FullMatrix &GMA, const Hexa &hexa,
                                     const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (displacement gradient)
int Hdiv_HexaP033_GramMat_NmlzBas(DiagMatrix &GM, const Hexa &hexa,
                                  const GaussQuad &GQH);

// Matrix-version, to be used for elasticity (strain)
int Hdiv_HexaP033_GramMatAvg_NmlzBas(FullMatrix &GMA, const Hexa &hexa,
                                     const GaussQuad &GQH);


////////////////////////////////////////////////////////////////////////////////
// Hdiv3d.h
