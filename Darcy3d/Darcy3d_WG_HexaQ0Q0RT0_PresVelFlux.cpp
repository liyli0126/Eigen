//
//  Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux.cpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#include "Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux.hpp"

#include <cmath>

#include "LinSys.h"
#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "HexaMesh.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"


int Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(Vector &NumerPresEm,
                                       FullMatrix &NumerVelCofRT0,
                                       FullMatrix &NumerFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Vector &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ)
{
  int labelFace[6], sign[6];
  PtVec3d EmCntr, FcCntr, nml;
  Quadri3d quadri, SQF[6];
  Vector cof(7), fluxtmp(6), vctmp(6);
  FullMatrix CDWG(7,6), TransCDWG(6,7);
  FullMatrix KtransMat(12,6), NumerFluxRT0(6,6), ProjMat(6,12);

  // Setup
  int Num0E = mesh.numberElements();

  // Extracting elementwise (constant) pressure
  NumerPresEm = sln.getSubVector(1,mesh.numberElements());

  //
  NumerVelCofRT0.resize(Num0E,6);
  for (int le=0; le<Num0E; ++le) {
    int labele = le + mesh.beginLabelElement();
    mesh.getElementFace(labele, labelFace);

    Hexa hexa = mesh.element(labele);
    hexa.enrich();

    Mat3 MatK = PermK[le];
    
    cof(1) = sln(labele);
    for (int k=2; k<=7; ++k)  cof(k) = sln(Num0E+labelFace[k-2]);
    
    CDWG = WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(hexa, GQH, GQQ);
    TransCDWG = transpose(CDWG);
    
    // std::cout << "CDWG=" << CDWG << "\n" << std::flush;
    
    Hdiv_Hexa01_KtransMat(KtransMat, MatK);
    Hdiv_Hexa01_ProjP1v3RT0(ProjMat, hexa, GQH);

    // std::cout << "KtransMat=" << KtransMat << "\n" << std::flush;
    // std::cout << "ProjMat=" << ProjMat << "\n" << std::flush;

    vctmp = ((ProjMat*KtransMat) * TransCDWG) * cof;
    for (int j=1; j<=6; ++j)  NumerVelCofRT0(labele,j) = vctmp(j);

    // Computing the normal fluxes of the 6 RT[0] bas.fxns. on 6 faces
    //
    for (int i=0; i<6; ++i)  SQF[i] = mesh.face(labelFace[i]);
    //
    EmCntr = hexa.trilinearmapping(0.5, 0.5, 0.5);
    for (int i=0; i<6; ++i) {
      sign[i] = 1;
      quadri = SQF[i];
      FcCntr = quadri.bilinearmapping(0.5, 0.5);
      nml = quadri.normal(0.5, 0.5);
      if (dotProduct(nml,FcCntr-EmCntr)<0)  sign[i] = -1;
    }
    //
    Hdiv_HexaRT0_NmlFlux_NmlzBas(NumerFluxRT0, hexa, SQF, sign, GQQ);

    //
    fluxtmp = NumerFluxRT0 * vctmp;
    for (int i=1; i<=6; ++i)  NumerFlux(labele,i) = fluxtmp(i);
  }
    //std::cout << "NumerFlux" << NumerFlux << std::endl;

  return(0);  // If successful
}



/*
#include "LinSys.h"
#include "matrix.h"
#include "vector.h"


#include "cell3d.h"
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "HexaMesh.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"
#include <Eigen/Dense>

int Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(Eigen::VectorXd &NumerPresEm,
                                       Eigen::MatrixXd &NumerVelCofRT0,
                                       Eigen::MatrixXd &NumerFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Eigen::VectorXd &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ)
{
    int labelFace[6], sign[6];
    PtVec3d EmCntr, FcCntr, nml;
    Quadri3d quadri, SQF[6];
    Eigen::VectorXd cof(7), fluxtmp(6), vctmp(6);
    Eigen::MatrixXd CDWG(7, 6), TransCDWG(6, 7);
    Eigen::MatrixXd KtransMat(12, 6), NumerFluxRT0(6, 6), ProjMat(6, 12);

    // Setup
    int Num0E = mesh.numberElements();

    // Extracting elementwise (constant) pressure
    NumerPresEm = sln.segment(1, Num0E);

    NumerVelCofRT0.resize(Num0E, 6);
    for (int le = 0; le < Num0E; ++le) {
        int labele = le + mesh.beginLabelElement();
        mesh.getElementFace(labele, labelFace);

        Hexa hexa = mesh.element(labele);
        hexa.enrich();

        Mat3 MatK = PermK[le];
        
        cof(0) = sln(labele); // Indexing adjusted for Eigen
        for (int k = 1; k <= 6; ++k) {
            cof(k) = sln(Num0E + labelFace[k - 1]); // 0-based index
        }

        CDWG = WG3d_HexaQ0Q0RT0_CofDiscWkGrad_NmlzBas(hexa, GQH, GQQ);
        TransCDWG = CDWG.transpose();

        Hdiv_Hexa01_KtransMat(KtransMat, MatK);
        Hdiv_Hexa01_ProjP1v3RT0(ProjMat, hexa, GQH);

        // Calculate velocity coefficients
        vctmp = (ProjMat * KtransMat * TransCDWG * cof);
        for (int j = 0; j < 6; ++j) { // Eigen uses 0-based indexing
            NumerVelCofRT0(labele, j) = vctmp(j);
        }

        // Computing the normal fluxes of the 6 RT[0] basis functions on 6 faces
        for (int i = 0; i < 6; ++i) {
            SQF[i] = mesh.face(labelFace[i]);
        }
        EmCntr = hexa.trilinearmapping(0.5, 0.5, 0.5);
        for (int i = 0; i < 6; ++i) {
            sign[i] = 1;
            quadri = SQF[i];
            FcCntr = quadri.bilinearmapping(0.5, 0.5);
            nml = quadri.normal(0.5, 0.5);
            if (dotProduct(nml, FcCntr - EmCntr) < 0) {
                sign[i] = -1;
            }
        }

        // Calculate normal fluxes
        Hdiv_HexaRT0_NmlFlux_NmlzBas(NumerFluxRT0, hexa, SQF, sign, GQQ);

        // Compute flux
        fluxtmp = NumerFluxRT0 * vctmp;
        for (int i = 0; i < 6; ++i) {
            NumerFlux(labele, i) = fluxtmp(i);
        }
    }

    return 0;  // If successful
}
*/
