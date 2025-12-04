//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.cpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#include "Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.hpp"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"


int Darcy3d_WG_HexaQ0Q0RT0_AsmSource1(Vector &GlbVecSource,
                                      double (*fxnf)(PtVec3d),
                                      const HexaMesh &mesh,
                                      const GaussQuad &GQH)
{
  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  int DOFs = Num0E + Num1C;
  GlbVecSource.resize(DOFs);

  // Assembling GlbRHS
  for (int labele=1; labele<=Num0E; ++labele) {
    Hexa hexa = mesh.element(labele);
    hexa.enrich();
    double fintgrl = 0;
    // double vol = 0;
    for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
      double xhat = GQH.CartesianCoordinate(k,0);
      double yhat = GQH.CartesianCoordinate(k,1);
      double zhat = GQH.CartesianCoordinate(k,2);
      double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
      jac = fabs(jac);  // JL20150422: DON'T USE fabs?
      PtVec3d qp(0,0,0);
      for (int j=0; j<GQH.numberVertices(); ++j) {
        qp = qp + GQH.baryCoordinate(k,j)*hexa.vertex(j);
      }
      // qp = hexa.trilinearmapping(xhat,yhat,zhat);
      double fval = fxnf(qp);
      fintgrl += GQH.weight(k) * jac * fval;
    }
    GlbVecSource(labele) = fintgrl;
  }

  return(0);  // If successful
}
