//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.cpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#include "Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.hpp"
#include <cmath>

#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"


int Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1(Vector &GlbVecDirichlet,
                                          Vector &GlbVecNeumann,
                                          double (*fxnpD)(PtVec3d),
                                          double (*fxnuN)(PtVec3d),
                                          const HexaMesh &mesh,
                                          const GaussQuad &GQQ)
{
  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  int DOFs = Num0E + Num1C;
  GlbVecDirichlet.resize(DOFs);
  GlbVecNeumann.resize(DOFs);

  // JL20150409: FOR NOW ASSUMING ALL DIRICHLET BOUNDARY CONDITIONS
  // Processing Dirichlet boundary conditions
  for (int labelc=1; labelc<=Num1C; ++labelc) {
    if (mesh.isBoundaryFace(labelc)>0) {
      Quadri3d quadri = mesh.face(labelc);
      double pDavg = 0, area = 0;
      for (int k=0; k<GQQ.numberQuadraturePoints(); ++k) {
        double xhat = GQQ.CartesianCoordinate(k,0);
        double yhat = GQQ.CartesianCoordinate(k,1);
        double jac = quadri.JacobianDeterminant(xhat, yhat);
        jac = fabs(jac);  // JL20150425: REALLY NEEDED
        PtVec3d qp(0,0,0);
        for (int j=0; j<GQQ.numberVertices(); ++j) {
          qp = qp + GQQ.baryCoordinate(k,j)*quadri.vertex(j);
        }
        double pval = fxnpD(qp);
        pDavg += GQQ.weight(k) * jac * pval;
        area += GQQ.weight(k) * jac;
      }
      pDavg /= area;
      GlbVecDirichlet(Num0E+labelc) = pDavg;
    }
  }

  // Processing Neumann boundary conditions
  // None

  return(0);  // If successful
}
