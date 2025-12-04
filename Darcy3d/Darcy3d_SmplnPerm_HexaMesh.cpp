// Darcy3d_SmplnPerm_HexaMesh.cpp
// Sampling permeability on a hexa. mesh using a Gaussian quadrature
// James Liu, Graham Harper, ColoState; 2014/07--2017/02

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "mat3.h"
#include "PtVec3d.h"


int Darcy3d_SmplnPerm_HexaMesh(Mat3 *PermK, Mat3 (*fxnMatK)(PtVec3d),
                               const HexaMesh &mesh,
                               const GaussQuad &GQH)
{
  int k, labele, le;
  double jac, vol, xhat, yhat, zhat;
  PtVec3d qp;
  Hexa hexa;
  Mat3 MatK, MatKavg;
  
  // std::cout << "In SmplnPerm_HexaMesh..." << "\n";
  
  for (le=0; le<mesh.numberElements(); ++le) {
    labele = le + mesh.beginLabelElement();
    hexa = mesh.element(labele);
    hexa.enrich();  // NOTE: enrich() needs to be done before proceeding
    
    MatKavg.reset();
    vol = 0;
    for (k=0; k<GQH.numberQuadraturePoints(); ++k) {
      xhat = GQH.CartesianCoordinate(k,0);
      yhat = GQH.CartesianCoordinate(k,1);
      zhat = GQH.CartesianCoordinate(k,2);
      jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
      jac = fabs(jac);  // JL20150422: It is "safer" to use fabs
      /*
       qp = PtVec3d(0,0,0);
       for (j=0; j<GQH.numberVertices(); ++j) {
       qp = qp + GQH.baryCoordinate(k,j) * hexa.vertex(j);
       }
       */
      // Alternative:
      qp = hexa.trilinearmapping(xhat, yhat, zhat);
      MatK = fxnMatK(qp);
      MatKavg = MatKavg + (GQH.weight(k) * jac) * MatK;
      vol += GQH.weight(k)*jac;
    }
    PermK[le] = (1/vol) * MatKavg;
  }
  
  return(0);  // if successful
}

// Darcy3d_SmplnPerm_HexaMesh.cpp