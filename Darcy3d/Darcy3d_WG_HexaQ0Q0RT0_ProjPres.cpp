// Darcy3d_WG_HexaQ0Q0RT0_ProjPres.cpp
// James Liu, Graham Harper, ColoState; 2014/07--2017/02

#include <cmath>
#include "vector.h"
#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"

// It is good to use std::abs
// But it is preferred to use fabs for float, double arguments
// It ensures the result stays the same if I accidently peel std:: off the abs


// Projecting pressure

int Darcy3d_WG_HexaQ0Q0RT0_ProjPres(Vector &ProjPresEm, Vector &ProjPresFc,
                                    double (*fxnp)(PtVec3d),
                                    double (*fxnpD)(PtVec3d),
                                    const HexaMesh &mesh,
                                    const GaussQuad &GQH, const GaussQuad &GQQ)
{
  int j, k, labelc, labele, Num1E, Num2C;
  double area, jac, pavg, pDavg, pval, vol, xhat, yhat, zhat;
  Hexa hexa;
  Quadri3d quadri;
  PtVec3d cntr;
  PtVec3d qp;

  // Setup
  Num1E = mesh.numberElements();
  Num2C = mesh.numberFaces();
  ProjPresEm.resize(Num1E);
  ProjPresFc.resize(Num2C);

  // Computing elementwise pressure averages
  for (labele=1; labele<=mesh.numberElements(); ++labele) {
    hexa = mesh.element(labele);
    hexa.enrich();
    pavg = 0;  vol = 0;
    for (k=0; k<GQH.numberQuadraturePoints(); ++k) {
      xhat = GQH.CartesianCoordinate(k,0);
      yhat = GQH.CartesianCoordinate(k,1);
      zhat = GQH.CartesianCoordinate(k,2);
      jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
      jac = fabs(jac);  // JL20150422: DON'T USE abs?
      qp = PtVec3d(0,0,0);
      for (j=0; j<GQH.numberVertices(); ++j) {
        qp = qp + GQH.baryCoordinate(k,j) * hexa.vertex(j);
      }
      // qp = hexa.trilinearmapping(xhat,yhat,zhat);
      pval = fxnp(qp);
      pavg = pavg + GQH.weight(k) * jac * pval;
      vol = vol + GQH.weight(k) * jac;
    }
    /*
     cntr = hexa.center();
     x = cntr.xCrd();
     y = cntr.yCrd();
     z = cntr.zCrd();
     pavg = fxnp(x,y,z);
    */
    ProjPresEm(labele) = pavg/vol;
  }

  // Computing facewise pressure averages
  for (labelc=1; labelc<=mesh.numberFaces(); ++labelc) {
    quadri = mesh.face(labelc);
    pDavg = 0;  area = 0;
    for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
      xhat = GQQ.CartesianCoordinate(k,0);
      yhat = GQQ.CartesianCoordinate(k,1);
      jac = quadri.JacobianDeterminant(xhat,yhat);
      jac = fabs(jac);  // JL20150422: DON'T USE abs?
      /*
      qp = PtVec3d(0,0,0);
      for (j=0; j<GQQ.numberVertices(); ++j) {
        qp = qp + GQQ.baryCoordinate(k,j) * quadri.vertex(j);
      }
      */
      qp = quadri.bilinearmapping(xhat,yhat);
      pval = fxnpD(qp);
      pDavg = pDavg + GQQ.weight(k) * jac * pval;
      area = area + GQQ.weight(k) * jac;
    }
    /*
    cntr = quadri.center();
    x = cntr.xCrd();
    y = cntr.yCrd();
    z = cntr.zCrd();
    pDavg = fxnpD(x,y,z);
    */
    ProjPresFc(labelc) = pDavg/area;
  }

  return 0;  // If successful
}
