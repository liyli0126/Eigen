// Darcy3d_WG_HexaQ0Q0RT0_Err.cpp
// JL20150629: Using quadratures independent of other subroutines
// James Liu, Graham Harper, ColoState; 2014/07--2017/02

#include <cmath>

#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"


int Darcy3d_WG_HexaQ0Q0RT0_Err(double &L2ErrPres, double &L2ErrVel,
                               double &L2ErrFlux, const HexaMesh &mesh,
                               const Vector &NumerPresEm,
                               const FullMatrix &NumerVelCofRT0,
                               double (*fxnp)(PtVec3d),
                               PtVec3d (*fxnu)(PtVec3d))
                               // const GaussQuad &GQH, const GaussQuad &GQQ)
{
  int i, j, k, labelc, labele, lc, le;
  int labelFace[6], sign[6];
  double ar, ferr, flux, jac, perr, pval, sa, sum, sumf, sump, sumv, vol, vtmp;
  double x, y, z, xc, yc, zc, xhat, yhat, zhat, X, Y, Z;
  PtVec3d cntr, nml, NumerVel, qp, vel, verr;
  PtVec3d w[7];  // w[0] not used
  Quadri3d quadri;
  Hexa hexa;

  // Using its own quadratures
  GaussQuad GQH, GQQ;
  GQH.setForBrick(9,9,9);
  GQQ.setForRectangle(9,9);

  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  w[1] = PtVec3d(1,0,0);
  w[2] = PtVec3d(0,1,0);
  w[3] = PtVec3d(0,0,1);

  // Computing L2-norm of errors in pressure, velocity
  Vector ErrPresEm(Num0E);
  Vector ErrVelEm(Num0E);
  for (le=0; le<Num0E; ++le) {
    labele = le + mesh.beginLabelElement();
    hexa = mesh.element(labele);
    hexa.enrich();
    cntr = hexa.center();
    xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
    sump = 0;
    sumv = 0;
    for (i=0; i<GQH.numberQuadraturePoints(); ++i) {
      xhat = GQH.CartesianCoordinate(i,0);
      yhat = GQH.CartesianCoordinate(i,1);
      zhat = GQH.CartesianCoordinate(i,2);
      jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
      // jac = fabs(jac);  // It is safer to use fabs
      qp = hexa.trilinearmapping(xhat, yhat, zhat);
      // Alternative:
      // qp = PtVec3d(0,0,0);
      // for (j=0; j<GQH.numberVertices(); ++j) {
      // qp = qp + GQH.baryCoordinate(i,j) * hexa.vertex(j);
      // }
      pval = fxnp(qp);
      perr = pval - NumerPresEm(labele);
      sump += GQH.weight(i) * jac * (perr*perr);
      //
      x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
      X = x - xc;  Y = y - yc;  Z = z - zc;
      w[4] = PtVec3d(X,0,0);
      w[5] = PtVec3d(0,Y,0);
      w[6]= PtVec3d(0,0,Z);
      NumerVel = PtVec3d(0,0,0);
      for (j=1; j<=6; ++j) {
        NumerVel = NumerVel + NumerVelCofRT0(labele,j) * w[j];
      }
      vel = fxnu(qp);
      verr = vel - NumerVel;
      vtmp = verr.l2norm();
      sumv += GQH.weight(i) * jac * (vtmp*vtmp);
    }
    ErrPresEm(labele) = sqrt(fabs(sump));
    ErrVelEm(labele) = sqrt(fabs(sumv));
  }
  L2ErrPres = ErrPresEm.l2norm();
  L2ErrVel = ErrVelEm.l2norm();

  // Computing element volumes 
  Vector volElt(Num0E);
  for (le=0; le<Num0E; ++le) {
    labele = le + mesh.beginLabelElement();
    hexa = mesh.element(labele);
    hexa.enrich();
    vol = 0;
    for (i=0; i<GQH.numberQuadraturePoints(); ++i) {
      xhat = GQH.CartesianCoordinate(i,0);
      yhat = GQH.CartesianCoordinate(i,1);
      zhat = GQH.CartesianCoordinate(i,2);
      jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
      // qp = hexa.trilinearmapping(xhat, yhat, zhat);
      vol += GQH.weight(i) * jac;
    }
    volElt(labele) = fabs(vol);
  }

  // Computing face areas 
  Vector areaFace(Num1C);
  for (lc=0; lc<Num1C; ++lc) {
	  labelc = lc + mesh.beginLabelFace();
	  quadri = mesh.face(labelc);
	  double ar = 0;
    for (i=0; i<GQQ.numberQuadraturePoints(); ++i) {
      xhat = GQQ.CartesianCoordinate(i,0);
      yhat = GQQ.CartesianCoordinate(i,1);
	    jac = quadri.JacobianDeterminant(xhat, yhat);
	    ar += GQQ.weight(i) * jac;
	  }
	  areaFace(labelc) = fabs(ar);
  }

  // Computing L2-norm of error in flux
  // Just use analytical and numerical velocity normal components
  sum = 0;
  for (le=0; le<Num0E; ++le) {
    labele = le + mesh.beginLabelElement();
    mesh.getElementFace(labele, labelFace);
    
    hexa = mesh.element(labele);
    hexa.enrich();
    
    cntr = hexa.trilinearmapping(0.5, 0.5, 0.5);
    xc = cntr.xCrd();  yc = cntr.yCrd();  zc = cntr.zCrd();
    
    for (i=0; i<6; ++i) {
      quadri = mesh.face(labelFace[i]);
      sumf = 0;
      for (k=0; k<GQQ.numberQuadraturePoints(); ++k) {
        xhat = GQQ.CartesianCoordinate(k,0);
        yhat = GQQ.CartesianCoordinate(k,1);
        jac = quadri.JacobianDeterminant(xhat, yhat);
        // jac = fabs(jac);  // JL20150422: DON'T USE fabs?
        qp = quadri.bilinearmapping(xhat, yhat);
        nml = quadri.normal(xhat, yhat);
        x = qp.xCrd();  y = qp.yCrd();  z = qp.zCrd();
        X = x - xc;  Y = y - yc;  Z = z - zc;
        w[4] = PtVec3d(X,0,0);
        w[5] = PtVec3d(0,Y,0);
        w[6] = PtVec3d(0,0,Z);
        NumerVel = PtVec3d(0,0,0);
        for (j=1; j<=6; ++j) {
          NumerVel = NumerVel + NumerVelCofRT0(labele,j) * w[j];
        }
        vel = fxnu(qp);
        verr = vel - NumerVel;
		    ferr = dotProduct(verr,nml);
        sumf += GQH.weight(i) * jac * (ferr*ferr);
      }
      sum += sumf * volElt(labele)/areaFace(labelFace[i]);
    }
  }
  L2ErrFlux = sqrt(fabs(sum));

  // Computing the discrete L2-norm of the pressure error
  /*
  Vector discL2ErrPresElt(Num0E);
  for (le=0; le<Num0E; ++le) {
    labele = le + mesh.beginLabelElement();
    hexa = mesh.element(labele);
    hexa.enrich();
    cntr = hexa.center();
    pval = fxnp(cntr);
    perr = pval - NumerPresEm(labele);
    discL2ErrPresElt(labele) = sqrt(volElt(labele)) * fabs(perr);
  }
  discL2ErrPres = discL2ErrPresElt.l2norm();
  */
  
  return 0;  // If successful
}

// Darcy3d_WG_HexaQ0Q0RT0_Err.cpp
