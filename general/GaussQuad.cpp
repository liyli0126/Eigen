// GaussQuad.cpp
// James Liu, Graham Harper, ColoState; 2007/01--2017/02

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "GaussQuad.h"


// GaussQuad: Destructor

GaussQuad::~GaussQuad()
{
  delete[] absc;
  delete[] wt;
  for (int i=0; i<dim; ++i)  delete[] crd[i];
  delete[] crd;
}


// GaussQuad: (Generalized) bary coordinate

double GaussQuad::baryCoordinate(int i, int j) const
{
  return absc[i*NumVrtx+j];
}


// GaussQuad: saving to a data file

void GaussQuad::save2file(char *filename) const
{
  int i, j;
  double sum;
  FILE *fp;
  
  if (NULL==(fp=fopen(filename, "w"))) {
    puts("Open data file failed.  Exit.");
    exit(-1);
  }
  
  fprintf(fp, "%3d  %2d\n", NumQuadPts, NumVrtx);

  for (i=0; i<NumQuadPts; ++i) {
    sum = 0;
    for (j=0; j<NumVrtx; ++j) {
      fprintf(fp, "%9.6f  ", absc[i*NumVrtx+j]);
      sum = sum + absc[i*NumVrtx+j];
    }
    fprintf(fp, "%6.3f  ", sum);
    fprintf(fp, "%9.6f\n", wt[i]);
  }

  fprintf(fp, "%3d  %2d\n", NumQuadPts, dim);

  for (i=0; i<NumQuadPts; ++i) {
    for (j=0; j<dim; ++j)
      fprintf(fp, "%9.6f  ", crd[i][j]);
    fprintf(fp, "\n");
  }

  fclose(fp);
  return;
}


// GaussQuad: setting up quadratures for intervals

void GaussQuad::setForInterval(int NumQuadPtsX)
{
  NumVrtx = 2;
  switch (NumQuadPtsX) {
    case 1:
      NumQuadPts = NumQuadPtsX;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2+0] = 0.5;  absc[0*2+1] = 0.5;
      // Setting up the weights
      wt[0] = 1;
      break;
    case 3:
      NumQuadPts = NumQuadPtsX;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2] = 0.88729833462074;  absc[0*2+1] = 0.11270166537926;
      absc[1*2] = 0.50000000000000;  absc[1*2+1] = 0.50000000000000;
      absc[2*2] = 0.11270166537926;  absc[2*2+1] = 0.88729833462074;
      // Setting up the weights
      wt[0] = 0.27777777777778;
      wt[1] = 0.44444444444444;
      wt[2] = 0.27777777777778;
      break;
    case 5:
      NumQuadPts = NumQuadPtsX;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2] = 0.95308992295000;  absc[0*2+1] = 0.04691007705000;
      absc[1*2] = 0.76923465505000;  absc[1*2+1] = 0.23076534495000;
      absc[2*2] = 0.50000000000000;  absc[2*2+1] = 0.50000000000000;
      absc[3*2] = 0.23076534495000;  absc[3*2+1] = 0.76923465505000;
      absc[4*2] = 0.04691007705000;  absc[4*2+1] = 0.95308992295000;
      // Setting up the weights
      wt[0] = 0.11846344250000;
      wt[1] = 0.23931433525000;
      wt[2] = 0.28444444445000;
      wt[3] = 0.23931433525000;
      wt[4] = 0.11846344250000;
      break;
    case 9:
      NumQuadPts = NumQuadPtsX;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2] = 0.98408011975400;  absc[0*2+1] = 0.01591988024600;
      absc[1*2] = 0.91801555366350;  absc[1*2+1] = 0.08198444633650;
      absc[2*2] = 0.80668571635050;  absc[2*2+1] = 0.19331428364950;
      absc[3*2] = 0.66212671170200;  absc[3*2+1] = 0.33787328829800;
      absc[4*2] = 0.50000000000000;  absc[4*2+1] = 0.50000000000000;
      absc[5*2] = 0.33787328829800;  absc[5*2+1] = 0.66212671170200;
      absc[6*2] = 0.19331428364950;  absc[6*2+1] = 0.80668571635050;
      absc[7*2] = 0.08198444633650;  absc[7*2+1] = 0.91801555366350;
      absc[8*2] = 0.01591988024600;  absc[8*2+1] = 0.98408011975400;
      // Setting up the weights
      wt[0] = 0.04063719418080;
      wt[1] = 0.09032408034750;
      wt[2] = 0.13030534820150;
      wt[3] = 0.15617353852000;
      wt[4] = 0.16511967750050;
      wt[5] = 0.15617353852000;
      wt[6] = 0.13030534820150;
      wt[7] = 0.09032408034750;
      wt[8] = 0.04063719418080;
      break;
    case 15:
      NumQuadPts = NumQuadPtsX;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2] = 0.006003740989757;  absc[0*2+1] = 0.993996259010243;
      absc[1*2] = 0.031363303799647;  absc[1*2+1] = 0.968636696200353;
      absc[2*2] = 0.075896708294786;  absc[2*2+1] = 0.924103291705214;
      absc[3*2] = 0.137791134319915;  absc[3*2+1] = 0.862208865680085;
      absc[4*2] = 0.214513913695731;  absc[4*2+1] = 0.785486086304269;
      absc[5*2] = 0.302924326461218;  absc[5*2+1] = 0.697075673538782;
      absc[6*2] = 0.399402953001283;  absc[6*2+1] = 0.600597046998717;
      absc[7*2] = 0.500000000000000;  absc[7*2+1] = 0.500000000000000;
      absc[8*2] = 0.600597046998717;  absc[8*2+1] = 0.399402953001283;
      absc[9*2] = 0.697075673538782;  absc[9*2+1] = 0.302924326461218;
      absc[10*2] = 0.785486086304269;  absc[10*2+1] = 0.214513913695731;
      absc[11*2] = 0.862208865680085;  absc[11*2+1] = 0.137791134319915;
      absc[12*2] = 0.924103291705214;  absc[12*2+1] = 0.075896708294786;
      absc[13*2] = 0.968636696200353;  absc[13*2+1] = 0.031363303799647;
      absc[14*2] = 0.993996259010243;  absc[14*2+1] = 0.006003740989757;
      // Setting up the weights
      wt[0] = 0.030753241996117;
      wt[1] = 0.070366047488108;
      wt[2] = 0.107159220467172;
      wt[3] = 0.139570677926154;
      wt[4] = 0.166269205816994;
      wt[5] = 0.186161000015562;
      wt[6] = 0.198431485327112;
      wt[7] = 0.202578241925561;
      wt[8] = 0.198431485327112;
      wt[9] = 0.186161000015562;
      wt[10] = 0.166269205816994;
      wt[11] = 0.139570677926154;
      wt[12] = 0.107159220467172;
      wt[13] = 0.070366047488108;
      wt[14] = 0.030753241996117;
      break;
    default:
      NumQuadPts = 1;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      absc[0*2+0] = 0.5;  absc[0*2+1] = 0.5;
      // Setting up the weights
      wt[0] = 1;
      break;
  }

  dim = 1;
  crd = new double*[NumQuadPts];
  for (int i=0; i<NumQuadPts; ++i) {
    crd[i] = new double[1];
    crd[i][0] = absc[i*2+1];
  }
  return;
}


// GaussQuad: setting up quadratures for rectangles

void GaussQuad::setForRectangle(int NumQuadPtsX, int NumQuadPtsY)
{
  int i, j, loc;
  GaussQuad GQX, GQY;
  GQX.setForInterval(NumQuadPtsX);
  GQY.setForInterval(NumQuadPtsY);
  
  NumVrtx = 4;
  NumQuadPts = GQX.numberQuadraturePoints() * GQY.numberQuadraturePoints();

  absc = new double[NumQuadPts*NumVrtx];
  wt = new double[NumQuadPts];

  dim = 2;
  crd = new double*[NumQuadPts];
  for (i=0; i<NumQuadPts; ++i)  crd[i] = new double[2];
  
  loc = 0;
  for (j=0; j<NumQuadPtsY; ++j) {
    for (i=0; i<NumQuadPtsX; ++i) {
      absc[4*loc+0] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,0);
      absc[4*loc+1] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,0);
      absc[4*loc+2] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,1);
      absc[4*loc+3] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,1);
      wt[loc] = GQX.weight(i) * GQY.weight(j);
      crd[loc][0] = GQX.baryCoordinate(i,1);
      crd[loc][1] = GQY.baryCoordinate(j,1);
      loc++;
    }
  }

  return;
}


// GaussQuad: setting up quadratures for triangles

void GaussQuad::setForTriangle(int NumQuadPtsT)
{
  NumVrtx = 3;

  switch (NumQuadPtsT) {
    case 1: {  // 1-point quadrature for triangles, Exact for poly. deg.<=1
      NumQuadPts = 1;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = 1.0/3;
      absc[0*3+0] = a;  absc[0*3+1] = a;  absc[0*3+2] = a;
      // Setting up the weights
      wt[0] = 1;
      break;
    }
    case 3: {  // 3-point quadrature for triangles, Exact for poly. deg.<=2
      NumQuadPts = 3;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a1 = 4.0/6;  double a2 = 1.0/6;
      absc[0*3+0] = a1;  absc[0*3+1] = a2;   absc[0*3+2] = a2;
      absc[1*3+0] = a2;  absc[1*3+1] = a1;   absc[1*3+2] = a2;
      absc[2*3+0] = a2;  absc[2*3+1] = a2;   absc[2*3+2] = a1;
      // Setting up the weights
      double c = 1.0/3;
      wt[0] = c;
      wt[1] = c;
      wt[2] = c;
      break;
    }
    case 7: {  // 7-point quadrature for triangles, Exact for poly. deg.<=5
      NumQuadPts = 7;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = 1/3.0;
      double b = (9+2*sqrt(15.0))/21;  // b = 0.797426985353087;
      double c =   (6-sqrt(15.0))/21;  // c = 0.101286507323456;
      double d =   (6+sqrt(15.0))/21;  // d = 0.470142064105115;
      double e = (9-2*sqrt(15.0))/21;  // e = 0.059715871789770;
      absc[0*3+0] = a;  absc[0*3+1] = a;  absc[0*3+2] = a;
      absc[1*3+0] = b;  absc[1*3+1] = c;  absc[1*3+2] = c;
      absc[2*3+0] = c;  absc[2*3+1] = b;  absc[2*3+2] = c;
      absc[3*3+0] = c;  absc[3*3+1] = c;  absc[3*3+2] = b;
      absc[4*3+0] = d;  absc[4*3+1] = d;  absc[4*3+2] = e;
      absc[5*3+0] = d;  absc[5*3+1] = e;  absc[5*3+2] = d;
      absc[6*3+0] = e;  absc[6*3+1] = d;  absc[6*3+2] = d;
      // Setting up the weights
      double u = 0.225;
      double v = (155.0-sqrt(15.0))/1200;
      double w = (155.0+sqrt(15.0))/1200;
      wt[0] = u;
      wt[1] = v;  wt[2] = v;  wt[3] = v;
      wt[4] = w;  wt[5] = w;  wt[6] = w;
      break;
    }
    case 13: {  // 13-point quadrature for triangles, Exact for poly. deg.<=7
      NumQuadPts = 13;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a0 = 1.0/3;
      double a = 0.479308067841923,  b = 0.260345966079039;
      double c = 0.869739794195568,  d = 0.065130102902216;
      double e = 0.638444188569809,  f = 0.312865496004875,  g = 0.048690315425316;
      absc[0*3+0] = a0;  absc[0*3+1] = a0;  absc[0*3+2] = a0;
      absc[1*3+0] = a;  absc[1*3+1] = b;  absc[1*3+2] = b;
      absc[2*3+0] = b;  absc[2*3+1] = a;  absc[2*3+2] = b;
      absc[3*3+0] = b;  absc[3*3+1] = b;  absc[3*3+2] = a;
      absc[4*3+0] = c;  absc[4*3+1] = d;  absc[4*3+2] = d;
      absc[5*3+0] = d;  absc[5*3+1] = c;  absc[5*3+2] = d;
      absc[6*3+0] = d;  absc[6*3+1] = d;  absc[6*3+2] = c;
      absc[7*3+0] = e;  absc[7*3+1] = f;  absc[7*3+2] = g;
      absc[8*3+0] = e;  absc[8*3+1] = g;  absc[8*3+2] = f;
      absc[9*3+0] = f;  absc[9*3+1] = e;  absc[9*3+2] = g;
      absc[10*3+0] = f;  absc[10*3+1] = g;  absc[10*3+2] = e;
      absc[11*3+0] = g;  absc[11*3+1] = e;  absc[11*3+2] = f;
      absc[12*3+0] = g;  absc[12*3+1] = f;  absc[12*3+2] = e;
      // Setting up the weights
      double s =-0.149570044467670,  t = 0.175615257433204;
      double u = 0.053347235608839,  v = 0.077113760890257;
      wt[0] = s;
      wt[1] = t;  wt[2] = t;  wt[3] = t;
      wt[4] = u;  wt[5] = u;  wt[6] = u;
      wt[7] = v;  wt[8] = v;  wt[9] = v;  wt[10] = v;  wt[11] = v;  wt[12] = v;
      break;
    }
    default: {  // 3-point quadrature for triangles, Exact for poly. deg.<=2
      NumQuadPts = 3;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a1 = 4.0/6;  double a2 = 1.0/6;
      absc[0*3+0] = a1;  absc[0*3+1] = a2;   absc[0*3+2] = a2;
      absc[1*3+0] = a2;  absc[1*3+1] = a1;   absc[1*3+2] = a2;
      absc[2*3+0] = a2;  absc[2*3+1] = a2;   absc[2*3+2] = a1;
      // Setting up the weights
      double c = 1.0/3;
      wt[0] = c;
      wt[1] = c;
      wt[2] = c;
      break;
    }
  }

  dim = 2;
  crd = new double*[NumQuadPts];
  for (int i=0; i<NumQuadPts; ++i) {
    crd[i] = new double[2];
    crd[i][0] = 0;
    crd[i][1] = 0;
  }

  return;
}


// GaussQuad: Setting up quadratures for bricks

void GaussQuad::setForBrick(int NumQuadPtsX, int NumQuadPtsY, int NumQuadPtsZ)
{
  int i, j, k, loc;
  
  GaussQuad GQX, GQY, GQZ;
  GQX.setForInterval(NumQuadPtsX);
  GQY.setForInterval(NumQuadPtsY);
  GQZ.setForInterval(NumQuadPtsZ);
  
  NumVrtx = 8;
  NumQuadPts = GQX.numberQuadraturePoints()
             * GQY.numberQuadraturePoints()
             * GQZ.numberQuadraturePoints();
  
  absc = new double[NumQuadPts*NumVrtx];
  wt = new double[NumQuadPts];

  dim = 3;
  crd = new double*[NumQuadPts];
  for (i=0; i<NumQuadPts; ++i)  crd[i] = new double[3];

  loc = 0;
  for (k=0; k<NumQuadPtsZ; ++k) {
    for (j=0; j<NumQuadPtsY; ++j) {
      for (i=0; i<NumQuadPtsX; ++i) {
        absc[8*loc+0] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,0) * GQZ.baryCoordinate(k,0);
        absc[8*loc+1] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,0) * GQZ.baryCoordinate(k,0);
        absc[8*loc+2] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,1) * GQZ.baryCoordinate(k,0);
        absc[8*loc+3] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,1) * GQZ.baryCoordinate(k,0);
        absc[8*loc+4] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,0) * GQZ.baryCoordinate(k,1);
        absc[8*loc+5] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,0) * GQZ.baryCoordinate(k,1);
        absc[8*loc+6] = GQX.baryCoordinate(i,1) * GQY.baryCoordinate(j,1) * GQZ.baryCoordinate(k,1);
        absc[8*loc+7] = GQX.baryCoordinate(i,0) * GQY.baryCoordinate(j,1) * GQZ.baryCoordinate(k,1);
        wt[loc] = GQX.weight(i) * GQY.weight(j) * GQZ.weight(k);
        crd[loc][0] = GQX.baryCoordinate(i,1);
        crd[loc][1] = GQY.baryCoordinate(j,1);
        crd[loc][2] = GQZ.baryCoordinate(k,1);
        loc++;
      }
    }
  }

  return;
}


// GaussQuad: setting up quadratures for tetrahedra

void GaussQuad::setForTetrahedron(int NumQuadPtsTe)
{
  NumVrtx = 4;

  switch (NumQuadPtsTe) {
    case 1: {  // 1-point quadrature for tetrahedra, Exact for poly. deg.<=1
      NumQuadPts = 1;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = 0.25;  // 0.25=1/4
      absc[0*4+0] = a;  absc[0*4+1] = a;  absc[0*4+2] = a;  absc[0*4+3] = a;
      // Setting up the weights
      wt[0] = 1;
      break;
    }
    case 4: { // 4-point quadrature for tetrahedra, Exact for poly. deg.<=2
      NumQuadPts = 4;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = (5+3*sqrt(5.0))/20;  // = 0. 5854101966249685
      double b = (5 - sqrt(5.0))/20;  // = 0. 1381966011250105
      absc[0*4+0] = a;  absc[0*4+1] = b;  absc[0*4+2] = b;  absc[0*4+3] = b;
      absc[1*4+0] = b;  absc[1*4+1] = a;  absc[1*4+2] = b;  absc[1*4+3] = b;
      absc[2*4+0] = b;  absc[2*4+1] = b;  absc[2*4+2] = a;  absc[2*4+3] = b;
      absc[3*4+0] = b;  absc[3*4+1] = b;  absc[3*4+2] = b;  absc[3*4+3] = a;
      // Setting up the weights
      double s = 1.0/4;
      wt[0] = s;  wt[1] = s;  wt[2] = s;  wt[3] = s;  wt[4] = s;
      break;
    }
    case 5: {  // 5-point quadrature for tetrahedra, Exact for poly. deg.<=3
      NumQuadPts = 5;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = 0.25,  b = 3.0/6,  c = 1.0/6;
      absc[0*4+0] = a;  absc[0*4+1] = a;  absc[0*4+2] = a;  absc[0*4+3] = a;
      absc[1*4+0] = b;  absc[1*4+1] = c;  absc[1*4+2] = c;  absc[1*4+3] = c;
      absc[2*4+0] = c;  absc[2*4+1] = b;  absc[2*4+2] = c;  absc[2*4+3] = c;
      absc[3*4+0] = c;  absc[3*4+1] = c;  absc[3*4+2] = b;  absc[3*4+3] = c;
      absc[4*4+0] = c;  absc[4*4+1] = c;  absc[4*4+2] = c;  absc[4*4+3] = b;
      // Setting up the weights
      double r = -4.0/5,  s = 9.0/20;
      wt[0] = r;
      wt[1] = s;  wt[2] = s;  wt[3] = s;  wt[4] = s;
      break;
    }
    case 11: {  // 11-point quadrature for tetrahedra, Exact for poly. deg.<=4
      NumQuadPts = 11;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double a = 0.25,  b = 11.0/14,  c = 1.0/14;
      double d = (1+sqrt(5.0/14))/4;  // = 0.3994035761667992
      double e = (1-sqrt(5.0/14))/4;  // = 0.1005964238332008
      absc[0*4+0] = a;  absc[0*4+1] = a;  absc[0*4+2] = a;  absc[0*4+3] = a;
      absc[1*4+0] = b;  absc[1*4+1] = c;  absc[1*4+2] = c;  absc[1*4+3] = c;
      absc[2*4+0] = c;  absc[2*4+1] = b;  absc[2*4+2] = c;  absc[2*4+3] = c;
      absc[3*4+0] = c;  absc[3*4+1] = c;  absc[3*4+2] = b;  absc[3*4+3] = c;
      absc[4*4+0] = c;  absc[4*4+1] = c;  absc[4*4+2] = c;  absc[4*4+3] = b;
      absc[5*4+0] = d;  absc[5*4+1] = d;  absc[5*4+2] = e;  absc[5*4+3] = e;
      absc[6*4+0] = d;  absc[6*4+1] = e;  absc[6*4+2] = d;  absc[6*4+3] = e;
      absc[7*4+0] = d;  absc[7*4+1] = e;  absc[7*4+2] = e;  absc[7*4+3] = d;
      absc[8*4+0] = e;  absc[8*4+1] = d;  absc[8*4+2] = d;  absc[8*4+3] = e;
      absc[9*4+0] = e;  absc[9*4+1] = d;  absc[9*4+2] = e;  absc[9*4+3] = d;
      absc[10*4+0] = e;  absc[10*4+1] = e;  absc[10*4+2] = d;  absc[10*4+3] = d;
      // Setting up the weights
      double r = -444.0/5625,  s = 2058.0/45000, t = 336.0/2250;
      wt[0] = r;
      wt[1] = s;  wt[2] = s;  wt[3] = s;  wt[4] = s;
      wt[5] = t;  wt[6] = t;  wt[7] = t;  wt[8] = t;  wt[9] = t;  wt[10] = t;
      break;
    }
    default: {
      NumQuadPts = 1;
      absc = new double[NumQuadPts*NumVrtx];
      wt = new double[NumQuadPts];
      // Setting up the abscissae
      double r = 0.25;  // 0.25=1/4
      absc[0*4+0] = r;  absc[0*4+1] = r;  absc[0*4+2] = r;  absc[0*4+3] = r;
      // Setting up the weights
      wt[0] = 1;
      break;
    }
  }

  dim = 3;
  crd = new double*[NumQuadPts];
  for (int i=0; i<NumQuadPts; ++i) {
    crd[i] = new double[3];
    crd[i][0] = 0;
    crd[i][1] = 0;
    crd[i][2] = 0;
  }

  return;
}

// GaussQuad.cpp