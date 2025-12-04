// HexaMesh.cpp
// Class for hexahedral meshes
// James Liu, Graham Harper, ColoState; 2007/01--2017/03

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <stdio.h>
#include <stdlib.h>

#include "cell3d.h"
#include "HexaMesh.h"
#include "PtVec3d.h"
#include "EqnBC_Poisson3d_Ex12coscoscos.h"

// Implementation of the friend function for the THex algorithm
// #include "THex.cpp"

// Implementation of the friend function for hexa.mesh. regular refinement 
// #include "HexaMeshRegRefi.cpp"

#ifndef PI
#define PI 3.141592653589793
#endif


// HexaMesh: Constructor #0: Default empty constructor
HexaMesh::HexaMesh()
{
  NumNds = 0;  NumFcs = 0;  NumEms = 0;  NumBndryPcs = 0;
  BgnLblNd = 1;  BgnLblFc = 1;  BgnLblEm = 1;  BgnLblBndry = 1;
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  nd = 0;  OutNmlBndry = 0;
  flag = 0;
  NX = 0;  NY = 0;  NZ = 0;
}


// Constructor #1: Need to call enrich() later
// For a brick domain: Prescribed perturbation with magnitude delta

HexaMesh::HexaMesh(double xa, double xb,
                   double yc, double yd,
                   double ze, double zf,
                   int nx, int ny, int nz, double delta)
{
  int i, j, k, ld;
  double hx, hy, hz;
  double *x, *y, *z;

  // Counts
  NumNds = (nx+1)*(ny+1)*(nz+1);
  NumFcs = nx*ny*(nz+1) + nx*nz*(ny+1) + ny*nz*(nx+1);
  NumEms = nx*ny*nz;
  NumBndryPcs = 6;

  // Default: These labels start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;

  x = new double[nx+1];
  hx = (xb-xa)/nx;
  x[0] = xa;
  for (i=1; i<=nx-1; ++i)  x[i] = x[i-1] + hx;
  x[nx] = xb;  // Avoiding rounding errors

  y = new double[ny+1];
  hy = (yd-yc)/ny;
  y[0] = yc;
  for (j=1; j<=ny-1; ++j)  y[j] = y[j-1] + hy;
  y[ny] = yd;  // Avoiding rounding errors

  z = new double[nz+1];
  hz = (zf-ze)/nz;
  z[0] = ze;
  for (k=1; k<=nz-1; ++k)  z[k] = z[k-1] + hz;
  z[nz] = zf;  // Avoiding rounding errors
  
  // Info of node coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        nd[ld] = PtVec3d(x[i],y[j],z[k]);
        ld++;
      }
    }
  }

  // Peturbating the interior nodes with magnitude (delta)
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        // if ((i==0 || i==nx) && (j=0 || j==ny) && (k==0 || k==nz))  continue;
        double xx = nd[ld].xCrd();
        double yy = nd[ld].yCrd();
        double zz = nd[ld].zCrd();
        xx = xx + 0.03*sin(xx*PI*3)*cos(yy*PI*3)*cos(zz*PI*3) * delta;
        yy = yy - 0.04*cos(xx*PI*3)*sin(yy*PI*3)*cos(zz*PI*3) * delta;
        zz = zz + 0.05*cos(xx*PI*3)*cos(yy*PI*3)*sin(zz*PI*3) * delta;
        nd[ld] = PtVec3d(xx,yy,zz);
        ld++;
      }
    }
  }

  // Freeing...
  delete[] x, y, z;

  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;

  flag = 1;  // Need to call enrich() later
  NX = nx;  NY = ny;  NZ = nz;  // Logically rectangular
}


// Constructor #2: Need to call enrich() later
// For a brick domain: Random perturbations in x-,y-,z-directions
HexaMesh::HexaMesh(double xa, double xb,
                   double yc, double yd,
                   double ze, double zf,
                   int nx, int ny, int nz,
                   double xdelta, double ydelta, double zdelta)
{
  int i, j, k, ld;
  double hx, hy, hz;
  double *x, *y, *z;
  
  // Seed the random number generator
  // srand((int)time(0));

  // Counts
  NumNds = (nx+1)*(ny+1)*(nz+1);
  NumFcs = nx*ny*(nz+1) + nx*nz*(ny+1) + ny*nz*(nx+1);
  NumEms = nx*ny*nz;
  NumBndryPcs = 6;
  
  // Default: These labels start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;
  
  x = new double[nx+1];
  hx = (xb-xa)/nx;
  x[0] = xa;
  for (i=1; i<=nx-1; ++i)  x[i] = x[i-1] + hx;
  x[nx] = xb;  // Avoiding rounding errors
  
  y = new double[ny+1];
  hy = (yd-yc)/ny;
  y[0] = yc;
  for (j=1; j<=ny-1; ++j)  y[j] = y[j-1] + hy;
  y[ny] = yd;  // Avoiding rounding errors
  
  z = new double[nz+1];
  hz = (zf-ze)/nz;
  z[0] = ze;
  for (k=1; k<=nz-1; ++k)  z[k] = z[k-1] + hz;
  z[nz] = zf;  // Avoiding rounding errors
  
  // Info of node coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        nd[ld] = PtVec3d(x[i],y[j],z[k]);
        ld++;
      }
    }
  }
  
  // Peturbating the interior nodes with parameters xdelta, ydelta, zdelta
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        double xx = nd[ld].xCrd();  // Note that if any of xdelta, ydelta, zdelta is equal to 1,
        double yy = nd[ld].yCrd();  // there is a chance that degenerate elements can be created,
        double zz = nd[ld].zCrd();  // but using 1/2 works well
        
        // Don't change x,y,z directions on x,y,z boundaries, respectively
        if (i!=0 && i!=nx) {
          xx = xx + xdelta * hx * (1. * rand() / RAND_MAX - 1./2);
        }
        if (j!=0 && j!=ny) {
          yy = yy + ydelta * hy * (1. * rand() / RAND_MAX - 1./2);
        }
        if (k!=0 && k!=nz) {
          zz = zz + zdelta * hz * (1. * rand() / RAND_MAX - 1./2);
        }
        nd[ld] = PtVec3d(xx,yy,zz);
        ld++;
      }
    }
  }

  // Freeing...
  delete[] x, y, z;

  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;

  flag = 1;  // Need to call enrich() later
  NX = nx;  NY = ny;  NZ = nz;  // Logically rectangular
}


// Constructor #3: Need to call enrich() later
// For a brick domain: Trapezoidal perturbations in y-direction
HexaMesh::HexaMesh(double xa, double xb, int nx,
                   double yc, double yd, int ny, double alpha,
                   double ze, double zf, int nz)
{
  int i, j, k, ld;
  double hx, hy, hz;
  double xx, yy, zz;
  double *x, *y, *z;
  
  // Counts
  NumNds = (nx+1)*(ny+1)*(nz+1);
  NumFcs = nx*ny*(nz+1) + nx*nz*(ny+1) + ny*nz*(nx+1);
  NumEms = nx*ny*nz;
  NumBndryPcs = 6;
  
  // Default: These labels start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;
  
  x = new double[nx+1];
  hx = (xb-xa)/nx;
  x[0] = xa;
  for (i=1; i<=nx-1; ++i)  x[i] = x[i-1] + hx;
  x[nx] = xb;  // Avoiding rounding errors
  
  y = new double[ny+1];
  hy = (yd-yc)/ny;
  y[0] = yc;
  for (j=1; j<=ny-1; ++j)  y[j] = y[j-1] + hy;
  y[ny] = yd;  // Avoiding rounding errors
  
  z = new double[nz+1];
  hz = (zf-ze)/nz;
  z[0] = ze;
  for (k=1; k<=nz-1; ++k)  z[k] = z[k-1] + hz;
  z[nz] = zf;  // Avoiding rounding errors
  
  // Info of node coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        nd[ld] = PtVec3d(x[i],y[j],z[k]);
        ld++;
      }
    }
  }
  
  // Perturbation in y-direction only
  for (int k=0; k<=nz; ++k) {
    for (int j=1; j<ny; j=j+2) {
      for (int i=0; i<=nx; ++i) {
        ld = k*(ny+1)*(nx+1) + j*(nx+1) + i;
        xx = nd[ld].xCrd();
        yy = nd[ld].yCrd();
        zz = nd[ld].zCrd();
        if (i%2==1)  yy = yy + hy*alpha;
        if (i%2==0)  yy = yy - hy*alpha;
        nd[ld] = PtVec3d(xx,yy,zz);
      }
    }
  }

  // Freeing...
  delete[] x, y, z;
  
  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;
  
  flag = 1;  // Need to call enrich() later
  NX = nx;  NY = ny;  NZ = nz;  // Logically rectangular
}


// HexaMesh: Constructor #4: Need to call enrich() later
// For a brick domain: Trapezoidal perturbations in x-,y-directions

HexaMesh::HexaMesh(double xa, double xb,
                   double yc, double yd,
                   double ze, double zf,
                   int nx, int ny, int nz,
                   double xdelta, double ydelta)
{
  int i, j, k, ld;
  double hx, hy, hz;
  double xx, yy, zz;
  double xtilde, ytilde, ztilde;
  double *x, *y, *z;
  
  // Counts
  NumNds = (nx+1)*(ny+1)*(nz+1);
  NumFcs = nx*ny*(nz+1) + nx*nz*(ny+1) + ny*nz*(nx+1);
  NumEms = nx*ny*nz;
  NumBndryPcs = 6;

  // Default: These labels start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;
  
  x = new double[nx+1];
  hx = (xb-xa)/nx;
  x[0] = xa;
  for (i=1; i<=nx-1; ++i)  x[i] = x[i-1] + hx;
  x[nx] = xb;  // Avoiding rounding errors
  
  y = new double[ny+1];
  hy = (yd-yc)/ny;
  y[0] = yc;
  for (j=1; j<=ny-1; ++j)  y[j] = y[j-1] + hy;
  y[ny] = yd;  // Avoiding rounding errors
  
  z = new double[nz+1];
  hz = (zf-ze)/nz;
  z[0] = ze;
  for (k=1; k<=nz-1; ++k)  z[k] = z[k-1] + hz;
  z[nz] = zf;  // Avoiding rounding errors
  
  // Info of node coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        nd[ld] = PtVec3d(x[i],y[j],z[k]);
        ld++;
      }
    }
  }
  
  // Peturbating the interior nodes with parameters xdelta, ydelta
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ny; ++j) {
      for (i=0; i<=nx; ++i) {
        ld = k*(ny+1)*(nx+1) + j*(nx+1) + i;
        xx = nd[ld].xCrd();  // Setting any of xdelta or ydelta equal to 1
        yy = nd[ld].yCrd();  // results in degenerate trapezoids (missing top side)
        zz = nd[ld].zCrd();  // but using 1/2 works well
        if (i!=0 && i!=nx) {
          xx = xx + hx / 2 * xdelta * (((i+k) % 2) * 2 - 1);
        }
        if (j!=0 && j!=ny) {
          yy = yy + hy / 2 * ydelta * (((j+k) % 2) * 2 - 1);
        }
        nd[ld] = PtVec3d(xx,yy,zz);
      }
    }
  }
  
  // Freeing...
  delete[] x, y, z;
  
  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;

  flag = 1;  // Need to call enrich() later
  NX = nx;  NY = ny;  NZ = nz;  // Logically rectangular
}


// HexaMesh: Constructor #5: Need to call enrich() later
// For a cylindrical sector (nonclose)

HexaMesh::HexaMesh(double ri, double ro,
                   double alpha, double beta,
                   double zb, double zt,
                   int nr, int ntheta, int nz)
{
  int i, j, k, ld;
  double hr, htheta, hz;
  double *r, *theta, *z;

  // Counts
  NumNds = (nr+1)*(ntheta+1)*(nz+1);
  NumFcs = nr*ntheta*(nz+1) + nr*(ntheta+1)*nz + (nr+1)*ntheta*nz;
  NumEms = nr*ntheta*nz;
  NumBndryPcs = 6;

  // Default: These labels start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;
  
  // Partition in the radial(r)-direction
  r = new double[nr+1];
  hr = (ro-ri)/nr;
  r[0] = ri;
  for (i=1; i<=nr-1; ++i)  r[i] = r[i-1] + hr;
  r[nr] = ro;  // Avoiding rounding errors
  
  // Partition in the annular(theta)-direction
  theta = new double[ntheta+1];
  htheta = (beta-alpha)/ntheta;
  theta[0] = alpha;
  for (j=1; j<=ntheta-1; ++j)  theta[j] = theta[j-1] + htheta;
  theta[ntheta] = beta;  // Avoiding rounding errors
  
  // Partition in the z-direction
  z = new double [nz+1];
  hz = (zt-zb)/nz;
  z[0] = zb;  // z-direction bottom
  for (k=1; k<=nz-1; ++k)  z[k] = z[k-1] + hz;
  z[nz] = zt;  // z-direction top
  
  // Geometrical: node (x,y,z) coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nz; ++k) {
    for (j=0; j<=ntheta; ++j) {
      for (i=0; i<=nr; ++i) {
        nd[ld] = PtVec3d(r[i]*cos(theta[j]),r[i]*sin(theta[j]),z[k]);
        ld++;
      }
    }
  }

  // Freeing...
  delete[] r, theta, z;

  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;

  flag = 1;  // Need to call enrich() later
  NX = nr;  NY = ntheta;  NZ = nz;  // Logically rectangular
}


// HexaMesh: Constructor #6: Need to call enrich() later
// For a spherical sector (nonclose)

HexaMesh::HexaMesh(double ri, double ro, int nr,
                   double alpha, double beta, int ntheta,
                   double gamma, double delta, int nphi)
{
  int i, j, k, ld;
  double hr, htheta, hphi;
  double *r, *theta, *phi;
  
  // Counts
  NumNds = (nr+1)*(ntheta+1)*(nphi+1);
  NumFcs = nr*ntheta*(nphi+1) + nr*(ntheta+1)*nphi + (nr+1)*ntheta*nphi;
  NumEms = nr*ntheta*nphi;
  NumBndryPcs = 6;

  // Default: These lables start at 1
  BgnLblNd = 1;
  BgnLblFc = 1;
  BgnLblEm = 1;
  BgnLblBndry = 1;

  // Partition in the radial(r)-direction
  r = new double[nr+1];
  hr = (ro-ri)/nr;
  r[0] = ri;
  for (i=1; i<=nr-1; ++i)  r[i] = r[i-1] + hr;
  r[nr] = ro;  // Avoiding rounding errors
  
  // Partition in the angular(theta)-direction
  theta = new double[ntheta+1];
  htheta = (beta-alpha)/ntheta;
  theta[0] = alpha;
  for (j=1; j<=ntheta-1; ++j)  theta[j] = theta[j-1] + htheta;
  theta[ntheta] = beta;  // Avoiding rounding errors
  
  // Partition in the phi-direction
  phi = new double [nphi+1];
  hphi = (delta-gamma)/nphi;
  phi[0] = gamma;  // phi-direction bottom
  for (k=1; k<=nphi-1; ++k)  phi[k] = phi[k-1] + hphi;
  phi[nphi] = delta;  // phi-direction top
  
  // Geometrical: node (x,y,z) coordinates
  nd = new PtVec3d[NumNds];
  ld = 0;
  for (k=0; k<=nphi; ++k) {
    for (j=0; j<=ntheta; ++j) {
      for (i=0; i<=nr; ++i) {
        nd[ld] = PtVec3d(r[i]*sin(phi[k])*cos(theta[j]),
                         r[i]*sin(phi[k])*sin(theta[j]),
                         r[i]*cos(phi[k]));
        ld++;
      }
    }
  }

  // Freeing...
  //delete[] r, theta, phi;

  // Further initialization
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  OutNmlBndry = 0;

  flag = 1;  // Need to call enrich() later
  NX = nr;  NY = ntheta;  NZ = nphi;  // Logically rectangular
}


// HexaMesh: enriching mesh info for logically rectangular mesh (flag=1)
// Generating topological info

void HexaMesh::enrich()
{
  int i, j, k, l, labelElementA, labelElementB, labelVertex[8];
  int lc, le, lea, leb;

  if (flag==2)  exit;
    
  // Generating info of element-vs-(8)nodes
  LblEmNd = new int[NumEms][8];
  le = 0;
  for (k=1; k<=NZ; ++k) {
    for (j=1; j<=NY; ++j) {
      for (i=1; i<=NX; ++i) {
        labelVertex[0] = (k-1)*(NY+1)*(NX+1) + (j-1)*(NX+1) + i;
        labelVertex[1] = labelVertex[0] + 1;
        labelVertex[2] = labelVertex[0] + NX+2;
        labelVertex[3] = labelVertex[0] + NX+1;
        labelVertex[4] = labelVertex[0] + (NY+1)*(NX+1);
        labelVertex[5] = labelVertex[4] + 1;
        labelVertex[6] = labelVertex[4] + NX+2;
        labelVertex[7] = labelVertex[4] + NX+1;
        for (l=0; l<8; ++l)  LblEmNd[le][l] = labelVertex[l];
        le++;
      }
    }
  }

  // Generating info of face-vs-(4)nodes
  LblFcNd = new int[NumFcs][4];
  
  // X-faces (faces perpendicular to X-axis)
  lc = 0;
  for (i=0; i<=NX; ++i) {
    for (k=1; k<=NZ; ++k) {
      for (j=1; j<=NY; ++j) {
        labelVertex[0] = (k-1)*(NY+1)*(NX+1) + (j-1)*(NX+1) + (i+1);
        labelVertex[1] = labelVertex[0] + (NX+1);
        labelVertex[2] = labelVertex[1] + (NX+1)*(NY+1);
        labelVertex[3] = labelVertex[0] + (NX+1)*(NY+1);
        for (l=0; l<4; ++l)  LblFcNd[lc][l] = labelVertex[l];
        lc++;
      }
    }
  }
  
  // Y-faces (faces perpendicular to Y-axis)
  lc = (NX+1)*NY*NZ;
  for (j=0; j<=NY; ++j) {
    for (k=1; k<=NZ; ++k) {
      for (i=1; i<=NX; ++i) {
        labelVertex[0] = (k-1)*(NY+1)*(NX+1) + j*(NX+1) + i;
        labelVertex[1] = labelVertex[0] + 1;
        labelVertex[2] = labelVertex[1] + (NX+1)*(NY+1);
        labelVertex[3] = labelVertex[0] + (NX+1)*(NY+1);
        for (l=0; l<4; ++l)  LblFcNd[lc][l] = labelVertex[l];
        lc++;
      }
    }
  }
  
  // Z-faces (faces perpendicular to Z-axis)
  lc = (NX+1)*NY*NZ + NX*(NY+1)*NZ;
  for (k=0; k<=NZ; ++k) {
    for (j=1; j<=NY; ++j) {
      for (i=1; i<=NX; ++i) {
        labelVertex[0] = k*(NY+1)*(NX+1) + (j-1)*(NX+1) + i;
        labelVertex[1] = labelVertex[0] + 1;
        labelVertex[2] = labelVertex[1] + (NX+1);
        labelVertex[3] = labelVertex[0] + (NX+1);
        for (l=0; l<4; ++l)  LblFcNd[lc][l] = labelVertex[l];
        lc++;
      }
    }
  }
  
  // Generating info about boundary faces
  IsBndryFc = new int[NumFcs];
  for (lc=0; lc<NumFcs; ++lc)  IsBndryFc[lc] = 0;
  // X-faces (faces perpendicular to X-axis)
  // i = 0;
  lc = 0;
  for (k=1; k<=NZ; ++k) {
    for (j=1; j<=NY; ++j) {
      IsBndryFc[lc] = 1;
      lc++;
    }
  }
  // i = NX;
  lc = NX*NY*NZ;
  for (k=1; k<=NZ; ++k) {
    for (j=1; j<=NY; ++j) {
      IsBndryFc[lc] = 2;
      lc++;
    }
  }
  // Y-faces (faces perpendicular to Y-axis)
  // j = 0;
  lc = (NX+1)*NY*NZ;
  for (k=1; k<=NZ; ++k) {
    for (i=1; i<=NX; ++i) {
      IsBndryFc[lc] = 3;
      lc++;
    }
  }
  // j = NY;
  lc = (NX+1)*NY*NZ + NX*NY*NZ;
  for (k=1; k<=NZ; ++k) {
    for (i=1; i<=NX; ++i) {
      IsBndryFc[lc] = 4;
      lc++;
    }
  }
  // Z-faces (faces perpendicular to Z-axis)
  // k = 0;
  lc = (NX+1)*NY*NZ + NX*(NY+1)*NZ;
  for (j=1; j<=NY; ++j) {
    for (i=1; i<=NX; ++i) {
      IsBndryFc[lc] = 5;
      lc++;
    }
  }
  // k = NZ;
  lc = (NX+1)*NY*NZ + NX*(NY+1)*NZ + NX*NY*NZ;
  for (j=1; j<=NY; ++j) {
    for (i=1; i<=NX; ++i) {
      IsBndryFc[lc] = 6;
      lc++;
    }
  }

  // Generating info of element-vs-(6)faces
  LblEmFc = new int[NumEms][6];
  le = 0;
  for (k=1; k<=NZ; ++k) {
    for (j=1; j<=NY; ++j) {
      for (i=1; i<=NX; ++i) {
        // X-faces: left, right
        LblEmFc[le][0] = (i-1)*NY*NZ + (k-1)*NY + j;
        LblEmFc[le][1] =     i*NY*NZ + (k-1)*NY + j;
        // Y-faces: back, front
        LblEmFc[le][2] = (j-1)*NX*NZ + (k-1)*NX + i + (NX+1)*NY*NZ;
        LblEmFc[le][3] =     j*NX*NZ + (k-1)*NX + i + (NX+1)*NY*NZ;
        // Z-faces: bottom, top
        LblEmFc[le][4] = (k-1)*NX*NY + (j-1)*NX + i + (NX+1)*NY*NZ + NX*(NY+1)*NZ;
        LblEmFc[le][5] =     k*NX*NY + (j-1)*NX + i + (NX+1)*NY*NZ + NX*(NY+1)*NZ;
        le++;
      }
    }
  }
  
  // Generating info of face-vs-(2)elements based on LblEmFc
  LblFcEm = new int[NumFcs][2];
  int* cntEmsFc = new int[NumFcs];
  for (lc=0; lc<NumFcs; ++lc) {
    LblFcEm[lc][0] = 0;
    LblFcEm[lc][1] = 0;
    cntEmsFc[lc] = 0;
  }
  for (le=0; le<NumEms; ++le) {
    for (l=0; l<6; ++l) {
      lc = LblEmFc[le][l] - BgnLblFc;
      LblFcEm[lc][cntEmsFc[lc]] = le + BgnLblEm;
      cntEmsFc[lc] = cntEmsFc[lc] + 1;
    }
  }
  delete[] cntEmsFc;

  // Generating info of number of face neighbors
  /*
  NumFcNbrs = new int[NumFcs];
  for (lc=0; lc<NumFcs; ++lc) {
    NumFcNbrs[lc] = 11;
    if (IsBndryFc[lc]>0)  NumFcNbrs[lc] = 6;
  }
  */

  // Generating info of face-vs-(6/11)faces based on LblEmFc
  LblFcFc = new int[NumFcs][11];
  // Setting to zero
  for (lc=0; lc<NumFcs; ++lc)
    for (l=0; l<11; ++l)
      LblFcFc[lc][l] = 0;
  for (lc=0; lc<NumFcs; ++lc) {
    labelElementA = LblFcEm[lc][0];
    labelElementB = LblFcEm[lc][1];
    for (l=0; l<6; ++l) {
      LblFcFc[lc][l] = LblEmFc[labelElementA-BgnLblEm][l];
    }
    if (labelElementB>0) {
      j = 0;
      for (l=0; l<6; ++l) {
        k = LblEmFc[labelElementB-BgnLblEm][l];
        if ((k-1)==lc)  continue;
        LblFcFc[lc][6+j] = k;
        j++;
      }
    }
  }

  // JL20160727: NOT NEEDED FOR NOW
  // Generating info of "which face" based on LblFcEm, LblEmFc
  /*
  WchFc = new int[NumFcs][2];
  for (lc=0; lc<NumFcs; ++lc) {
    WchFc[lc][0] = 0;
    WchFc[lc][1] = 0;
  }
  for (lc=0; lc<NumFcs; ++lc) {
    lea = LblFcEm[lc][0] - BgnLblEm;
    for (j=0; j<6; ++j) {
      if (LblEmFc[lea][j]==lc+BgnLblFc) {
        WchFc[lc][0] = j+1;
        break;
      }
    }
    if (IsBndryFc[lc]==0) {
      leb = LblFcEm[lc][1] - BgnLblEm;
      for (j=0; j<6; ++j) {
        if (LblEmFc[leb][j]==lc+BgnLblFc) {
          WchFc[lc][1] = j+1;
          break;
        }
      }
    }
  }
  */

  flag = 2;
}


// HexaMesh: Destructor

HexaMesh::~HexaMesh()
{

  NumNds = 0;  NumFcs = 0;  NumEms = 0;  NumBndryPcs = 0;
  BgnLblNd = 0;  BgnLblFc = 0;  BgnLblEm = 0;  BgnLblBndry = 0;
  IsBndryFc = 0;  LblFcNd = 0;  LblFcFc = 0;  LblFcEm = 0;
  LblEmNd = 0;  LblEmFc = 0;
  nd = 0;  OutNmlBndry = 0;
  flag = 0;
  NX = 0;  NY = 0;  NZ = 0;
}


// HexaMesh:

int HexaMesh::numberBoundaryFaces() const
{
  int numBndryFcs = 0;
  for (int lc=0; lc<NumFcs; ++lc) {
    if (IsBndryFc[lc])  numBndryFcs++;
  }
  return numBndryFcs;
}


// HexaMesh: get a face

Quadri3d HexaMesh::face(int labelc) const
{
  int lc = labelc - BgnLblFc;
  PtVec3d A = nd[LblFcNd[lc][0]-BgnLblNd];
  PtVec3d B = nd[LblFcNd[lc][1]-BgnLblNd];
  PtVec3d C = nd[LblFcNd[lc][2]-BgnLblNd];
  PtVec3d D = nd[LblFcNd[lc][3]-BgnLblNd];
  return Quadri3d(A,B,C,D);
}


// JL20160708: AWKWARD!!! TO BE REVISED
// HexaMesh: Get an element

Hexa HexaMesh::element(int labele) const
{
  PtVec3d P0, P1, P2, P3, P4, P5, P6, P7;
  int ie = labele - BgnLblEm;
  P0 = nd[LblEmNd[ie][0]-BgnLblNd];
  P1 = nd[LblEmNd[ie][1]-BgnLblNd];
  P2 = nd[LblEmNd[ie][2]-BgnLblNd];
  P3 = nd[LblEmNd[ie][3]-BgnLblNd];
  P4 = nd[LblEmNd[ie][4]-BgnLblNd];
  P5 = nd[LblEmNd[ie][5]-BgnLblNd];
  P6 = nd[LblEmNd[ie][6]-BgnLblNd];
  P7 = nd[LblEmNd[ie][7]-BgnLblNd];
  return Hexa(P0,P1,P2,P3,P4,P5,P6,P7);
}


// HexaMesh: Get face-vs-(4)nodes info

void HexaMesh::getFaceNode(int labelc, int labelVertex[4]) const
{
  int lc = labelc - BgnLblFc;
  for (int j=0; j<4; ++j)  labelVertex[j] = LblFcNd[lc][j];
  return;
}


// HexaMesh: Get face-vs-(11)faces (neighbors) info

void HexaMesh::getFaceFace(int labelc, int labelFaceNeighbor[11]) const
{
  int lc = labelc - BgnLblFc;
  for (int j=0; j<11; ++j) {
    labelFaceNeighbor[j] = LblFcFc[lc][j];
  }
  return;
}


// HexaMesh: Get face-vs-(2)elements info

void HexaMesh::getFaceElement(int labelc,
                              int &labelElementA, int &labelElementB) const
{
  int lc = labelc - BgnLblFc;
  labelElementA = LblFcEm[lc][0];
  labelElementB = LblFcEm[lc][1];
  return;
}


// HexaMesh: Get element-vs-(8)nodes info

void HexaMesh::getElementNode(int labele, int labelVertex[8]) const
{
  int ie = labele - BgnLblEm;
  for (int j=0; j<8; ++j)  labelVertex[j] = LblEmNd[ie][j];
  return;
}


// HexaMesh: Get element-vs-(6)faces info

void HexaMesh::getElementFace(int labele, int labelFace[6]) const
{
  int ie = labele - BgnLblEm;
  for (int j=0; j<6; ++j)  labelFace[j] = LblEmFc[ie][j];
  return;
}


// HexaMesh: diameter

double HexaMesh::diameter() const
{
  double diam = 0;
  for (int labele=1; labele<=NumEms; ++labele) {
    Hexa hexa = element(labele);
    double d = hexa.diameter();
    if (d>diam)  diam = d;
  }
  return diam;
}


void HexaMesh::save2file(const char *filename, const Vector& scalarData, const FullMatrix& velocity) const {
    std::ofstream ofs("/Users/yinglili/Desktop/Numerical Solution of Velocity u (n=16).vtk");
    if (!ofs) {
        std::cerr << "Open data file failed. Exit." << std::endl;
        exit(-1);
    }

    // Write VTK header
    ofs << "# vtk DataFile Version 3.0" << std::endl;
    ofs << "HexaMesh data and Scalar pressure" << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write points (nodes)
    ofs << "POINTS " << NumNds << " float" << std::endl;
    for (int id = 0; id < NumNds; ++id) {
        ofs << nd[id].xCrd() << " " << nd[id].yCrd() << " " << nd[id].zCrd() << std::endl;
    }

    // Write cells (elements)
    ofs << "CELLS " << NumEms << " " << (NumEms * 9) << std::endl; // 9 includes the count
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << "8 "; // Number of points for each hexahedron
        for (int j = 0; j < 8; ++j) {
            ofs << LblEmNd[ie][j] - BgnLblNd << " "; // Adjust for 0-indexing
        }
        ofs << std::endl;
    }

    // Write cell types
    ofs << "CELL_TYPES " << NumEms << std::endl;
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << "12" << std::endl; // VTK_HEXAHEDRON
    }

    // Write cell data (scalar values)
    ofs << "CELL_DATA " << NumEms << std::endl;
    ofs << "SCALARS pressure float 1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << scalarData[ie] << std::endl; // Write the scalar data for each element
    }
    
    // Write point data for velocity
      ofs << "POINT_DATA " << NumNds << std::endl;  // Point data corresponds to nodes
      ofs << "VECTORS velocity float" << std::endl;  // Name of the vector
      for (int id = 0; id < NumNds; ++id) {
          PtVec3d pt = nd[id];  // Get the point for the current node
          PtVec3d velocityVec = fxnu(pt);  // Calculate velocity
          ofs << velocityVec.xCrd() << " " << velocityVec.yCrd() << " " << velocityVec.zCrd() << std::endl; // 3D velocity
      }

      // Write point data for pressure gradient
      ofs << "POINT_DATA " << NumNds << std::endl;  // Point data for nodes
      ofs << "VECTORS pressure_gradient float" << std::endl;  // Name of the vector
      for (int id = 0; id < NumNds; ++id) {
          PtVec3d pt = nd[id];  // Get the point for the current node
          PtVec3d pressureGradient = fxnpg(pt);  // Calculate pressure gradient
          ofs << pressureGradient.xCrd() << " " << pressureGradient.yCrd() << " " << pressureGradient.zCrd() << std::endl; // 3D pressure gradient
      }

      ofs.close();
    /*
    // Write point data for velocity
    ofs << "POINT_DATA " << NumNds << std::endl;  // Point data corresponds to nodes
    ofs << "VECTORS velocity float" << std::endl;  // Name of the vector
    for (int id = 0; id < NumNds; ++id) {
        ofs << velocity(id, 0) << " " << velocity(id, 1) << " " << velocity(id, 2) << std::endl; // 3D velocity
    }

    ofs.close();
*/
    
}


/*
void HexaMesh::save2file(const char *filename, const Vector& scalarData) const {
    
    std::ofstream ofs("/Users/yinglili/Desktop/Numerical_Solution of scalarData.vtk");
    
    if (!ofs) {
        std::cerr << "Open data file failed. Exit." << std::endl;
        exit(-1);
    }

    // Write VTK header
    ofs << "# vtk DataFile Version 3.0" << std::endl;
    ofs << "HexaMesh data and Scalar pressure" << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET STRUCTURED_GRID" << std::endl;
    
    // Write dimensions
    ofs << "DIMENSIONS " << (NX + 1) << " " << (NY + 1) << " " << (NZ + 1) << std::endl;


    // Write points (nodes)
    ofs << "POINTS " << NumNds << " float" << std::endl;
    for (int id = 0; id < NumNds; ++id) {
        ofs << nd[id].xCrd() << " " << nd[id].yCrd() << " " << nd[id].zCrd() << std::endl;
    }
    
    // Write cell types
    //ofs << "CELL_TYPES " << "12" << std::endl;
    ofs << "CELL_TYPES " << NumEms << std::endl;
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << "12" << std::endl; // VTK_HEXAHEDRON
    }


    // Write connectivity (element vs. nodes)
    ofs << "CELLS " << NumEms << " " << (NumEms * 9) << std::endl; // 9 includes the count
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << "8 "; // Number of points for each hexahedron
        for (int j = 0; j < 8; ++j) {
            ofs << LblEmNd[ie][j] - BgnLblNd << " "; // Adjust for 0-indexing
        }
        ofs << std::endl;
    }
   
    // Write cell data (scalar values)
    ofs << "CELL_DATA " << NumEms << std::endl;
    ofs << "SCALARS pressure float 1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for (int ie = 0; ie < NumEms; ++ie) {
        ofs << scalarData[ie] << std::endl; // Write the scalar data for each element
    }

    ofs.close();
}
*/
