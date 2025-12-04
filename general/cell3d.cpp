// cell3d.cpp
// Classes for 3-dim cells: Just geometric stuff:
//   bricks, tetrahedra, pyramids, prisms, hexahdra,
//   triangles-in-3d, rectangles-in-3d, etc.
// James Liu, ColoState; 2007/01--2018/01

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "vector.h"

#include "cell3d.h"
// #include "GaussQuad.h"  // Needed for Hexa::volume()
#include "mat3.h"
#include "PtVec3d.h"

using namespace std;


// The determinant related to any 3 given 3-dim points

double det3P(const PtVec3d P1, const PtVec3d P2, const PtVec3d P3)
{
  double x1 = P1.xCrd();
  double y1 = P1.yCrd();
  double z1 = P1.zCrd();

  double x2 = P2.xCrd();
  double y2 = P2.yCrd();
  double z2 = P2.zCrd();

  double x3 = P3.xCrd();
  double y3 = P3.yCrd();
  double z3 = P3.zCrd();

  double det = x1*y2*z3 + x2*y3*z1 + x3*y1*z2;
  det = det - x3*y2*z1 - x1*y3*z2 - x2*y1*z3;

  return det;
}


////////////////////////////////////////////////////////////////////////////////
// For class Brick
////////////////////////////////////////////////////////////////////////////////

// Brick: A constructor

Brick::Brick(PtVec3d P0, PtVec3d P1)
{
  vrtx[0] = P0;  vrtx[1] = P1;
}


// Brick: Another constructor

Brick::Brick(PtVec3d *Pt)
{
  vrtx[0] = Pt[0];  vrtx[1] = Pt[1];
}


// Brick: Copy constructor

Brick::Brick(const Brick &B)
{
  vrtx[0] = B.vrtx[0];  vrtx[1] = B.vrtx[1];
}


// Brick: Copy assignment

Brick &Brick::operator=(const Brick &B)
{
  if (this != &B) {
    vrtx[0] = B.vrtx[0];
    vrtx[1] = B.vrtx[1];
  }
  return *this;
}


// Brick: vertex

PtVec3d Brick::vertex(int i) const
{
  double x0, x1, y0, y1, z0, z1;
  x0 = vrtx[0].xCrd();  x1 = vrtx[1].xCrd();
  y0 = vrtx[0].yCrd();  y1 = vrtx[1].yCrd();
  z0 = vrtx[0].zCrd();  z1 = vrtx[1].zCrd();
  if (i==0)  return PtVec3d(x0,y0,z0);
  if (i==1)  return PtVec3d(x1,y0,z0);
  if (i==2)  return PtVec3d(x1,y1,z0);
  if (i==3)  return PtVec3d(x0,y1,z0);
  if (i==4)  return PtVec3d(x0,y0,z1);
  if (i==5)  return PtVec3d(x1,y0,z1);
  if (i==6)  return PtVec3d(x1,y1,z1);
  if (i==7)  return PtVec3d(x0,y1,z1);
  return PtVec3d(0,0,0);
}


// Brick: volume

double Brick::volume() const
{
  double x0 = vrtx[0].xCrd();
  double y0 = vrtx[0].yCrd();
  double z0 = vrtx[0].zCrd();

  double x1 = vrtx[1].xCrd();
  double y1 = vrtx[1].yCrd();
  double z1 = vrtx[1].zCrd();

  return (x1-x0)*(y1-y0)*(z1-z0);
}


// Brick: center

PtVec3d Brick::center() const
{
  return 0.5*(vrtx[0]+vrtx[1]);
}


////////////////////////////////////////////////////////////////////////////////
// For class Tri3d
////////////////////////////////////////////////////////////////////////////////

// Tri3d: A constructor

Tri3d::Tri3d(PtVec3d A, PtVec3d B, PtVec3d C)
{
  vrtx[0] = A;  vrtx[1] = B;  vrtx[2] = C;
}


// Tri3d: Another constructor

Tri3d::Tri3d(PtVec3d P[3])
{
  for (int i=0; i<3; ++i)  vrtx[i] = P[i];
}


// Tri3d: Copy constructor

Tri3d::Tri3d(const Tri3d &T)
{
  for (int i=0; i<3; ++i)  vrtx[i] = T.vrtx[i];
}


// Tri3d: Copy assignment

Tri3d &Tri3d::operator=(const Tri3d &T)
{
  if (this != &T) {
    for (int i=0; i<3; ++i)  vrtx[i] = T.vrtx[i];
  }
  return *this;
}


// Tri3d: area

double Tri3d::area() const
{
  // Approach 1
  double a = dist(vrtx[1],vrtx[2]);
  double b = dist(vrtx[0],vrtx[2]);
  double c = dist(vrtx[0],vrtx[1]);
  double s = 0.5*(a+b+c);
  return sqrt(s*(s-a)*(s-b)*(s-c));
  // Approach 2
  /*
  PtVec3d v = crossProd(vrtx[1]-vrtx[0], vrtx[2]-vrtx[0]);
  return 0.5*v.l2norm();
  */
}


// Tri3d: center

PtVec3d Tri3d::center() const
{
  return (1.0/3)*(vrtx[0]+vrtx[1]+vrtx[2]);
}


// Tri3d: diameter
/*
 double Tri3d::diameter() const
 {
 double diam;
 vector<double> edge(3);
 edge[0] = dist(P0, P1);
 edge[1] = dist(P1, P2);
 edge[2] = dist(P2, P0);
 diam = *(max_element(edge.begin(), edge.end()));
 return diam;
 }
 */


// Tri3d: (unit) normal (vector)

PtVec3d Tri3d::normal() const
{
  PtVec3d u = vrtx[1] - vrtx[0];
  PtVec3d v = vrtx[2] - vrtx[0];
  PtVec3d w = crossProduct(u, v);
  double a = w.l2norm();
  return (1/a)*w;
}


// Tri3d: The mapping from the reference triangle

PtVec3d Tri3d::mapping(double xhat, double yhat) const
{
  return (1-xhat-yhat)*vrtx[0] + xhat*vrtx[1] + yhat*vrtx[2];
}


////////////////////////////////////////////////////////////////////////////////
// For class Rect3d
////////////////////////////////////////////////////////////////////////////////

// Rect3d: A constructor

Rect3d::Rect3d(PtVec3d A, PtVec3d B, PtVec3d C, PtVec3d D)
{
  vrtx[0] = A;  vrtx[1] = B;  vrtx[2] = C;  vrtx[3] = D;
}


// Rect3d: Copy constructor

Rect3d::Rect3d(const Rect3d &R)
{
  for (int i=0; i<4; ++i)  vrtx[i] = R.vrtx[i];
}


// Rect3d: Copy assignment

Rect3d &Rect3d::operator=(const Rect3d &R)
{
  if (this != &R) {
    for (int i=0; i<4; ++i)  vrtx[i] = R.vrtx[i];
  }
  return *this;
}


// Rect3d:: area
// Assuming it is INDEED a rectangle

double Rect3d::area() const
{
  double a = dist(vrtx[0],vrtx[1]);
  double b = dist(vrtx[0],vrtx[3]);
  return a*b;
}


// Rect3d: center
// Assuming it is INDEED a rectangle

PtVec3d Rect3d::center() const
{
  return 0.5*(vrtx[0]+vrtx[2]);
}


////////////////////////////////////////////////////////////////////////////////
// For class Quadri3d
////////////////////////////////////////////////////////////////////////////////

// Quadri3d: A constructor

Quadri3d::Quadri3d(PtVec3d A, PtVec3d B, PtVec3d C, PtVec3d D)
{
  vrtx[0] = A;  vrtx[1] = B;  vrtx[2] = C;  vrtx[3] = D;
}


// Quadri3d: Another constructor

Quadri3d::Quadri3d(PtVec3d P[4])
{
  for (int i=0; i<4; ++i)  vrtx[i] = P[i];
}


// Quadri3d: Yet another constructor
/*
Quadri3d::Quadri3d(PtVec3d *Pt)
{
  for (int i=0; i<4; ++i)  vrtx[i] = Pt[i];
}
*/


// Quadri3d: Copy constructor

Quadri3d::Quadri3d(const Quadri3d &Q)
{
  for (int i=0; i<4; ++i)  vrtx[i] = Q.vrtx[i];
}


// Quadri3d: Copy assignment

Quadri3d &Quadri3d::operator=(const Quadri3d &Q)
{
  if (this != &Q) {
    for (int i=0; i<4; ++i)  vrtx[i] = Q.vrtx[i];
  }
  return *this;
}


// JL20150331: TO BE FINISHED
// Quadri3d:: area
/*
 double Quadri3d::area() const
 {
 return 1.0;
 }
 */


// Quadri3d: center

PtVec3d Quadri3d::center() const
{
  PtVec3d a = vrtx[1] - vrtx[0];
  PtVec3d b = vrtx[3] - vrtx[0];
  PtVec3d c = (vrtx[2]-vrtx[0]) - (a+b);
  return vrtx[0] + 0.5*a + 0.5*b + 0.25*c;
  // JL20150425: The following is also correct
  // return 0.25*(vrtx[0]+vrtx[1]+vrtx[2]+vrtx[3]);  // 0.25=1/4, average
}


// JL20171104: For a flat quadrilateral only
// Quadri3d: area

double Quadri3d::area() const
{
  PtVec3d tmp1 = crossProduct(vrtx[0]-vrtx[1],vrtx[0]-vrtx[2]);
  PtVec3d tmp2 = crossProduct(vrtx[0]-vrtx[2],vrtx[0]-vrtx[3]);
  return 0.5*tmp1.l2norm() + 0.5*tmp2.l2norm();
}


// JL20160723: TO BE FINISHED LATER
// Quadri3d: area, general case, need quadrilature
/*
double Quadri3d::area(const GaussQuad &GQQ) const
*/


// Quadri3d: bilinearmapping from [0,1]^2 for (\hat{x},\hat{y})

PtVec3d Quadri3d::bilinearmapping(double xhat, double yhat) const
{
  PtVec3d a = vrtx[1] - vrtx[0];
  PtVec3d b = vrtx[3] - vrtx[0];
  PtVec3d c = (vrtx[2]-vrtx[0]) - (a+b);
  return vrtx[0] + xhat*a + yhat*b + xhat*yhat*c;
}


// JL20150409: NEEDS DOUBLE-CHECK!!
// Quadri3d: Jacobian determinant for (xhat,yhat)

double Quadri3d::JacobianDeterminant(double xhat, double yhat) const
{
  PtVec3d a = vrtx[1] - vrtx[0];
  PtVec3d b = vrtx[3] - vrtx[0];
  PtVec3d c = (vrtx[2]-vrtx[0]) - (a+b);
  PtVec3d w = crossProduct(a+yhat*c, b+xhat*c);
  return w.l2norm();
}


// Quadri3d: (unit) normal (vector)

PtVec3d Quadri3d::normal(double xhat, double yhat) const
{
  PtVec3d a = vrtx[1] - vrtx[0];
  PtVec3d b = vrtx[3] - vrtx[0];
  PtVec3d c = (vrtx[2]-vrtx[0]) - (a+b);
  PtVec3d w = crossProduct(a+yhat*c, b+xhat*c);
  double s = w.l2norm();
  // if (s==0) std::cout << "s=" << s << "\n";
  return (1.0/s)*w;
}


////////////////////////////////////////////////////////////////////////////////
// For class Tetra
////////////////////////////////////////////////////////////////////////////////

// Tetra: A constructor

Tetra::Tetra(PtVec3d P0, PtVec3d P1, PtVec3d P2, PtVec3d P3)
{
  vrtx[0] = P0;  vrtx[1] = P1;  vrtx[2] = P2;  vrtx[3] = P3;
}


// Tetra: Another constructor

Tetra::Tetra(PtVec3d *Pt)
{
  const int numVrtz = 4;
  for (int i=0; i<numVrtz; ++i)  vrtx[i] = Pt[i];
}


// Tetra: Copy constructor

Tetra::Tetra(const Tetra &Te)
{
  const int numVrtz = 4;
  for (int i=0; i<numVrtz; ++i)  vrtx[i] = Te.vrtx[i];
}


// Tetra: Copy assignment

Tetra &Tetra::operator=(const Tetra &Te)
{
  if (this != &Te) {
    const int numVrtz = 4;
    for (int i=0; i<numVrtz; ++i)  vrtx[i] = Te.vrtx[i];
  }
  return *this;
}


// Tetra: enriching

void Tetra::enrich()
{
  double x[4], y[4], z[4];
  PtVec3d P[4], U[4], V[4], W[4];
  
  // The determinant for the tetra
  double dd = 6 * volume();
  
  // Primary points/vectors
  for (int i=0; i<=3; ++i) {
    P[i] = vrtx[i];
    x[i] = P[i].xCrd();
    y[i] = P[i].yCrd();
    z[i] = P[i].zCrd();
  }
  
  // Auxiliary points/vectors
  for (int i=0; i<4; ++i) {
    U[i] = PtVec3d(1,y[i],z[i]);
    V[i] = PtVec3d(1,x[i],z[i]);
    W[i] = PtVec3d(1,x[i],y[i]);
  }
  
  // For coefficients
  cof[0][0] =  det3P(P[1],P[2],P[3])/dd;
  cof[0][1] = -det3P(U[1],U[2],U[3])/dd;
  cof[0][2] =  det3P(V[1],V[2],V[3])/dd;
  cof[0][3] = -det3P(W[1],W[2],W[3])/dd;
  
  // For coefficients
  cof[1][0] = -det3P(P[0],P[2],P[3])/dd;
  cof[1][1] =  det3P(U[0],U[2],U[3])/dd;
  cof[1][2] = -det3P(V[0],V[2],V[3])/dd;
  cof[1][3] =  det3P(W[0],W[2],W[3])/dd;
  
  // For coefficients
  cof[2][0] =  det3P(P[0],P[1],P[3])/dd;
  cof[2][1] = -det3P(U[0],U[1],U[3])/dd;
  cof[2][2] =  det3P(V[0],V[1],V[3])/dd;
  cof[2][3] = -det3P(W[0],W[1],W[3])/dd;
  
  // For coefficients
  cof[3][0] = -det3P(P[0],P[1],P[2])/dd;
  cof[3][1] =  det3P(U[0],U[1],U[2])/dd;
  cof[3][2] = -det3P(V[0],V[1],V[2])/dd;
  cof[3][3] =  det3P(W[0],W[1],W[2])/dd;
}


// Tetra: volume

double Tetra::volume() const
{
  double d0 = det3P(vrtx[1], vrtx[2], vrtx[3]);
  double d1 = det3P(vrtx[0], vrtx[2], vrtx[3]);
  double d2 = det3P(vrtx[0], vrtx[1], vrtx[3]);
  double d3 = det3P(vrtx[0], vrtx[1], vrtx[2]);
  return fabs(d0-d1+d2-d3)/6.0;
}


// Tetra: diameter

double Tetra::diameter() const
{
  PtVec3d edge;
  Vector edgelen(6);
  
  edge = vrtx[0] - vrtx[1];
  edgelen(1) = edge.l2norm();

  edge = vrtx[0] - vrtx[2];
  edgelen(2) = edge.l2norm();

  edge = vrtx[0] - vrtx[3];
  edgelen(3) = edge.l2norm();

  edge = vrtx[1] - vrtx[2];
  edgelen(4) = edge.l2norm();
  
  edge = vrtx[1] - vrtx[3];
  edgelen(5) = edge.l2norm();

  edge = vrtx[2] - vrtx[3];
  edgelen(6) = edge.l2norm();
  
  return edgelen.l0norm();
}


// JL20150428: CODE COULD BE OPTIMIZED, TO BE CHECKED LATER!
// Tetra: the outward unit normal vector for the face opp. to vertex i(=0,1,2,3)

PtVec3d Tetra::normal(int i) const
{
  double a = 1;
  PtVec3d u, v, w;
  switch (i) {
    case 0:
      u = vrtx[2] - vrtx[1];
      v = vrtx[3] - vrtx[1];
      break;
    case 1:
      u = vrtx[3] - vrtx[0];
      v = vrtx[2] - vrtx[0];
      break;
    case 2:
      u = vrtx[1] - vrtx[0];
      v = vrtx[3] - vrtx[0];
      break;
    case 3:
      u = vrtx[2] - vrtx[0];
      v = vrtx[1] - vrtx[0];
      break;
  }
  w = crossProduct(u,v);
  a = w.l2norm();
  return (1/a)*w;
}


// Tetra: the face oppoiste to vertex i (0<=i<3)

Tri3d Tetra::face(int i) const
{
  PtVec3d A, B, C;
  switch (i) {
    case 0:
      A = vrtx[1];  B = vrtx[2];  C = vrtx[3];
      break;
    case 1:
      A = vrtx[0];  B = vrtx[2];  C = vrtx[3];
      break;
    case 2:
      A = vrtx[0];  B = vrtx[1];  C = vrtx[3];
      break;
    case 3:
      A = vrtx[0];  B = vrtx[1];  C = vrtx[2];
      break;
  }
  return Tri3d(A,B,C);
}


// Tetra:

void Tetra::getMaxMinCrds(double &xmax, double &xmin,
                          double &ymax, double &ymin,
                          double &zmax, double &zmin) const
{
  const int numVrtz = 4;
  double x[numVrtz], y[numVrtz], z[numVrtz];
  
  for (int i=0; i<numVrtz; ++i) {
    x[i] = vrtx[i].xCrd();
    y[i] = vrtx[i].yCrd();
    z[i] = vrtx[i].zCrd();
  }
  
  xmax = *(max_element(x, x + numVrtz));
  ymax = *(max_element(y, y + numVrtz));
  zmax = *(max_element(z, z + numVrtz));
  
  xmin = *(min_element(x, x + numVrtz));
  ymin = *(min_element(y, y + numVrtz));
  zmin = *(min_element(z, z + numVrtz));
}


// Tetra:

void Tetra::getMaxMinEdges(double &maxEdge, double &minEdge)
{
  vector<double> edge(6);
  edge[0] = dist(vrtx[0], vrtx[1]);
  edge[1] = dist(vrtx[0], vrtx[2]);
  edge[2] = dist(vrtx[1], vrtx[2]);
  edge[3] = dist(vrtx[0], vrtx[3]);
  edge[4] = dist(vrtx[1], vrtx[3]);
  edge[5] = dist(vrtx[2], vrtx[3]);
  maxEdge = *(max_element(edge.begin(), edge.end()));
  minEdge = *(min_element(edge.begin(), edge.end()));
  return;
}


/*************************************************************
 // Tetra: the mapping FROM the reference tetrahedron
 
 PtVec3d Tetra::mapping(double a, double b, double c) const
 {
 double x[4], y[4], z[4], xx, yy, zz;
 
 for (int i=0; i<4; ++i) {
 x[i] = vrtx[i].xCrd();
 y[i] = vrtx[i].yCrd();
 z[i] = vrtx[i].zCrd();
 }
 
 xx = (x[1]-x[0])*a + (x[2]-x[0])*b + (x[3]-x[0])*c + x[0];
 yy = (y[1]-y[0])*a + (y[2]-y[0])*b + (y[3]-y[0])*c + y[0];
 zz = (z[1]-z[0])*a + (z[2]-z[0])*b + (z[3]-z[0])*c + z[0];
 
 return PtVec3d(xx, yy, zz);
 }
 
 
 // Tetra: the inverse mapping TO the reference tetrahedron
 
 PtVec3d Tetra::invMapping(PtVec3d P) const
 {
 double D, D1, D2, D3;
 
 D = det3P(vrtx[1]-vrtx[0], vrtx[2]-vrtx[0], vrtx[3]-vrtx[0]);
 D1 = det3P(P-vrtx[0], vrtx[2]-vrtx[0], vrtx[3]-vrtx[0]);
 D2 = det3P(vrtx[1]-vrtx[0], P-vrtx[0], vrtx[3]-vrtx[0]);
 D3 = det3P(vrtx[1]-vrtx[0], vrtx[2]-vrtx[0], P-vrtx[0]);
 
 if (D==0)  {
 cout << "Determinant 0!\n";
 exit(EXIT_FAILURE);
 }
 
 return PtVec3d(D1/D, D2/D, D3/D);
 }
 
 
 int Tetra::isInTetra(PtVec3d P) const
 {
 int yn = 0;
 
 PtVec3d Q = invMapping(P);
 double a = Q.xCrd();
 double b = Q.yCrd();
 double c = Q.zCrd();
 
 if (a>=0 && b>=0 && c>=0 && (a+b+c)<=1)  yn = 1;
 
 return yn;
 }
 
 
 void Tetra::getBaryCrds(double *baryCrd, PtVec3d P) const
 {
 Tetra T0(P, vrtx[1], vrtx[2], vrtx[3]);
 Tetra T1(vrtx[0], P, vrtx[2], vrtx[3]);
 Tetra T2(vrtx[0], vrtx[1], P, vrtx[3]);
 Tetra T3(vrtx[0], vrtx[1], vrtx[2], P);
 
 baryCrd[0] = T0.volume()/volume();
 baryCrd[1] = T1.volume()/volume();
 baryCrd[2] = T2.volume()/volume();
 baryCrd[3] = T3.volume()/volume();
 
 return;
 }
 
 
 Mat3 Tetra::invJacobian() const
 {
 Mat3 Jac('c', vrtx[1]-vrtx[0], vrtx[2]-vrtx[0], vrtx[3]-vrtx[0]);
 return Jac.inverse();
 }
*************************************************************/


////////////////////////////////////////////////////////////////////////////////
// For class Hexa
////////////////////////////////////////////////////////////////////////////////

// Hexa: Default constructor

Hexa::Hexa()
{
  const int numVrtz = 8;
  for (int i=0; i<numVrtz; ++i)  vrtx[i] = PtVec3d(0,0,0);
}


// Hexa: A constructor

Hexa::Hexa(PtVec3d P0, PtVec3d P1, PtVec3d P2, PtVec3d P3,
           PtVec3d P4, PtVec3d P5, PtVec3d P6, PtVec3d P7)
{
  vrtx[0] = P0;  vrtx[1] = P1;  vrtx[2] = P2;  vrtx[3] = P3;
  vrtx[4] = P4;  vrtx[5] = P5;  vrtx[6] = P6;  vrtx[7] = P7;
}


// Hexa: Yet another constructor

Hexa::Hexa(PtVec3d P[8])
{
  vrtx[0] = P[0];  vrtx[1] = P[1];  vrtx[2] = P[2];  vrtx[3] = P[3];
  vrtx[4] = P[4];  vrtx[5] = P[5];  vrtx[6] = P[6];  vrtx[7] = P[7];
}


// Hexa: Copy constructor

Hexa::Hexa(const Hexa &H)
{
  const int numVrtz = 8;
  for (int i=0; i<numVrtz; ++i)  vrtx[i] = H.vrtx[i];
}


// Hexa: Copy assignment

Hexa &Hexa::operator=(const Hexa &H)
{
  if (this != &H) {
    const int numVrtz = 8;
    for (int i=0; i<numVrtz; ++i)  vrtx[i] = H.vrtx[i];
  }
  return *this;
}


// Hexa: enriching...

void Hexa::enrich()
{
  // vectors for edges, faces, volume
  veca = vrtx[1] - vrtx[0];  // For xhat
  vecb = vrtx[3] - vrtx[0];  // For yhat
  vecc = vrtx[4] - vrtx[0];  // For zhat
  vecd = (vrtx[7] - vrtx[0]) - (vecb + vecc);  // No xhat
  vece = (vrtx[5] - vrtx[0]) - (vecc + veca);  // No yhat
  vecf = (vrtx[2] - vrtx[0]) - (veca + vecb);  // No zhat
  // For all xhat, yhat, zhat: deviation from being a parallelopiped
  vecg = (vrtx[6]-vrtx[0]) - ((veca+vecb+vecc)+(vecd+vece+vecf));
  return;
}


// Hexa: center

PtVec3d Hexa::center() const
{
  PtVec3d cntr(0,0,0);
  for (int i=0; i<8; ++i)  cntr = cntr + vrtx[i];
  cntr = 0.125*cntr;  // 0.125 = 1/8
  return cntr;
  // JL20150604: WHICH ONE IS CORRECT?
  // return trilinearmapping(0.5, 0.5, 0.5);
}


// Hexa: diameter
// Assuming the 8 vertices are correctly oriented

double Hexa::diameter() const
{
  PtVec3d diag;
  Vector diaglen(4);

  diag = vrtx[0] - vrtx[6];
  diaglen[0] = diag.l2norm();
  
  diag = vrtx[1] - vrtx[7];
  diaglen[1] = diag.l2norm();
  
  diag = vrtx[2] - vrtx[4];
  diaglen[2] = diag.l2norm();

  diag = vrtx[3] - vrtx[5];
  diaglen[3] = diag.l2norm();

  return diaglen.l0norm();
}


// Hexa: volume
// Assuming the 8 vertices are correctly oriented
/*
double Hexa::volume(const GaussQuad &GQH) const
{
  double vol = 0;
  for (int k=0; k<GQH.numberQuadraturePoints(); ++k) {
    double xhat = GQH.CartesianCoordinate(k,0);
    double yhat = GQH.CartesianCoordinate(k,1);
    double zhat = GQH.CartesianCoordinate(k,2);
    double jac = hexa.JacobianDeterminant(xhat, yhat, zhat);
    jac = fabs(jac);  // JL20150422: DON'T USE fabs?
    vol = vol + GQH.weight(k) * jac;
  }
}
*/


/*
 void Hexa::getMaxMinCrds(double &xmax, double &xmin,
 double &ymax, double &ymin, double &zmax, double &zmin)
 {
 const int numVrtz = 8;
 double x[numVrtz], y[numVrtz], z[numVrtz];
 
 for (int i=0; i<numVrtz; ++i) {
 x[i] = vrtx[i].xCrd();
 y[i] = vrtx[i].yCrd();
 z[i] = vrtx[i].zCrd();
 }
 
 xmax = *(max_element(x, x + numVrtz));
 ymax = *(max_element(y, y + numVrtz));
 zmax = *(max_element(z, z + numVrtz));
 
 xmin = *(min_element(x, x + numVrtz));
 ymin = *(min_element(y, y + numVrtz));
 zmin = *(min_element(z, z + numVrtz));
 }
 */


/*
 void Hexa::getMaxMinEdges(double &maxEdge, double &minEdge)
 {
 vector<double> edge(12);
 edge[0] = dist(vrtx[0], vrtx[1]);
 edge[1] = dist(vrtx[1], vrtx[2]);
 edge[2] = dist(vrtx[2], vrtx[3]);
 edge[3] = dist(vrtx[0], vrtx[3]);
 edge[4] = dist(vrtx[4], vrtx[5]);
 edge[5] = dist(vrtx[5], vrtx[6]);
 edge[6] = dist(vrtx[6], vrtx[7]);
 edge[7] = dist(vrtx[4], vrtx[7]);
 edge[8] = dist(vrtx[0], vrtx[4]);
 edge[9] = dist(vrtx[1], vrtx[5]);
 edge[10] = dist(vrtx[2], vrtx[6]);
 edge[11] = dist(vrtx[3], vrtx[7]);
 maxEdge = *(max_element(edge.begin(), edge.end()));
 minEdge = *(min_element(edge.begin(), edge.end()));
 return;
 }
 */


// JL20150404: UNFORTUNATELY UNTRUE FOR GENERAL HEXAHEDRA !!
// Hexa: Volume as the union of 5 tetrahedra
/*
 double Hexa::volume5() const
 {
 double vol = 0.0;
 vector<Tetra> HT = decmp5();
 for (int i=0; i<HT.size(); ++i)  vol += HT[i].volume();
 return vol;
 }
 */


// JL20150404: UNFORTUNATELY UNTRUE FOR GENERAL HEXAHEDRA !!
// Hexa: Volume as the union of 6 tetrahedra
/*
 double Hexa::volume6() const
 {
 double vol = 0.0;
 vector<Tetra> HT = decmp6();
 for (int i=0; i<HT.size(); ++i)  vol += HT[i].volume();
 return vol;
 }
 */


// JL20150404: UNFORTUNATELY UNTRUE FOR GENERAL HEXAHEDRA !!
// Hexa: Decomposition into 5 tetrahedra
/*
 vector<Tetra> Hexa::decmp5() const
 {
 vector<Tetra> HT(5);
 HT[0] = Tetra(vrtx[0], vrtx[1], vrtx[3], vrtx[4]);
 HT[1] = Tetra(vrtx[1], vrtx[2], vrtx[3], vrtx[6]);
 HT[2] = Tetra(vrtx[1], vrtx[3], vrtx[4], vrtx[6]);
 HT[3] = Tetra(vrtx[1], vrtx[4], vrtx[5], vrtx[6]);
 HT[4] = Tetra(vrtx[3], vrtx[4], vrtx[6], vrtx[7]);
 return HT;
 }
 */


// JL20150404: UNFORTUNATELY UNTRUE FOR GENERAL HEXAHEDRA !!
// Hexa: Decomposition into 6 tetrahedra
/*
 vector<Tetra> Hexa::decmp6() const
 {
 vector<Tetra> HT(6);
 HT[0] = Tetra(vrtx[0], vrtx[1], vrtx[2], vrtx[4]);
 HT[1] = Tetra(vrtx[1], vrtx[2], vrtx[4], vrtx[5]);
 HT[2] = Tetra(vrtx[2], vrtx[4], vrtx[5], vrtx[6]);
 HT[3] = Tetra(vrtx[0], vrtx[2], vrtx[3], vrtx[4]);
 HT[4] = Tetra(vrtx[2], vrtx[3], vrtx[4], vrtx[6]);
 HT[5] = Tetra(vrtx[3], vrtx[4], vrtx[6], vrtx[7]);
 return HT;
 }
 */


/*
 int Hexa::isInHexa(PtVec3d P)
 {
 int j=0;
 vector<Tetra> HT = decmp5();
 for (int i=0; i<HT.size(); ++i) {
 if (HT[i].isInTetra(P)) {
 j = i+1;
 break;
 }
 }
 return j;
 }
 */


// Hexa: getting one of the 6 faces: 0<=i<=5

Quadri3d Hexa::face(int i) const
{
  PtVec3d A, B, C, D;
  switch (i) {
    case 0:  // x-left face
      A = vrtx[0];  B = vrtx[4];  C = vrtx[7];  D = vrtx[3];
      break;
    case 1:  // x-right face
      A = vrtx[1];  B = vrtx[2];  C = vrtx[6];  D = vrtx[5];
      break;
    case 2:  // y-back face
      A = vrtx[0];  B = vrtx[1];  C = vrtx[5];  D = vrtx[4];
      break;
    case 3:  // y-front face
      A = vrtx[2];  B = vrtx[3];  C = vrtx[7];  D = vrtx[6];
      break;
    case 4:  // z-bottom face
      A = vrtx[0];  B = vrtx[3];  C = vrtx[2];  D = vrtx[1];
      break;
    case 5:  // z-top face
      A = vrtx[4];  B = vrtx[5];  C = vrtx[6];  D = vrtx[7];
      break;
    default:
      A = PtVec3d(0,0,0);  B = A;  C = A;  D = A;
      break;
  }
  return Quadri3d(A,B,C,D);
}


// Hexa: trilinear mapping from [0,1]^3 for (\hat{x},\hat{y},\hat{z})

PtVec3d Hexa::trilinearmapping(double xhat, double yhat, double zhat) const
{
  return vrtx[0] + xhat*veca + yhat*vecb + zhat*vecc
  + yhat*zhat*vecd + zhat*xhat*vece + xhat*yhat*vecf
  + xhat*yhat*zhat*vecg;
}


// Hexa: xhat covariant vector for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3
// i.e., the tangential vector of the trilinear mapping with respect to xhat

PtVec3d Hexa::covariantVector_xhat(double xhat, double yhat, double zhat) const
{
  return veca + yhat*vecf + zhat*vece + (yhat*zhat)*vecg;
}


// Hexa: yhat covariant vector for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3
// i.e., the tangential vector of the trilinear mapping with respect to yhat

PtVec3d Hexa::covariantVector_yhat(double xhat, double yhat, double zhat) const
{
  return vecb + zhat*vecd + xhat*vecf + (zhat*xhat)*vecg;
}


// Hexa: zhat covariant vector for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3
// i.e., the tangential vector of the trilinear mapping with respect to zhat

PtVec3d Hexa::covariantVector_zhat(double xhat, double yhat, double zhat) const
{
  return vecc + xhat*vece + yhat*vecd + (xhat*yhat)*vecg;
}


// Hexa: Jacobian determinant for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3

double Hexa::JacobianDeterminant(double xhat, double yhat, double zhat) const
{
  PtVec3d TanX = covariantVector_xhat(xhat, yhat, zhat);
  PtVec3d TanY = covariantVector_yhat(xhat, yhat, zhat);
  PtVec3d TanZ = covariantVector_zhat(xhat, yhat, zhat);
  return dotProduct(crossProduct(TanX,TanY),TanZ);
}


// Hexa: Jacobian matrix for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3

Mat3 Hexa::JacobianMatrix(double xhat, double yhat, double zhat) const
{
  PtVec3d TanX = covariantVector_xhat(xhat, yhat, zhat);
  PtVec3d TanY = covariantVector_yhat(xhat, yhat, zhat);
  PtVec3d TanZ = covariantVector_zhat(xhat, yhat, zhat);
  return Mat3('c', TanX, TanY, TanZ);
}


// Hexa: Piola matrix for (\hat{x},\hat{y},\hat{z}) \in [0,1]^3

Mat3 Hexa::PiolaMatrix(double xhat, double yhat, double zhat) const
{
  PtVec3d TanX = covariantVector_xhat(xhat, yhat, zhat);
  PtVec3d TanY = covariantVector_yhat(xhat, yhat, zhat);
  PtVec3d TanZ = covariantVector_zhat(xhat, yhat, zhat);
  double det_1 = 1.0/dotProduct(crossProduct(TanX,TanY),TanZ);
  return det_1 * Mat3('c',TanX,TanY,TanZ);
}


