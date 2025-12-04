// cell3d.h
// Classes for 3-dim cells: Just geometric stuff:
//   bricks, tetrahedra, pyramids, prisms, hexahdra,
//   triangles-in-3d, rectangles-in-3d, etc.
// James Liu, ColoState; 2007/01--2018/01

#ifndef Cell3d_H
#define Cell3d_H

#include "mat3.h"
// #include "GaussQuad.h"  // Needed for Quadri3d::area(), Hexa::volume()
#include "PtVec3d.h"

// enum Cell3d {brick=3, tetra, prism, pyra, hexa, tri3d, rect3d};

double det3P(const PtVec3d P1, const PtVec3d P2, const PtVec3d P3);


////////////////////////////////////////////////////////////////////////////////

// Bricks

class Brick {
private:
  PtVec3d vrtx[2];  // Lower-back-left (x0,y0,z0), Upper-front-right (x1,y1,z1)
public:
  Brick() {;}  // Default constructor
  Brick(PtVec3d P0, PtVec3d P1);  // A constructor
  Brick(PtVec3d *Pt);  // Another constructor
  Brick(const Brick &B);  // Copy constructor
  ~Brick() {;}  // Destructor
  Brick &operator=(const Brick &B);  // Copy assignment
  PtVec3d vertex(int i) const;  // 0<=i<8
  double volume() const;
  PtVec3d center() const;
  void getVertices(PtVec3d &P0, PtVec3d &P1) const {P0=vrtx[0]; P1=vrtx[1];}
};

////////////////////////////////////////////////////////////////////////////////

// Triangles in 3-dim

class Tri3d {
private:
  PtVec3d vrtx[3];
public:
  Tri3d() {;}  // Default constructor
  Tri3d(PtVec3d A, PtVec3d B, PtVec3d C);  // A constructor
  Tri3d(PtVec3d P[3]);  // Another constructor
  Tri3d(const Tri3d &T);  // Copy constructor
  ~Tri3d() {;}  // Destructor
  Tri3d &operator=(const Tri3d &T);  // Copy assignment
  PtVec3d vertex(int j) const {return vrtx[j];}  // 0<=j<3
  double area() const;
  PtVec3d center() const;
  // double diameter() const;
  PtVec3d normal() const;
  PtVec3d mapping(double xhat, double yhat) const;
};


// Rectangles in 3-dim

class Rect3d {
private:
  PtVec3d vrtx[4];
public:
  Rect3d() {;}  // Default constructor
  Rect3d(PtVec3d A, PtVec3d B, PtVec3d C, PtVec3d D);  // A constructor
  Rect3d(const Rect3d &R);  // Copy constructor
  ~Rect3d() {;}  // Destructor
  Rect3d &operator=(const Rect3d &R);  // Copy assignment
  double area() const;
  PtVec3d center() const;
  // double diameter() const;
  // PtVec3d normal() const;
  // PtVec3d mapping(double xhat, double yhat) const;
};


// Quadrilaterals in 3-dim

class Quadri3d {
private:
  PtVec3d vrtx[4];
  PtVec3d veca, vecb, vecc;
public:
  Quadri3d() {;}  // Default constructor
  Quadri3d(PtVec3d A, PtVec3d B, PtVec3d C, PtVec3d D);  // A constructor
  Quadri3d(PtVec3d P[4]);  // Another constructor
  // Quadri3d(PtVec3d *Pt);  // Yet another constructor
  Quadri3d(const Quadri3d &Q);  // Copy constructor
  ~Quadri3d() {;}  // Destructor
  Quadri3d &operator=(const Quadri3d &Q);  // Copy assignment
  PtVec3d vertex(int i) const {return vrtx[i];}  // 0<=i<4
  PtVec3d center() const;
  double area() const;  // JL20171104: For flat quadrilateral only
  // double diameter() const;
  PtVec3d bilinearmapping(double xhat, double yhat) const;
  double JacobianDeterminant(double xhat, double yhat) const;
  PtVec3d normal(double xhat, double yhat) const;
};

////////////////////////////////////////////////////////////////////////////////

// Tetrahedra

class Tetra {
private:
  PtVec3d vrtx[4];
  double cof[4][4];
public:
  Tetra() {;}  // Default constructor
  Tetra(PtVec3d P0, PtVec3d P1, PtVec3d P2, PtVec3d P3);  // A constructor
  Tetra(PtVec3d *Pt);  // Another constructor
  Tetra(const Tetra &Te);  // Copy constructor
  ~Tetra() {;}  // Destructor
  Tetra &operator=(const Tetra &Te);  // Copy assignment
  void enrich();  // enriching...
  double volume() const;
  double diameter() const;
  PtVec3d vertex(int i) const {return vrtx[i];}  // 0<=i<=3
  PtVec3d center() const {return (1.0/4)*(vrtx[0]+vrtx[1]+vrtx[2]+vrtx[3]);}
  PtVec3d normal(int i) const;  // 0<=i<=3, unit, for face opp. to vertex i
  Tri3d face(int i) const;  // 0<=i<=3, for the face opposite to vertex i
  double coefficient(int i, int j) const {return cof[i][j];}  // 0<=i,j<=3
  void getMaxMinCrds(double &xmax, double &xmin,
                     double &ymax, double &ymin,
                     double &zmax, double &zmin) const;
  void getMaxMinEdges(double &maxEdge, double &minEdge);
  /*
  PtVec3d mapping(double a, double b, double c) const;
  PtVec3d invMapping(PtVec3d P) const;
  int isInTetra(PtVec3d P) const;
  void getBaryCrds(double *baryCrd, PtVec3d P) const;
  Mat3 invJacobian() const;
  */
};

////////////////////////////////////////////////////////////////////////////////

// Hexahedra

class Hexa {
private:
  PtVec3d vrtx[8];
  PtVec3d veca, vecb, vecc, vecd, vece, vecf, vecg;
public:
  Hexa();  // Default constructor
  Hexa(PtVec3d P0, PtVec3d P1, PtVec3d P2, PtVec3d P3,
       PtVec3d P4, PtVec3d P5, PtVec3d P6, PtVec3d P7);  // A constructor
  Hexa(PtVec3d P[8]);  // Yet another constructor
  Hexa(const Hexa &H);  // Copy constructor
  ~Hexa() {;}  // Destructor
  Hexa &operator=(const Hexa&);  // Copy assignment
  void enrich();  // enriching...
  PtVec3d center() const;
  double diameter() const;
  // double volume(const GaussQuad &GQH);  // FOR GENERAL HEXAHEDRA
  // double volume5() const;  // FOR FLAT HEXAHEDRA ONLY
  // double volume6() const;  // FOR FLAT HEXAHEDRA ONLY
  // vector<Tetra> decmp5() const;  // FOR FLAT HEXAHEDRA ONLY
  // vector<Tetra> decmp6() const;  // FOR FLAT HEXAHEDRA ONLY
  // int isInHexa(PtVec3d P);
  PtVec3d vertex(int k) const {return vrtx[k];}  // 0<=k<8
  Quadri3d face(int i) const;  // 0<=i<6
  PtVec3d trilinearmapping(double xhat, double yhat, double zhat) const;
  PtVec3d covariantVector_xhat(double xhat, double yhat, double zhat) const;
  PtVec3d covariantVector_yhat(double xhat, double yhat, double zhat) const;
  PtVec3d covariantVector_zhat(double xhat, double yhat, double zhat) const;
  double JacobianDeterminant(double xhat, double yhat, double zhat) const;
  Mat3 JacobianMatrix(double xhat, double yhat, double zhat) const;
  Mat3 PiolaMatrix(double xhat, double yhat, double zhat) const;
};

////////////////////////////////////////////////////////////////////////////////

// Pyramids
/*
class Pyramid {
private:
  PtVec3d vrtx[5];
public:
  Pyramid();  // Default constructor
  Pyramid(PtVec3d P0, PtVec3d P1, PtVec3d P2, PtVec3d P3, PtVec3d P4);
  Pyramid(PtVec3d *Pt);  // Another constructor
  Pyramid(const Pyramid &Py);  // Copy constructor
  ~Pyramid() {;}  // Destructor
  Pyramid &operator=(const Pyramid&);  // Copy assignment
  double volume();
  void getMaxMinCrds(double &xmax, double &xmin,
  double &ymax, double &ymin, double &zmax, double &zmin);
  void getMaxMinEdges(double &maxEdge, double &minEdge);
  vector<Tetra> decmp2();
  int isInPyramid(PtVec3d P);
};
*/


// Prism
/*
class Prism {
private:
  PtVec3d vrtx[6];
public:
  Prism();  // Default constructor
  Prism(PtVec3d P0, PtVec3d P1, PtVec3d P2,
  PtVec3d P3, PtVec3d P4, PtVec3d P5);  // A constructor
  Prism(PtVec3d *Pt);  // Another constructor
  Prism(const Prism &Pr);  // Copy constructor
  ~Prism() {;}  // Destructor
  Prism &operator=(const Prism&);  // Copy assignment
  double volume();
  void getMaxMinCrds(double &xmax, double &xmin,
  double &ymax, double &ymin, double &zmax, double &zmin);
  void getMaxMinEdges(double &maxEdge, double &minEdge);
  vector<Tetra> decmp3();
  int isInPrism(PtVec3d P);
};
*/

#endif  // Cell3d_H

// cell3d.h
