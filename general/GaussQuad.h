// GaussQuad.h
// James Liu, Graham Harper, ColoState; 2007/01--2017/02

#ifndef GAUSSQUAD_H
#define GAUSSQUAD_H

class GaussQuad {
private:
  int NumQuadPts;
  int NumVrtx;
  double *absc;  // abscissa
  double *wt;  // weight
  int dim;  // dimension of the reference element this quadrature is set for
  double **crd;  // The Cartesian coordinates of the standard quadrature points
public:
  // Default constructor
  GaussQuad() {NumQuadPts=0; NumVrtx=0; absc=0; wt=0; dim=0; crd=0;}
  ~GaussQuad();  // Destructor
  int numberQuadraturePoints() const {return NumQuadPts;}
  int numberVertices() const {return NumVrtx;}
  double baryCoordinate(int i, int j) const;  // 0<=i<NumQuadPts, 0<=j<NumVrtx
  double weight(int i) const {return wt[i];}  // 0<=i<NumQuadPts
  int dimension() const {return dim;}  // dimension of the reference element
  double CartesianCoordinate(int i, int j) const {return crd[i][j];}  // CAUTION!!!
  void save2file(char *filename) const;  // saving to a data file;
  void setForInterval(int NumQuadPtsX);
  void setForRectangle(int NumQuadPtsX, int NumQuadPtsY);
  void setForTriangle(int NumQuadPtsT);
  void setForBrick(int NumQuadPtsX, int NumQuadPtsY, int NumQuadPtsZ);
  void setForTetrahedron(int NumQuadPtsTe);
};

#endif  // GAUSSQUAD_H