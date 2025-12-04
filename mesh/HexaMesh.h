// HexaMesh.h
// Class for hexahedral meshes
// James Liu, Graham Harper, ColoState; 2007/01--2017/03

#ifndef HEXAMESH_H
#define HEXAMESH_H

#include "cell3d.h"
#include "PtVec3d.h"
// #include "TetraMesh.h"

class HexaMesh {
private:
  int NumNds, NumFcs, NumEms, NumBndryPcs;
  int BgnLblNd, BgnLblFc, BgnLblEm, BgnLblBndry;  // Default: 1
  int *IsBndryFc, (*LblFcNd)[4], (*LblFcFc)[11], (*LblFcEm)[2];
  int(*LblEmNd)[8], (*LblEmFc)[6];
  PtVec3d *nd, *OutNmlBndry;
  int flag;  // 0,1,2 for private use
  int NX, NY, NZ;  // For private use in enrich() for logically rectangular mesh
public:
  // Constructor #0: Default empty constructor
  HexaMesh();
  //
  // Constructor #1: Need to call enrich() later
  // For a brick domain: Prescribed perturbation with magnitude delta
  HexaMesh(double xa, double xb,
           double yc, double yd,
           double ze, double zf,
           int nx, int ny, int nz,
           double delta);
  //
  // Constructor #2: Need to call enrich() later
  // For a brick domain: Random perturbations in x-,y-,z-directions
  HexaMesh(double xa, double xb,
           double yc, double yd,
           double ze, double zf,
           int nx, int ny, int nz,
           double xdelta, double ydelta, double zdelta);
  //
  // Constructor #3: Need to call enrich() later
  // For a brick domain: Trapezoidal perturbations in y-direction
  HexaMesh(double xa, double xb, int nx,
           double yc, double yd, int ny, double alpha,
           double ze, double zf, int nz);
  //
  // Constructor #4: Need to call enrich() later
  // For a brick domain: Trapezoidal perturbations in x-,y-directions
  HexaMesh(double xa, double xb,
           double yc, double yd,
           double ze, double zf,
           int nx, int ny, int nz,
           double xdelta, double ydelta);
  //
  // Constructor #5: Need to call enrich() later
  // For a cylindrical sector (nonclose)
  HexaMesh(double ri, double ro,
           double alpha, double beta,
           double zb, double zt,
           int nr, int ntheta, int nz);
  //
  // Constructor #6: Need to call enrich() later
  // For a spherical sector (nonclose)
  HexaMesh(double ri, double ro, int nr,
           double alpha, double beta, int ntheta,
           double gamma, double delta, int nphi);
  //
  ~HexaMesh();  // Destructor
  //
  // enriching mesh info for logically rectangular mesh (flag=1)
  void enrich();
  //
  int numberNodes() const {return NumNds;}
  int numberFaces() const {return NumFcs;}
  int numberElements() const {return NumEms;}
  //
  int beginLabelNode() const {return BgnLblNd;}
  int beginLabelFace() const {return BgnLblFc;}
  int beginLabelElement() const {return BgnLblEm;}
  int endLabelNode() const {return BgnLblNd+NumNds-1;}
  int endLabelFace() const {return BgnLblFc+NumFcs-1;}
  int endLabelElement() const {return BgnLblEm+NumEms-1;}
  //
  int isBoundaryFace(int labelc) const {return IsBndryFc[labelc-BgnLblFc];}
  int numberBoundaryFaces() const;
  // PtVec3d outerNormalBoundary(int labelc) const {return OutNmlBndry[labelc-1];}
  //
  PtVec3d node(int labeld) const {return nd[labeld-BgnLblNd];}  // 1<=labeld
  Quadri3d face(int labelc) const;  // 1<=labelc
  Hexa element(int labele) const;  // 1<=labele
  //
  void getFaceNode(int labelc, int labelVertex[4]) const;
  void getFaceFace(int labelc, int labelFaceNeighbor[11]) const;
  void getFaceElement(int labelc, int &labelElementA, int &labelElementB) const;
  void getElementNode(int labele, int labelVertex[8]) const;  // 1<=labele
  void getElementFace(int labele, int labelFace[6]) const;
  //
  double diameter() const;
  //
  // void fillNodeInfo(int numberNodes, double (*crd)[3], int *boundaryNodeMark);
  // void fillFaceInfo(int numberFaces, int (*faceNode)[4],
  //                   int (*faceElement)[2], int *boundaryFaceMark);
  // void fillElementInfo(int numberElements, int (*elementNode)[8]);
  //
    void save2file(const char *filename, const Vector& scalarData, const FullMatrix& velocity) const;
  //
// NOTE: Function "THex" is separately implemented in "THex.cpp"   
// friend int THex(HexaMesh &hexamesh, const TetraMesh &tetramesh);
// NOTE: Function "HexaMeshRegRefi" is separately implemented in "HexaMeshRegRefi.cpp"
// friend int HexaMeshRegRefi(HexaMesh &newmesh, HexaMesh &oldmesh);
};

#endif  // HEXAMESH_H
// HexaMesh.h
