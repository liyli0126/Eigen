//
//  Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.cpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#include "Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.hpp"
#include "matrix.h"
#include "vector.h"

#include "HexaMesh.h"

int Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1(SparseMatrix &GlbMat, Vector &GlbRHS,
                                       Vector &GlbVecSource,
                                       Vector &GlbVecDirichlet,
                                       Vector &GlbVecNeumann,
                                       const HexaMesh &mesh)
{
  int labelElementA, labelElementB, labelFaceNeighbor[11];

  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  int DOFs = Num0E + Num1C;
  GlbRHS.resize(DOFs);
  
  // Incorporating the Neumann boundary condition into the global lin.sys.
  GlbRHS = GlbVecSource - GlbVecNeumann;
  // GlbRHS.save2file("GlbRHS.dat");

  // Approach "harsh": Easy
  // Adjusting the glb.disc.lin.sys. according to the Dirichlet bndry.cond.
  /*
  for (labelc=1; labelc<=Num1C; ++labelc) {
    if (mesh.isBoundaryFace(labelc)>0) {
      mesh.getFaceElement(labelc, labelElementA, labelElementB);
      mesh.getFaceFace(labelc,labelFaceNeighbor);
      // For GlbRHS, easy setting ...
      GlbRHS(Num0E+labelc) = DirichletVec(Num0E+labelc);
      // For GlbMat: the very row
      GlbMat.unitfyRow(Num0E+labelc);
    }
  }
  */

  // Approach "gentle": Maintaining symmetry
  // Adjusting the glb.disc.lin.sys. according to the Dirichlet bndry.cond.
  GlbRHS = GlbRHS - GlbMat * GlbVecDirichlet;
  for (int labelc=1; labelc<=Num1C; ++labelc) {
    if (mesh.isBoundaryFace(labelc)>0) {
      mesh.getFaceElement(labelc, labelElementA, labelElementB);
      mesh.getFaceFace(labelc,labelFaceNeighbor);
      GlbRHS(Num0E+labelc) = GlbVecDirichlet(Num0E+labelc);
      GlbMat.unitfyRow(Num0E+labelc);
      GlbMat.setEntry(labelElementA, Num0E+labelc, 0);
      for (int j=0; j<6; ++j) {
        if (labelFaceNeighbor[j]==labelc)  continue;
        GlbMat.setEntry(Num0E+labelFaceNeighbor[j], Num0E+labelc, 0);
      }
    }
  }

  return(0);  // If successful
}
