//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.cpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#include "Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.hpp"
#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "HexaMesh.h"
#include "mat3.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"

int Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1(SparseMatrix &GlbMat,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      const GaussQuad &GQH, const GaussQuad &GQQ)
{
  int i, ii, j, jj, labelc, labele, le;
  int labelElementA, labelElementB;
  int labelFace[6], labelFaceNeighbor[11];
  double EGGM[7][7];
  Hexa hexa;
  Mat3 MatK;
  FullMatrix EltMatK(7,7);

  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  int DOFs = Num0E + Num1C;
  int NumRows = DOFs;
  int NumCols = NumRows;
  int *NumEntrsInRow = new int[NumRows];
  int **ColIdx = new int*[NumRows];

  // Setting up GlbMat: elements: (1)element + (6)faces interaction
  for (ii=0; ii<Num0E; ++ii) {
    labele = mesh.beginLabelElement() + ii;
    mesh.getElementFace(labele, labelFace);
    NumEntrsInRow[ii] = 7;
    ColIdx[ii] = new int[NumEntrsInRow[ii]];
    ColIdx[ii][0] = labele;
    for (j=0; j<6; ++j)  ColIdx[ii][1+j] = Num0E + labelFace[j];
  }

  // Setting up GlbMat: faces: (1 or 2) elements + (6 or 11) faces interaction
  // std::cout << "Displaying boundary faces: ";
  for (ii=Num0E; ii<Num0E+Num1C; ++ii) {
    labelc = (ii-Num0E) + mesh.beginLabelFace();
    mesh.getFaceElement(labelc, labelElementA, labelElementB);
    mesh.getFaceFace(labelc, labelFaceNeighbor);
    if (mesh.isBoundaryFace(labelc)==0) {
      NumEntrsInRow[ii] = 13;  // 13=2+11
      ColIdx[ii] = new int[NumEntrsInRow[ii]];
      ColIdx[ii][0] = labelElementA;
      ColIdx[ii][1] = labelElementB;
      for (j=0; j<11; ++j)
        ColIdx[ii][2+j] = Num0E + labelFaceNeighbor[j];
    }
    if (mesh.isBoundaryFace(labelc)>0) {
      // std::cout << labelc << "  " << std::flush;
      NumEntrsInRow[ii] = 7;  // 7=1+6
      ColIdx[ii] = new int[NumEntrsInRow[ii]];
      ColIdx[ii][0] = labelElementA;
      for (j=0; j<6; ++j)
        ColIdx[ii][1+j] = Num0E + labelFaceNeighbor[j];
    }
  }
  // std::cout << "done!\n";

  // Setting up GlbMat: resize & release
  GlbMat.resize(NumRows,NumCols,NumEntrsInRow,ColIdx);
  for (ii=0; ii<NumRows; ++ii)  delete[] ColIdx[ii];
  delete[] ColIdx;
  delete[] NumEntrsInRow;

  // Assembling the global matrices
  std::cout << "le/1000= ";
  for (labele=1; labele<=mesh.numberElements(); ++labele) {
    // std::cout << "labele=" << labele << "\n";
    le = labele - 1;
    if (le%1000==0)  std::cout << le/1000 << "  " << std::flush;
    mesh.getElementFace(labele,labelFace);
    hexa = mesh.element(labele);
    hexa.enrich();

    // JL20150502: THIS IS THE WAY TO GO
    MatK = PermK[labele-1];
    // std::cout << MatK << "\n";

    EltMatK = WG3d_HexaQ0Q0RT0_EltGradGradMatK(hexa, MatK, GQH, GQQ);
    // std::cout << EltMatK << "\n";
    // std::cout << "\n";

    for (i=0; i<7; ++i)
      for (int j=0; j<7; ++j)
        EGGM[i][j] = EltMatK(i+1,j+1);

    i = 0;  ii = labele;
    j = 0;  jj = labele;
    GlbMat.addEntry(ii, jj, EGGM[i][j]);
    for (j=0; j<6; ++j) {
      jj = Num0E + labelFace[j];
      GlbMat.addEntry(ii, jj, EGGM[i][j+1]);
    }
    for (i=0; i<6; ++i) {
      ii = Num0E + labelFace[i];
      j = 0;  jj = labele;
      GlbMat.addEntry(ii, jj, EGGM[i+1][j]);
      for (j=0; j<6; ++j) {
        jj = Num0E + labelFace[j];
        GlbMat.addEntry(ii, jj, EGGM[i+1][j+1]);
      }
    }
  }
  std::cout << "done!\n";

  return(0);  // If successful
}

// Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.cpp
