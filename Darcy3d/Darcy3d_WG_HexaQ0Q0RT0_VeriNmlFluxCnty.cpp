// Darcy3d_WG_HexaQ0Q0RT0_VeriNmlFluxCnty.cpp
// James Liu, Graham Harper, ColoState; 2014/07--2016/12

#include <cmath>

#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"


int Darcy3d_WG_HexaQ0Q0RT0_VeriNmlFluxCnty(const HexaMesh &mesh,
                                           const FullMatrix &NmlFlux)
{
  int labelFace[6];
  
  // Setup
  int Num0E = mesh.numberElements();
  int Num1C = mesh.numberFaces();
  Vector FluxDscp(Num1C);
  
  // Checking
  for (int le=0; le<Num0E; ++le) {
    int labele = le + mesh.beginLabelElement();
    mesh.getElementFace(labele, labelFace);
    for (int j=0; j<6; ++j) {
      FluxDscp(labelFace[j]) += NmlFlux(labele,j+1);
    }
  }

  // Modification for boundary faces
  for (int lc=0; lc<Num1C; ++lc) {
    int labelc = lc + mesh.beginLabelFace();
    if (mesh.isBoundaryFace(labelc)>0)  FluxDscp(labelc) = 0;
  }

  std::cout << "Max. flux discrepancy= " << FluxDscp.l0norm() << "\n";

  return 0;  // If successful
}
