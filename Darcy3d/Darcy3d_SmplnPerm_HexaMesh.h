// Darcy3d_SmplnPerm_HexaMesh.h
// Darcy3d: Sampling permeability on a hexahedral mesh by a Gaussian quadrature
// James Liu, Graham Harper, ColoState; 2014/07--2017/02

#include "GaussQuad.h"
#include "HexaMesh.h"
#include "mat3.h"
#include "PtVec3d.h"

int Darcy3d_SmplnPerm_HexaMesh(Mat3 *PermK, Mat3 (*fxnMatK)(PtVec3d),
                               const HexaMesh &mesh, const GaussQuad &GQH);

// Darcy3d_SmplnPerm_HexaMesh.h