//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.hpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#ifndef Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1_hpp
#define Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1_hpp

#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "HexaMesh.h"
#include "mat3.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"

int Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1(SparseMatrix &GlbMat,
                                      const HexaMesh &mesh, Mat3 *PermK,
                                      const GaussQuad &GQH, const GaussQuad &GQQ);

#endif /* Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1_hpp */
