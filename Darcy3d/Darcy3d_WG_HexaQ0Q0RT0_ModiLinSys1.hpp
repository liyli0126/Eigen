//
//  Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.hpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#ifndef Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1_hpp
#define Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1_hpp

#include "matrix.h"
#include "vector.h"

#include "HexaMesh.h"

int Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1(SparseMatrix &GlbMat, Vector &GlbRHS,
                                       Vector &GlbVecSource,
                                       Vector &GlbVecDirichlet,
                                       Vector &GlbVecNeumann,
                                       const HexaMesh &mesh);

#endif /* Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1_hpp */
