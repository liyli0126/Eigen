//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.hpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#ifndef Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1_hpp
#define Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1_hpp
#include <cmath>

#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"


int Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1(Vector &GlbVecDirichlet,
                                          Vector &GlbVecNeumann,
                                          double (*fxnpD)(PtVec3d),
                                          double (*fxnuN)(PtVec3d),
                                          const HexaMesh &mesh,
                                          const GaussQuad &GQQ);

#endif /* Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1_hpp */
