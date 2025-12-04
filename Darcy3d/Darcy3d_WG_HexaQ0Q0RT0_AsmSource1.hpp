//
//  Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.hpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#ifndef Darcy3d_WG_HexaQ0Q0RT0_AsmSource1_hpp
#define Darcy3d_WG_HexaQ0Q0RT0_AsmSource1_hpp

#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "HexaMesh.h"
#include "PtVec3d.h"


int Darcy3d_WG_HexaQ0Q0RT0_AsmSource1(Vector &GlbVecSource,
                                      double (*fxnf)(PtVec3d),
                                      const HexaMesh &mesh,
                                      const GaussQuad &GQH);

#endif /* Darcy3d_WG_HexaQ0Q0RT0_AsmSource1_hpp */
