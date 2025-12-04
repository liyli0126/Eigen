//
//  Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux.hpp
//  Darcy+Eigen
//
//  Created by Yingli Li on 9/26/24.
//

#ifndef Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux_hpp
#define Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux_hpp

#include <cmath>

#include "LinSys.h"
#include "matrix.h"
#include "vector.h"

#include "cell3d.h"
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "HexaMesh.h"
#include "PtVec3d.h"
#include "WG3d_Hexa.h"


int Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(Vector &NumerPresEm,
                                       FullMatrix &NumerVelCofRT0,
                                       FullMatrix &NumerFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Vector &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);


/*
#include "GaussQuad.h"
#include "Hdiv3d.h"
#include "HexaMesh.h"
#include "WG3d_Hexa.h"

#include <Eigen/Dense>

int Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(Eigen::VectorXd &NumerPresEm,
                                       Eigen::MatrixXd &NumerVelCofRT0,
                                       Eigen::MatrixXd &NumerFlux,
                                       const HexaMesh &mesh,
                                       Mat3 *PermK, const Eigen::VectorXd &sln,
                                       const GaussQuad &GQH,
                                       const GaussQuad &GQQ);
*/

#endif /* Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux_hpp */
