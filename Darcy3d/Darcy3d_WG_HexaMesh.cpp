// Darcy3d_WG_HexaMesh.cpp
// James Liu, Graham Harper, ColoState; 2014/07--2017/11

// Solving Darcy in 3d on a hexahedral mesh:
// -- Using WG(P0,P0;AT0) elements:     Good for convex hexahedron with flat faces
// -- Using WG(Q0,Q0;RT[0]) elements:   Good if not too much distortion
// -- Using WG(Q0,Q1;P1^3) elements:    Unsuccessful
// -- Using WG(P1,P0;P0^3)+stabilizer:  Theoretically good


////////////////////////////////////////////////////////////////////////////////
// For solving Darcy eqn. in 3d by the lowest order WG(P0,P0;AT0)Hexa elements
////////////////////////////////////////////////////////////////////////////////

// The single-matrix approach: 4 functions with suffix "1"
/*
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.cpp"

// The Schur-complemnt approach: 4 functions
#include "Darcy3d_WG_HexaP0P0AT0_AsmBndryConds.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_AsmGlbMats.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_AsmSource.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_ModiLinSys.cpp"
// Common functions
#include "Darcy3d_WG_HexaP0P0AT0_Err.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_PresVelFlux.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_ProjPres.cpp"
#include "Darcy3d_WG_HexaP0P0AT0_VeriNmlFluxCnty.cpp"
 */
 
////////////////////////////////////////////////////////////////////////////////
// For solving Darcy eqn. in 3d by the lowest order WG(Q0,Q0;RT[0])Hexa elements
////////////////////////////////////////////////////////////////////////////////

// The single-matrix approach: 4 functions with suffix "1"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.cpp"
// The Schur-complemnt approach: 4 functions
//#include "Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds.cpp"
//#include "Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMats.cpp"
//#include "Darcy3d_WG_HexaQ0Q0RT0_AsmSource.cpp"
//#include "Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys.cpp"
// Common functions 
#include "Darcy3d_WG_HexaQ0Q0RT0_Err.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_ProjPres.cpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_VeriNmlFluxCnty.cpp"

////////////////////////////////////////////////////////////////////////////////
// For solving Darcy in 3d by WG(P1,P0;P0^3)Hexa elements with a stabilizer
////////////////////////////////////////////////////////////////////////////////
/*
// The single-matrix approach: 4 functions with suffix "1"
#include "Darcy3d_WG_HexaP1P0P03_AsmBndryConds1.cpp"
#include "Darcy3d_WG_HexaP1P0P03_AsmGlbMat1.cpp"
#include "Darcy3d_WG_HexaP1P0P03_AsmSource1.cpp"
#include "Darcy3d_WG_HexaP1P0P03_ModiLinSys1.cpp"
// The Schur-complemnt approach: 4 functions
#include "Darcy3d_WG_HexaP1P0P03_AsmBndryConds.cpp"
#include "Darcy3d_WG_HexaP1P0P03_AsmGlbMats.cpp"
#include "Darcy3d_WG_HexaP1P0P03_AsmSource.cpp"
#include "Darcy3d_WG_HexaP1P0P03_ModiLinSys.cpp"
// Common functions
*/
// Darcy3d_WG_HexaMesh.cpp
