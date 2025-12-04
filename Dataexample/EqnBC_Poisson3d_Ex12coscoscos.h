// EqnBC_Poisson3d_Ex12coscoscos.h
// James Liu, ColoState; 2014/07--2016/05

#ifndef EQNBC_POISSON3D_EX12COSCOSCOS_H
#define EQNBC_POISSON3D_EX12COSCOSCOS_H

#include "mat3.h"
#include "PtVec3d.h"


Mat3 fxnMatK(PtVec3d pt);

double fxnf(PtVec3d pt);

double fxnp(PtVec3d pt);

double fxnpD(PtVec3d pt);

double fxnuN(PtVec3d pt);

PtVec3d fxnpg(PtVec3d pt);

PtVec3d fxnu(PtVec3d pt);

#endif  // EQNBC_POISSON3D_EX12COSCOSCOS_H

// EqnBC_Poisson3d_Ex12coscoscos.h