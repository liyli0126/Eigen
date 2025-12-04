#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <Eigen/Core>
#include <fstream>  // For file operations

#include "GaussQuad.h"
#include "Hdiv3d.h"
//#include "Darcy3d_WG_HexaMesh.hpp"
#include "HexaMesh.h"
//#include "Write_VTK_File.hpp"

#include "EqnBC_Poisson3d_Ex12coscoscos.h"

#include "Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1.hpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1.hpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_AsmSource1.hpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1.hpp"
#include "Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux.hpp"
#include "Darcy3d_SmplnPerm_HexaMesh.h"


int main(int argc, const char * argv[])
{
    //char filename1a[64] = "Darcy3d_WG_NumerPres.vtk";
    //char filename2a[64] = "Darcy3d_WG_ProjPres.vtk";
   // char varname1a[64] = "NumerPres";
   // char varname2a[64] = "ProjPres";

    clock_t TimeStart = clock();

    GaussQuad GQH, GQQ;
    GQH.setForBrick(3, 3, 3);
    GQQ.setForRectangle(3, 3);

    // Setting up a 3D brick domain: The unit cube
    double xa = 0, xb = 1;
    double yc = 0, yd = 1;
    double ze = 0, zf = 1;

    // Generating a logically rectangular hexahedral mesh
    int n = 16;
    HexaMesh mesh(xa, xb, yc, yd, ze, zf, n, n, n, 1.0);
    mesh.enrich();

    std::cout << "#elements= " << mesh.numberElements() << "  ";
    std::cout << "#faces= " << mesh.numberFaces() << "\n";

    // Sampling permeability
    std::cout << "Sampling permeability...\n";
    Mat3 *PermK = new Mat3[mesh.numberElements()];
    Darcy3d_SmplnPerm_HexaMesh(PermK, fxnMatK, mesh, GQH);
    std::cout << "\n";

    // Assembling and solving
    std::cout << "Assembling and solving...\n";
    int DOFs = mesh.numberElements() + mesh.numberFaces();
    SparseMatrix GlbMat;
    Vector GlbVecSource, GlbVecDirichlet, GlbVecNeumann, GlbRHS, sln(DOFs);
    
    Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1(GlbMat, mesh, PermK, GQH, GQQ);
    Darcy3d_WG_HexaQ0Q0RT0_AsmSource1(GlbVecSource, fxnf, mesh, GQH);
    Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1(GlbVecDirichlet, GlbVecNeumann, fxnpD, fxnuN, mesh, GQQ);
    Darcy3d_WG_HexaQ0Q0RT0_ModiLinSys1(GlbMat, GlbRHS, GlbVecSource, GlbVecDirichlet, GlbVecNeumann, mesh);
    
    int itr = 0;
    slvSpaSpdSysCG(sln, GlbMat, GlbRHS, itr, 5000, 1e-18, 1e-18);
    std::cout << "itr=" << itr << "\n";
    std::cout << "\n";

    // Postprocessing
    Vector NumerPresEm(mesh.numberElements());
    FullMatrix NumerVelCofRT0(mesh.numberElements(), 6);
    FullMatrix NumerFlux(mesh.numberElements(), 6);
    
    Darcy3d_WG_HexaQ0Q0RT0_PresVelFlux(NumerPresEm, NumerVelCofRT0, NumerFlux, mesh, PermK, sln, GQH, GQQ);
    
    
   // Vector ProjPresEm(mesh.numberElements());
    //Vector ProjPresFc(mesh.numberFaces());
   // Darcy3d_WG_HexaQ0Q0RT0_ProjPres(ProjPresEm, ProjPresFc, fxnp, fxnpD, mesh, GQH, GQQ);
    
    
    
    std::cout << "NumerPresEm max = " << NumerPresEm.l0norm() << "\n";
//std::cout << "ProjPresEm  max = " << ProjPresEm.l0norm() << "\n";

    // Output & saving in VTK format
   mesh.save2file("Numerical_Solution.vtk", NumerPresEm, NumerVelCofRT0);


    // Freeing resources
    delete[] PermK;

    // Finishing
    clock_t TimeEnd = clock();
    std::cout << "Time taken = "
              << static_cast<double>(TimeEnd - TimeStart) / CLOCKS_PER_SEC << " seconds \n";

    return 0;  // If successful
}





/*

 // Setup Gauss Quadrature
 void setupGaussQuadrature(GaussQuad &GQH, GaussQuad &GQQ) {
     GQH.setForBrick(3, 3, 3);
     GQQ.setForRectangle(3, 3);
 }

 // Create Hexahedral Mesh
 HexaMesh createMesh(double xa, double xb, double yc, double yd, double ze, double zf, int n, double delta) {
     HexaMesh mesh(xa, xb, yc, yd, ze, zf, n, n, n, delta);
     mesh.enrich();
     return mesh;
 }

 // Print Mesh Information
 void printMeshInfo(const HexaMesh &mesh) {
     std::cout << "#elements= " << mesh.numberElements() << "  ";
     std::cout << "#faces= " << mesh.numberFaces() << "\n";
 }

// Assemble Global Matrix
int assembleGlobalMatrix(Eigen::SparseMatrix<double> &GlbMat, Eigen::VectorXd &GlbVecSource,
                         Eigen::VectorXd &GlbVecDirichlet, Eigen::VectorXd &GlbVecNeumann,
                         const HexaMesh &mesh, const GaussQuad &GQH, const GaussQuad &GQQ,
                         Eigen::Matrix3d *PermK) {
    
    int AsmBndry = Darcy3d_WG_HexaQ0Q0RT0_AsmBndryConds1(GlbVecDirichlet, GlbVecNeumann, fxnpD, fxnuN, mesh, GQQ);
    if (AsmBndry != 0) {
        std::cerr << "Error assembling boundary conditions: " << AsmBndry << std::endl;
        return AsmBndry;
    }
    
    int asmGlbMat = Darcy3d_WG_HexaQ0Q0RT0_AsmGlbMat1(GlbMat, mesh, PermK, GQH, GQQ);
    if (asmGlbMat != 0) {
        std::cerr << "Error assembling global matrix: " << asmGlbMat << std::endl;
        return asmGlbMat;
    }

    Darcy3d_WG_HexaQ0Q0RT0_AsmSource1(GlbVecSource, fxnf, mesh, GQH);
    return 0; // Success
}

// Main Function
int main(int argc, const char * argv[]) {
    GaussQuad GQH, GQQ;
    setupGaussQuadrature(GQH, GQQ);

    double xa = 0, xb = 1, yc = 0, yd = 1, ze = 0, zf = 1;
    int n = 4;
    double delta = 1;
    HexaMesh mesh = createMesh(xa, xb, yc, yd, ze, zf, n, delta);
    printMeshInfo(mesh);
    
    // After creating the mesh, print the internal counts
    std::cout << "Number of Elements: " << mesh.numberElements() << std::endl;
    std::cout << "Number of Faces: " << mesh.numberFaces() << std::endl;
     
    // Initialize global matrix and vectors
    int DOFs = mesh.numberElements() + mesh.numberFaces();
    Eigen::SparseMatrix<double> GlbMat(DOFs, DOFs);
    Eigen::VectorXd GlbVecSource = Eigen::VectorXd::Zero(DOFs);
    Eigen::VectorXd GlbVecDirichlet = Eigen::VectorXd::Zero(DOFs);
    Eigen::VectorXd GlbVecNeumann = Eigen::VectorXd::Zero(DOFs);
    Eigen::VectorXd GlbRHS = Eigen::VectorXd::Zero(DOFs);
    Eigen::VectorXd sln(DOFs);
    
    // Define PermK as an identity matrix
    int numElements = mesh.numberElements();
    Eigen::Matrix3d *PermK = new Eigen::Matrix3d[numElements];
    for (int i = 0; i < numElements; ++i) {
        PermK[i] = Eigen::Matrix3d::Identity();
    }

    // Assemble the global matrix
    int assemblyResult = assembleGlobalMatrix(GlbMat, GlbVecSource, GlbVecDirichlet, GlbVecNeumann, mesh, GQH, GQQ, PermK);
    if (assemblyResult != 0) {
        delete[] PermK; // Clean up before returning
        return assemblyResult;
    }

    // Solve the system
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(GlbMat);
    sln = solver.solve(GlbRHS);
    std::cout << "Solution obtained.\n";

    // Postprocessing
    Eigen::VectorXd NumerPresEm(mesh.numberElements());
    Eigen::MatrixXd NumerVelCofRT0(mesh.numberElements(), 6);
    Eigen::MatrixXd NumerFlux(mesh.numberElements(), 6);
    
    // Compute pressure, velocity, and flux
    int pres = Darcy3d_WG_HexaP0P0RT0_PresVelFlux(NumerPresEm, NumerVelCofRT0, NumerFlux, mesh, sln, GQH, GQQ);
    if (pres != 0) {
        std::cerr << "Error computing pressure, velocity, and flux: " << pres << std::endl;
        return pres;
    }
 
    // Output the pressure results
    std::cout << "Pressure Results:\n";
    for (int i = 0; i < NumerPresEm.size(); ++i) {
         std::cout << "Element " << i << ": Pressure = " << NumerPresEm[i] << "\n";
    }

    // Output & saving in VTK format
    // Write_VTK_File(mesh, NumerPresEm, filename1a, varname1a);

    delete[] PermK; // Clean up

    return 0;  // If successful
}
*/
