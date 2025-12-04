//
//  Write_VTK_File.hpp
//  Possion_3D_test1
//
//  Created by Yingli Li on 9/17/24.
//
/*
#ifndef Write_VTK_File_hpp
#define Write_VTK_File_hpp

#include "HexaMesh.hpp"
#include <iostream>
#include <fstream>
#include <vector>

// Function to write the mesh and scalar data to a .vtk file
void Write_VTK_File(const HexaMesh& mesh, const Eigen::VectorXd& scalarData, const std::string& filename, const std::string& varname);

#endif /* Write_VTK_File_hpp */

