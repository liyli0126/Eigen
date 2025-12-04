//
//  Write_VTK_File.cpp
//  Possion_3D_test1
//
//  Created by Yingli Li on 9/17/24.
//


/*
#include "Write_VTK_File.hpp"
#include "HexaMesh.hpp"
#include <iostream>
#include <fstream>
#include <vector>

// Function to write the mesh and scalar data to a .vtk file
void Write_VTK_File(const HexaMesh& mesh, const Eigen::VectorXd& scalarData, const std::string& filename, const std::string& varname) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    int nx = mesh.getNX();
    int ny = mesh.getNY();
    int nz = mesh.getNZ();  

    int numNodes = (nx + 1) * (ny + 1) * (nz + 1);
    int numCells = nx * ny * nz;

    file << "# vtk DataFile Version 3.0\n";
    file << "Hexahedral mesh\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << (nx + 1) << " " << (ny + 1) << " " << (nz + 1) << "\n";
    file << "POINTS " << numNodes << " float\n";

    // Write points
    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                Eigen::Vector3d node = mesh.node(i + j * (nx + 1) + k * (nx + 1) * (ny + 1));
                file << node.x() << " " << node.y() << " " << node.z() << "\n";
            }
        }
    }

    file << "CELLS " << numCells << " " << (numCells * 9) << "\n";

    // Write cells
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cellId = i + nx * (j + ny * k);
                file << "8 "
                     << (i + (j + k * (ny + 1))) << " "
                     << (i + 1 + (j + k * (ny + 1))) << " "
                     << (i + 1 + ((j + 1) + k * (ny + 1))) << " "
                     << (i + ((j + 1) + k * (ny + 1))) << " "
                     << (i + ((j) + (k + 1) * (ny + 1))) << " "
                     << (i + 1 + ((j) + (k + 1) * (ny + 1))) << " "
                     << (i + 1 + ((j + 1) + (k + 1) * (ny + 1))) << " "
                     << (i + ((j + 1) + (k + 1) * (ny + 1))) << "\n";
            }
        }
    }

    file << "CELL_TYPES " << numCells << "\n";

    // Write cell types
    for (int i = 0; i < numCells; ++i) {
        file << "12\n";  // VTK_HEXAHEDRON
    }

    file << "POINT_DATA " << numNodes << "\n";
    file << "SCALARS " << varname << " float 1\n";
    file << "LOOKUP_TABLE default\n";

    // Write scalar data
    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                int nodeId = i + (j + k * (ny + 1));
                file << scalarData[nodeId] << "\n";
            }
        }
    }

    file.close();
}
*/
