#include "write_VTK.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
void writeVTK(const std::vector<double>& x,
              const std::vector<double>& y,
              const std::vector<double>& z,
              const Eigen::VectorXd& V,
              const std::string& filename) {

    // Open the output VTK file
    std::ofstream vtkfile("/Users/yinglili/Desktop/Numerical_Solution.vtk");
    if (!vtkfile) {
        std::cerr << "Error opening VTK file for writing: " << filename << std::endl;
        return;
    }

    // Write the header
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << "VTK output" << std::endl;
    vtkfile << "ASCII" << std::endl;

    // Write grid points (for visualization points)
    vtkfile << "DATASET STRUCTURED_POINTS" << std::endl;
    vtkfile << "DIMENSIONS " << x.size() << " " << y.size() << " " << z.size() << std::endl;
    vtkfile << "ORIGIN " << x[0] << " " << y[0] << " " << z[0] << std::endl;
    vtkfile << "SPACING " << x[1] - x[0] << " " << y[1] - y[0] << " " << z[1] - z[0] << std::endl;

    // Write grid data (for visualization points)
    vtkfile << "POINT_DATA " << x.size() * y.size() * z.size() << std::endl;
    vtkfile << "SCALARS NumericalSolutionPoints double 1" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (std::size_t k = 0; k < z.size(); ++k) {
        for (std::size_t j = 0; j < y.size(); ++j) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                std::size_t index = (k * y.size() + j) * x.size() + i;
                vtkfile << V[index] << std::endl;
            }
        }
    }

    // Write grid data (for visualization cells)
    vtkfile << "CELL_DATA " << static_cast<int>((x.size() - 1) * (y.size() - 1) * (z.size() - 1)) << std::endl;
    vtkfile << "SCALARS NumericalSolutionCells double 1" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (std::size_t k = 0; k < z.size() - 1; ++k) {
        for (std::size_t j = 0; j < y.size() - 1; ++j) {
            for (std::size_t i = 0; i < x.size() - 1; ++i) {
                // Calculate indices for the corners of the cell
                std::size_t index1 = (k * y.size() + j) * x.size() + i;
                std::size_t index2 = index1 + 1;
                std::size_t index3 = index1 + x.size();
                std::size_t index4 = index3 + 1;

                // Output the cell with its value
                vtkfile << V[index1] << std::endl;
                vtkfile << V[index2] << std::endl;
                vtkfile << V[index3] << std::endl;
                vtkfile << V[index4] << std::endl;
            }
        }
    }


       // Write point coordinates (for visualization points)
       vtkfile << "POINTS " << x.size() * y.size() * z.size() << " double" << std::endl;
       for (int k = 0; k < z.size(); ++k) {
           for (int j = 0; j < y.size(); ++j) {
               for (int i = 0; i < x.size(); ++i) {
                   vtkfile << x[i] << " " << y[j] << " " << z[k] << std::endl;
               }
           }
       }
/*
    // Write grid data (for visualization cells)
    vtkfile << "CELLS " << (x.size() - 1) * (y.size() - 1) * (z.size() - 1) << " " << 5 * (x.size() - 1) * (y.size() - 1) * (z.size() - 1) << std::endl;
    for (std::size_t k = 0; k < z.size() - 1; ++k) {
        for (std::size_t j = 0; j < y.size() - 1; ++j) {
            for (std::size_t i = 0; i < x.size() - 1; ++i) {
                // Calculate indices for the corners of the cell
                std::size_t index1 = (k * y.size() + j) * x.size() + i;
                std::size_t index2 = index1 + 1;
                std::size_t index3 = index1 + x.size();
                std::size_t index4 = index3 + 1;

                // Output the cell vertices
                vtkfile << "4 " << index1 << " " << index2 << " " << index4 << " " << index3 << std::endl;
            }
        }
    }

       // Write cell types (for visualization cells)
       vtkfile << "CELL_TYPES " << (x.size() - 1) * (y.size() - 1) * (z.size() - 1) << std::endl;
       for (int k = 0; k < z.size() - 1; ++k) {
           for (int j = 0; j < y.size() - 1; ++j) {
               for (int i = 0; i < x.size() - 1; ++i) {
                   vtkfile << "9" << std::endl;  // VTK_HEXAHEDRON cell type
               }
           }
       }
 */
    // Close the file
    vtkfile.close();

    std::cout << "VTK file " << filename << " written successfully." << std::endl;
}

/*
void writeVTK(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const Eigen::VectorXd& V, const std::string& filepath) {
    std::ofstream file("/Users/yinglili/Desktop/Numerical_Solution.vtk");
    //std::ofstream vtkFile("/Users/yinglili/Desktop/Exact_Solution.vtk");
    
    
    auto size = x.size();
    
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Exact Solution" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << size << " " << size << " " << size << std::endl;

    // Write points
    file << "POINTS " << size * size * size << " float" << std::endl;
    for (double z_val : z) {
        for (double y_val : y) {
            for (double x_val : x) {
                file << x_val << " " << y_val << " " << z_val << std::endl;
            }
        }
    }

    // Write cell data
    file << "CELL_DATA " << (size - 1) * (size - 1) * (size - 1) << std::endl;
    file << "SCALARS numerical_solution float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < size - 1; ++k) {
        for (int j = 0; j < size - 1; ++j) {
            for (int i = 0; i < size - 1; ++i) {
                std::size_t index = (k * (size-1) + j) * (size-1) + i;
                file << V[index];
            }
        }
    }

    file.close();
}
*/

/*
void writeVTK(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const Eigen::VectorXd& V, const std::string& filepath) {
    std::ofstream file("/Users/yinglili/Desktop/Numerical_Solution.vtk");

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        return;
    }

    std::size_t size_x = x.size();
    std::size_t size_y = y.size();
    std::size_t size_z = z.size();
    std::size_t num_cells_x = size_x - 1;
    std::size_t num_cells_y = size_y - 1;
    std::size_t num_cells_z = size_z - 1;
    std::size_t num_cells = num_cells_x * num_cells_y * num_cells_z;

    if (V.size() != size_x * size_y * size_z) {
        std::cerr << "Error: Size of V (" << V.size() << ") does not match the number of points ("
                  << size_x * size_y * size_z << ")." << std::endl;
        return;
    }

    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Numerical Solution" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << num_cells_x << " " << num_cells_y << " " << num_cells_z << std::endl;

    // Write points
    file << "POINTS " << size_x * size_y * size_z << " float" << std::endl;
    for (std::size_t k = 0; k < size_z; ++k) {
        for (std::size_t j = 0; j < size_y; ++j) {
            for (std::size_t i = 0; i < size_x; ++i) {
                file << x[i] << " " << y[j] << " " << z[k] << std::endl;
            }
        }
    }

    // Write cell data
    file << "CELL_DATA " << num_cells << std::endl;
    file << "SCALARS numerical_solution float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (std::size_t k = 0; k < num_cells_z; ++k) {
        for (std::size_t j = 0; j < num_cells_y; ++j) {
            for (std::size_t i = 0; i < num_cells_x; ++i) {
                std::size_t index = (k * num_cells_y + j) * num_cells_x + i;
                file << V[index] << std::endl;
            }
        }
    }

    file.close();
}
*/
