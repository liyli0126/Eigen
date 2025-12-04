//
//  write_VTK.hpp
//  Possion_3D_test1
//
//  Created by Yingli Li on 7/15/24.
//

#ifndef WRITE_VTK_HPP
#define WRITE_VTK_HPP

#include <vector>
#include <Eigen/Dense>
#include <string>

void writeVTK(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const Eigen::VectorXd& V, const std::string& filepath);

#endif // WRITE_VTK_HPP

