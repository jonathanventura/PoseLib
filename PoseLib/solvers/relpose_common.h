#pragma once

#include <Eigen/Dense>
#include <vector>

namespace poselib {

int relpose(const Eigen::Matrix<double, 9, Eigen::Dynamic> &epipolar_constraints,
            std::vector<Eigen::Matrix3d> *essential_matrices);

}; // namespace poselib
