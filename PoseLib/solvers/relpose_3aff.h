#pragma once

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Computes the fundamental matrix from three affine correspondences.
int relpose_3aff(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2, const std::vector<Eigen::Matrix2d> &A,
                std::vector<Eigen::Matrix3d> *fundamental_matrices);

}; // namespace poselib