#pragma once

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Computes the surface normal from one affine correspondence.
int normal_1aff(const Eigen::Vector3d &x, const Eigen::Vector3d &X, const Eigen::Matrix2d &A, const Eigen::Matrix3d &R, const Eigen::Vector3d &t, Eigen::Vector3d *output );

}; // namespace poselib