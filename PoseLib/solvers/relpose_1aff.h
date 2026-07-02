#pragma once

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Computes the essential matrix from one affine correspondence and point normal.
int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n, const Eigen::Matrix2d &A,
                std::vector<Eigen::Matrix3d> *essential_matrices);
int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n, const Eigen::Matrix2d &A,
                std::vector<CameraPose> *output);

}; // namespace poselib