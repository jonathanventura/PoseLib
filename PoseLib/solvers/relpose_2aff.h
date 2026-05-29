#pragma once

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Computes the essential matrix from two affine correspondences.
//    Raposo and Barreto, Theory and Practice fo Structure-from-Motion using Affine Correspondences, CVPR 2016
int relpose_2aff(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2, const std::vector<Eigen::Matrix2d> &A,
                std::vector<Eigen::Matrix3d> *essential_matrices);
int relpose_2aff(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2, const std::vector<Eigen::Matrix2d> &A,
                std::vector<CameraPose> *output);

}; // namespace poselib