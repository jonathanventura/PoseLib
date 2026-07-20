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

// Computes the essential matrix from one affine correspondence and two point normals.
int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                std::vector<Eigen::Matrix3d> *essential_matrices);
int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                std::vector<CameraPose> *output);

// Computes the essential matrix and shared focal from one affine correspondence and two point normals.
int relpose_1aff_focal(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                ImagePairVector *output);

}; // namespace poselib