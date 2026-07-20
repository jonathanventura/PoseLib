#pragma once

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Solves for relative pose with one unknown focal length from 2 affine correspondences
// The solver is created using Larsson et al. CVPR 2017
int relpose_2aff_shared_focal(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2, const std::vector<Eigen::Matrix2d> &A,
                             ImagePairVector *out_image_pairs);
}; // namespace poselib