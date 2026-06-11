#pragma once

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>

namespace poselib {

int p1ac(const Eigen::Vector3d &x, const Eigen::Vector3d &X, const Eigen::Vector3d &n, const Eigen::Matrix2d &A, std::vector<CameraPose> *output);

} // namespace poselib
