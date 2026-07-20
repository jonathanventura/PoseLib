#pragma once

#include <Eigen/Dense>
#include <vector>

namespace poselib {

int relpose_shared_focal(const Eigen::Matrix<double, 9, Eigen::Dynamic> &epipolar_constraints, const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                             ImagePairVector *out_image_pairs);

}; // namespace poselib
