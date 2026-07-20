#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include "relpose_focal_common.h"

#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <stdio.h>

namespace poselib {

int relpose_6pt_shared_focal(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                             ImagePairVector *out_image_pairs) {

    // Compute nullspace to epipolar constraints
    Eigen::Matrix<double, 9, 6> epipolar_constraints;
    for (size_t i = 0; i < 6; ++i) {
        epipolar_constraints.col(i) << x1[i](0) * x2[i], x1[i](1) * x2[i], x1[i](2) * x2[i];
    }
    return relpose_shared_focal(epipolar_constraints, x1, x2, out_image_pairs);
}
} // namespace poselib