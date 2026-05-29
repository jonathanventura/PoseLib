#include "relpose_5pt.h"
#include "relpose_common.h"

#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include <Eigen/Dense>

namespace poselib {

int relpose_5pt(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                std::vector<Eigen::Matrix3d> *essential_matrices) {

    // Compute epipolar constraints
    Eigen::Matrix<double, 9, 5> epipolar_constraints;
    for (int i = 0; i < 5; ++i) {
        epipolar_constraints.col(i) << x1[i](0) * x2[i], x1[i](1) * x2[i], x1[i](2) * x2[i];
    }
    return relpose(epipolar_constraints, essential_matrices);
}

int relpose_5pt(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                std::vector<CameraPose> *output) {
    std::vector<Eigen::Matrix3d> essential_matrices;
    int n_sols = relpose_5pt(x1, x2, &essential_matrices);

    output->clear();
    output->reserve(n_sols);
    for (int i = 0; i < n_sols; ++i) {
        motion_from_essential(essential_matrices[i], x1, x2, output);
    }

    return output->size();
}

} // namespace poselib