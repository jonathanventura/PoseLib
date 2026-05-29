#include "relpose_5pt.h"
#include "relpose_common.h"

#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include <Eigen/Dense>

namespace poselib {

int relpose_2aff(const std::vector<Eigen::Vector3d> &u, const std::vector<Eigen::Vector3d> &v, const std::vector<Eigen::Matrix2d> &A,
                std::vector<Eigen::Matrix3d> *essential_matrices) {

    // Compute epipolar constraints
    Eigen::Matrix<double, 9, 6> epipolar_constraints;
    for (int i = 0; i < 2; ++i) {
        const double x1 = u[i][0]/u[i][2];
        const double x2 = u[i][1]/u[i][2];
        const double y1 = v[i][0]/v[i][2];
        const double y2 = v[i][1]/v[i][2];
        const double a1 = A[i](0,0);
        const double a2 = A[i](1,0);
        const double a3 = A[i](0,1);
        const double a4 = A[i](1,1);
        epipolar_constraints.col(3*i+0) << x1 * y1, x1 * y2, x1, x2 * y1, x2 * y2, x2, y1, y2, 1;
        epipolar_constraints.col(3*i+1) << a3 * x1, a4 * x1, 0, y1 + a3 * x2, y2 + a4 * x2, 1, a3, a4, 0;
        epipolar_constraints.col(3*i+2) << y1 + a1 * x1, y2 + a2 * x1, 1, a1 * x2, a2 * x2, 0, a1, a2, 0;
    }
    return relpose(epipolar_constraints, essential_matrices);
}

int relpose_2aff(const std::vector<Eigen::Vector3d> &u, const std::vector<Eigen::Vector3d> &v, const std::vector<Eigen::Matrix2d> &A,
                std::vector<CameraPose> *output) {
    std::vector<Eigen::Matrix3d> essential_matrices;
    int n_sols = relpose_2aff(u, v, A, &essential_matrices);

    output->clear();
    output->reserve(n_sols);
    for (int i = 0; i < n_sols; ++i) {
        motion_from_essential(essential_matrices[i], u, v, output);
    }

    return output->size();
}

} // namespace poselib