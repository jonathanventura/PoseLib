#include "relpose_5pt.h"
#include "relpose_common.h"

#include "p1ac.h"

#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include <Eigen/Dense>

namespace poselib {

int relpose_1aff(const Eigen::Vector3d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &n, const Eigen::Matrix2d &A,
                std::vector<CameraPose> *output) {
    int n_sols = p1ac(v, u, n, A, output);

    for ( int i = 0; i < n_sols; i++ )
    {
        (*output)[i].t /= (*output)[i].t.norm();
    }

    return n_sols;
}

int relpose_1aff(const Eigen::Vector3d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &n, const Eigen::Matrix2d &A,
                std::vector<Eigen::Matrix3d> *essential_matrices) {

    std::vector<CameraPose> output;

    int n_sols = relpose_1aff(u,v,n,A,&output);

    essential_matrices->clear();
    essential_matrices->reserve(n_sols);
    for ( int i = 0; i < n_sols; i++ )
    {
        Eigen::Matrix3d E;
        essential_from_motion(output[i], &E);
        essential_matrices->emplace_back(E);
    }

    return n_sols;
}

} // namespace poselib