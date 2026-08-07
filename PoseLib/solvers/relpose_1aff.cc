#include "relpose_5pt.h"
#include "relpose_common.h"

#include "p1ac.h"

#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include <Eigen/Dense>
#include "Polynomial/Polynomial.hpp"
#include <unsupported/Eigen/Polynomials>

namespace poselib {

int relpose_1aff(const Eigen::Vector3d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &n, const Eigen::Matrix2d &A,
                std::vector<CameraPose> *output) {
    std::vector<CameraPose> p1ac_output;

    int n_sols = p1ac(v, u, n, A, &p1ac_output);

    output->clear();
    for ( int i = 0; i < n_sols; i++ )
    {
        CameraPose pose = p1ac_output[i];
        pose.t.normalize();
        if ( check_cheirality(pose, u, v ) )
        {
            output->emplace_back(pose);
        }
    }

    return output->size();
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


int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                std::vector<CameraPose> *output) {

    const Eigen::Vector2d pref(x1(0)/x1(2),x1(1)/x1(2));
    const Eigen::Vector2d pquery(x2(0)/x2(2),x2(1)/x2(2));
    const Eigen::Vector3d nref = n1.normalized();
    const Eigen::Vector3d r = n2.normalized();

    // determine rotation from n1 to n2
    const Eigen::Matrix3d Rn = Eigen::Quaterniond::FromTwoVectors(n1,n2).toRotationMatrix();

    // build linear constraints
    Eigen::Matrix<double,4,6> M;

    // version using only diagonal of A
    M << 1, 0, -pquery(0), pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) - pquery(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2), 2*Rn(2,2)*r(1) - 2*Rn(1,2)*r(2) + pquery(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) - pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) - pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1)), Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1) - pquery(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 1, -pquery(1), pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) - pquery(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2), 2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pquery(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0)), Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1) - pquery(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 0, A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), nref(0)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) - pquery(0)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) + A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), (2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - pquery(0)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0)))) - A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), nref(0)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - pquery(0)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - Rn(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 0, A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), nref(1)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) - pquery(0)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), (2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - pquery(0)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0)))) - A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), nref(1)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - pquery(0)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - Rn(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1));

    // const double sqrtdetA = sqrt(A.determinant());
    // const double cref = sqrt(2)/2;
    // const double sref = cref;
    // const double cquery = (A(0,0)*cref + A(0,1)*sref)/sqrtdetA;
    // const double squery = (A(1,0)*cref + A(1,1)*sref)/sqrtdetA;
    // version using sref, cref, cquery, squery, sqrtdetA
    // M << 1, 0, -pquery(0), pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) - pquery(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2), 2*Rn(2,2)*r(1) - 2*Rn(1,2)*r(2) + pquery(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) - pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) - pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1)), Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1) - pquery(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 1, -pquery(1), pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) - pquery(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2), 2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pquery(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0)), Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1) - pquery(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 0, -cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), cref*(pquery(0)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - nref(0)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2))) + sref*(pquery(0)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - nref(1)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2))) - cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), cref*(nref(0)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - (2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + pquery(0)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sref*(nref(1)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - (2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + pquery(0)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), cref*(pquery(0)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(0)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) + Rn(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) + sref*(pquery(0)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(1)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) + Rn(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 0, -sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), cref*(pquery(1)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2))) + sref*(pquery(1)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2))) - sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), cref*((2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))) + pquery(1)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sref*((2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))) + pquery(1)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), cref*(pquery(1)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(0)*(Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) + Rn(1,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) + sref*(pquery(1)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(1)*(Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) + Rn(1,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1));
    
    // G-J elimination
    const Eigen::Matrix<double,4,3> C = M.block<4,4>(0,0).partialPivLu().solve(M.block<4,3>(0,3));

    // solve quadratic polynomial
    Eigen::Vector3d poly(C(3,2),C(3,1),C(3,0));
    Eigen::PolynomialSolver<double,2> poly_solver(poly);
    std::vector<double> q_solns;
    poly_solver.realRoots(q_solns);

    int nsolns = q_solns.size();
    for ( int i = 0; i < nsolns; i++ )
    {
        double q = q_solns[i];
        double theta = 2*atan(q);

        Eigen::Matrix3d Rq = Eigen::Quaterniond(Eigen::AngleAxisd(theta,r)).toRotationMatrix();

        // Eigen::Matrix<double,6,4> Mt;
        // Mt << 1, 0, -pquery(0), Rn(2,2)*(2*q*r(1) + 2*q*q*r(0)*r(2)) - pref(1)*(Rn(0,1)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,1)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,1)*(2*q*r(1) + 2*q*q*r(0)*r(2))) - Rn(0,2)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(1,2)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - pref(0)*(Rn(0,0)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,0)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,0)*(2*q*r(1) + 2*q*q*r(0)*r(2))) + pquery(0)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2))),
        // 0, 1, -pquery(1), Rn(0,2)*(2*q*r(2) + 2*q*q*r(0)*r(1)) - pref(1)*(Rn(1,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,1)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,1)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - Rn(1,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - pref(0)*(Rn(1,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,0)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,0)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - Rn(2,2)*(2*q*r(0) - 2*q*q*r(1)*r(2)) + pquery(1)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2))),
        // 0, 0, A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(0,0)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,0)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,0)*(2*q*r(1) + 2*q*q*r(0)*r(2))) - pquery(0)*((nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) - nref(0)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2)))) - nref(0)*(pref(0)*(Rn(0,0)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,0)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,0)*(2*q*r(1) + 2*q*q*r(0)*r(2))) + pref(1)*(Rn(0,1)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,1)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,1)*(2*q*r(1) + 2*q*q*r(0)*r(2))) + Rn(0,2)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,2)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,2)*(2*q*r(1) + 2*q*q*r(0)*r(2))) - A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2))),
        // 0, 0, A(1,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(1,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,0)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,0)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - pquery(1)*((nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) - nref(0)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2)))) - nref(0)*(pref(0)*(Rn(1,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,0)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,0)*(2*q*r(0) - 2*q*q*r(1)*r(2))) + pref(1)*(Rn(1,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,1)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,1)*(2*q*r(0) - 2*q*q*r(1)*r(2))) + Rn(1,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,2)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,2)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - A(1,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2))),
        // 0, 0, A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(0,1)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,1)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,1)*(2*q*r(1) + 2*q*q*r(0)*r(2))) - pquery(0)*((nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) - nref(1)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2)))) - nref(1)*(pref(0)*(Rn(0,0)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,0)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,0)*(2*q*r(1) + 2*q*q*r(0)*r(2))) + pref(1)*(Rn(0,1)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,1)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,1)*(2*q*r(1) + 2*q*q*r(0)*r(2))) + Rn(0,2)*(2*q*q*r(1)*r(1) + 2*q*q*r(2)*r(2) - q*q - 1) + Rn(1,2)*(2*q*r(2) - 2*q*q*r(0)*r(1)) - Rn(2,2)*(2*q*r(1) + 2*q*q*r(0)*r(2))) - A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2))),
        // 0, 0, A(1,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(1,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,1)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,1)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - pquery(1)*((nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) - nref(1)*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2)))) - nref(1)*(pref(0)*(Rn(1,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,0)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,0)*(2*q*r(0) - 2*q*q*r(1)*r(2))) + pref(1)*(Rn(1,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,1)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,1)*(2*q*r(0) - 2*q*q*r(1)*r(2))) + Rn(1,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(2)*r(2) - q*q - 1) - Rn(0,2)*(2*q*r(2) + 2*q*q*r(0)*r(1)) + Rn(2,2)*(2*q*r(0) - 2*q*q*r(1)*r(2))) - A(1,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(Rn(2,0)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,0)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,0)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + pref(1)*(Rn(2,1)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,1)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,1)*(2*q*r(0) + 2*q*q*r(1)*r(2))) + Rn(2,2)*(2*q*q*r(0)*r(0) + 2*q*q*r(1)*r(1) - q*q - 1) + Rn(0,2)*(2*q*r(1) - 2*q*q*r(0)*r(2)) - Rn(1,2)*(2*q*r(0) + 2*q*q*r(1)*r(2)));
        // compute null vector of Mt
        // Eigen::JacobiSVD< Eigen::Matrix<double,6,4> > svdMt(Mt, Eigen::ComputeFullV);
        // Eigen::Matrix<double,4,1> th = svdMt.matrixV().col(3);
        // Eigen::Vector3d t = th.head(3)/th(3);

        Eigen::Vector3d mon( q*q, q, 1 );
        Eigen::Vector3d t = -C.block<3,3>(0,0)*mon;

        Eigen::Matrix3d Rsoln = Rq*Rn;
        Eigen::Vector3d tsoln = t.normalized();

        CameraPose pose(Rsoln,tsoln);

        if ( check_cheirality(pose, x1, x2 ) )
        {
            output->emplace_back(pose);
        }
    }

    return output->size();
}

int relpose_1aff(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                std::vector<Eigen::Matrix3d> *essential_matrices) {

    std::vector<CameraPose> output;

    int n_sols = relpose_1aff(x1,x2,n1,n2,A,&output);

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

int relpose_1aff_focal(const Eigen::Vector3d &x1, const Eigen::Vector3d &x2, const Eigen::Vector3d &n1, const Eigen::Vector3d &n2, const Eigen::Matrix2d &A,
                ImagePairVector *output)
{
    const Eigen::Vector2d pref(x1(0)/x1(2),x1(1)/x1(2));
    const Eigen::Vector2d pquery(x2(0)/x2(2),x2(1)/x2(2));
    const Eigen::Vector3d nref = n1.normalized();
    const Eigen::Vector3d nquery = n2.normalized();

    // determine rotation from n1 to n2
    const Eigen::Matrix3d Rn = Eigen::Quaterniond::FromTwoVectors(n1,n2).toRotationMatrix();

    Eigen::Matrix4d C0, C1, C2;
    C2 << 0, pquery(0)*((nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) - nref(0)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))) + A(0,0)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(0)*(nref(0)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))) - (2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1))) - A(0,0)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)))*(nref(0)*pref(0) + nref(1)*pref(1)), A(0,0)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - pquery(0)*(nref(0)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(0)*pref(0) + nref(1)*pref(1))),
    0, pquery(1)*((nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) - nref(0)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))) + A(1,0)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(1)*(nref(0)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))) - (2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1))) - A(1,0)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)))*(nref(0)*pref(0) + nref(1)*pref(1)), A(1,0)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - pquery(1)*(nref(0)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(0)*pref(0) + nref(1)*pref(1))),
    0, pquery(0)*((nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)) - nref(1)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))) + A(0,1)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(0)*(nref(1)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))) - (2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1))) - A(0,1)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)))*(nref(0)*pref(0) + nref(1)*pref(1)), A(0,1)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - pquery(0)*(nref(1)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(0)*pref(0) + nref(1)*pref(1))),
    0, pquery(1)*((nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)) - nref(1)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))) + A(1,1)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)))*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(1)*(nref(1)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))) - (2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1))) - A(1,1)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)))*(nref(0)*pref(0) + nref(1)*pref(1)), A(1,1)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - pquery(1)*(nref(1)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(0)*pref(0) + nref(1)*pref(1)));
    C1 << A(0,0)*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(0)*(nref(2)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) - nref(0)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2))) - (nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,0)*nquery(0)*nquery(1) - Rn(0,0)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(0)*nquery(2)) + nref(0)*(pref(0)*(2*Rn(1,0)*nquery(0)*nquery(1) - Rn(0,0)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(0)*nquery(2)) + pref(1)*(2*Rn(1,1)*nquery(0)*nquery(1) - Rn(0,1)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(0)*nquery(2))) + A(0,0)*(nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)) + A(0,0)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2))), (2*Rn(1,0)*nquery(2) - 2*Rn(2,0)*nquery(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(pref(0)*(2*Rn(1,0)*nquery(2) - 2*Rn(2,0)*nquery(1)) + pref(1)*(2*Rn(1,1)*nquery(2) - 2*Rn(2,1)*nquery(1))) - pquery(0)*(nref(2)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) - nref(0)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))) - A(0,0)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(0,0)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))), pquery(0)*(Rn(2,0)*nref(2) - Rn(2,2)*nref(0)) + nref(0)*(Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - Rn(0,0)*(nref(0)*pref(0) + nref(1)*pref(1)) + A(0,0)*nref(2)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) + A(0,0)*Rn(2,2)*(nref(0)*pref(0) + nref(1)*pref(1)),
    A(1,0)*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(1)*(nref(2)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) - nref(0)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2))) - (nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*nquery(0)*nquery(1) - Rn(1,0)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(1)*nquery(2)) + nref(0)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(1) - Rn(1,0)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(1) - Rn(1,1)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(1)*nquery(2))) + A(1,0)*(nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)) + A(1,0)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2))), nref(0)*(pref(0)*(2*Rn(0,0)*nquery(2) - 2*Rn(2,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(2) - 2*Rn(2,1)*nquery(0))) - pquery(1)*(nref(2)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) - nref(0)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))) - (2*Rn(0,0)*nquery(2) - 2*Rn(2,0)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(1,0)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(1,0)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))), pquery(1)*(Rn(2,0)*nref(2) - Rn(2,2)*nref(0)) + nref(0)*(Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) - Rn(1,0)*(nref(0)*pref(0) + nref(1)*pref(1)) + A(1,0)*nref(2)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) + A(1,0)*Rn(2,2)*(nref(0)*pref(0) + nref(1)*pref(1)),
    A(0,1)*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(0)*(nref(2)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)) - nref(1)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2))) - (nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,1)*nquery(0)*nquery(1) - Rn(0,1)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(0)*nquery(2)) + nref(1)*(pref(0)*(2*Rn(1,0)*nquery(0)*nquery(1) - Rn(0,0)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(0)*nquery(2)) + pref(1)*(2*Rn(1,1)*nquery(0)*nquery(1) - Rn(0,1)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(0)*nquery(2))) + A(0,1)*(nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)) + A(0,1)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2))), (2*Rn(1,1)*nquery(2) - 2*Rn(2,1)*nquery(1))*(nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(pref(0)*(2*Rn(1,0)*nquery(2) - 2*Rn(2,0)*nquery(1)) + pref(1)*(2*Rn(1,1)*nquery(2) - 2*Rn(2,1)*nquery(1))) - pquery(0)*(nref(2)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)) - nref(1)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))) - A(0,1)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(0,1)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))), pquery(0)*(Rn(2,1)*nref(2) - Rn(2,2)*nref(1)) + nref(1)*(Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - Rn(0,1)*(nref(0)*pref(0) + nref(1)*pref(1)) + A(0,1)*nref(2)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) + A(0,1)*Rn(2,2)*(nref(0)*pref(0) + nref(1)*pref(1)),
    A(1,1)*(nref(0)*pref(0) + nref(1)*pref(1)), pquery(1)*(nref(2)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2)) - nref(1)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2))) - (nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*nquery(0)*nquery(1) - Rn(1,1)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(1)*nquery(2)) + nref(1)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(1) - Rn(1,0)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(1) - Rn(1,1)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(1)*nquery(2))) + A(1,1)*(nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)) + A(1,1)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(0)*nquery(2) - Rn(2,0)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,0)*nquery(1)*nquery(2)) + pref(1)*(2*Rn(0,1)*nquery(0)*nquery(2) - Rn(2,1)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,1)*nquery(1)*nquery(2))), nref(1)*(pref(0)*(2*Rn(0,0)*nquery(2) - 2*Rn(2,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(2) - 2*Rn(2,1)*nquery(0))) - pquery(1)*(nref(2)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0)) - nref(1)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))) - (2*Rn(0,1)*nquery(2) - 2*Rn(2,1)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(1,1)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0))*(nref(0)*pref(0) + nref(1)*pref(1)) - A(1,1)*nref(2)*(pref(0)*(2*Rn(0,0)*nquery(1) - 2*Rn(1,0)*nquery(0)) + pref(1)*(2*Rn(0,1)*nquery(1) - 2*Rn(1,1)*nquery(0))), pquery(1)*(Rn(2,1)*nref(2) - Rn(2,2)*nref(1)) + nref(1)*(Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) - Rn(1,1)*(nref(0)*pref(0) + nref(1)*pref(1)) + A(1,1)*nref(2)*(Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) + A(1,1)*Rn(2,2)*(nref(0)*pref(0) + nref(1)*pref(1));
    C0 << A(0,0)*nref(2), nref(0)*(2*Rn(1,2)*nquery(0)*nquery(1) - Rn(0,2)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,2)*nquery(0)*nquery(2)) - nref(2)*(2*Rn(1,0)*nquery(0)*nquery(1) - Rn(0,0)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(0)*nquery(2)) + A(0,0)*nref(2)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)), nref(2)*(2*Rn(1,0)*nquery(2) - 2*Rn(2,0)*nquery(1)) - nref(0)*(2*Rn(1,2)*nquery(2) - 2*Rn(2,2)*nquery(1)) - A(0,0)*nref(2)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0)), Rn(0,2)*nref(0) - Rn(0,0)*nref(2) + A(0,0)*Rn(2,2)*nref(2),
    A(1,0)*nref(2), nref(0)*(2*Rn(0,2)*nquery(0)*nquery(1) - Rn(1,2)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,2)*nquery(1)*nquery(2)) - nref(2)*(2*Rn(0,0)*nquery(0)*nquery(1) - Rn(1,0)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,0)*nquery(1)*nquery(2)) + A(1,0)*nref(2)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)), nref(0)*(2*Rn(0,2)*nquery(2) - 2*Rn(2,2)*nquery(0)) - nref(2)*(2*Rn(0,0)*nquery(2) - 2*Rn(2,0)*nquery(0)) - A(1,0)*nref(2)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0)), Rn(1,2)*nref(0) - Rn(1,0)*nref(2) + A(1,0)*Rn(2,2)*nref(2),
    A(0,1)*nref(2), nref(1)*(2*Rn(1,2)*nquery(0)*nquery(1) - Rn(0,2)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,2)*nquery(0)*nquery(2)) - nref(2)*(2*Rn(1,1)*nquery(0)*nquery(1) - Rn(0,1)*(2*nquery(1)*nquery(1) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(0)*nquery(2)) + A(0,1)*nref(2)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)), nref(2)*(2*Rn(1,1)*nquery(2) - 2*Rn(2,1)*nquery(1)) - nref(1)*(2*Rn(1,2)*nquery(2) - 2*Rn(2,2)*nquery(1)) - A(0,1)*nref(2)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0)), Rn(0,2)*nref(1) - Rn(0,1)*nref(2) + A(0,1)*Rn(2,2)*nref(2),
    A(1,1)*nref(2), nref(1)*(2*Rn(0,2)*nquery(0)*nquery(1) - Rn(1,2)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,2)*nquery(1)*nquery(2)) - nref(2)*(2*Rn(0,1)*nquery(0)*nquery(1) - Rn(1,1)*(2*nquery(0)*nquery(0) + 2*nquery(2)*nquery(2) - 1) + 2*Rn(2,1)*nquery(1)*nquery(2)) + A(1,1)*nref(2)*(2*Rn(0,2)*nquery(0)*nquery(2) - Rn(2,2)*(2*nquery(0)*nquery(0) + 2*nquery(1)*nquery(1) - 1) + 2*Rn(1,2)*nquery(1)*nquery(2)), nref(1)*(2*Rn(0,2)*nquery(2) - 2*Rn(2,2)*nquery(0)) - nref(2)*(2*Rn(0,1)*nquery(2) - 2*Rn(2,1)*nquery(0)) - A(1,1)*nref(2)*(2*Rn(0,2)*nquery(1) - 2*Rn(1,2)*nquery(0)), Rn(1,2)*nref(1) - Rn(1,1)*nref(2) + A(1,1)*Rn(2,2)*nref(2);

    Eigen::Matrix<double,8,8> CA, CB;
    CA << Eigen::Matrix4d::Zero(), Eigen::Matrix4d::Identity(),
    -C0, -C1;
    CB << Eigen::Matrix4d::Identity(), Eigen::Matrix4d::Zero(),
    Eigen::Matrix4d::Zero(), C2;

    Eigen::GeneralizedEigenSolver< Eigen::Matrix<double,8,8> > ges(CA,CB);

    output->clear();

    for ( int i = 0; i < 8; i++ )
    {
        if ( ges.betas()[i] == 0 || ges.alphas()[i].imag() != 0 ) continue;

        double invf = ges.alphas()[i].real() / ges.betas()[i];

        Eigen::Matrix<double,8,1> soln = ges.eigenvectors().col(i).real();

        double t3 = soln(0)/soln(3);
        double q = soln(2)/soln(3);

        double theta = 2*atan(q);

        Eigen::Matrix3d Rq = Eigen::Quaterniond(Eigen::AngleAxisd(theta,nquery)).toRotationMatrix();
        Eigen::Matrix3d R = Rq*Rn;

        Eigen::Matrix<double,2,3> Ct;
        Ct << 1, 0, R(0,2) - invf*pquery(0)*(R(2,2) + t3 + R(2,0)*invf*pref(0) + R(2,1)*invf*pref(1)) + R(0,0)*invf*pref(0) + R(0,1)*invf*pref(1),
        0, 1, R(1,2) - invf*pquery(1)*(R(2,2) + t3 + R(2,0)*invf*pref(0) + R(2,1)*invf*pref(1)) + R(1,0)*invf*pref(0) + R(1,1)*invf*pref(1);

        // compute null vector of Ct
        Eigen::JacobiSVD< Eigen::Matrix<double,2,3> > svdCt(Ct, Eigen::ComputeFullV);
        Eigen::Vector3d t12sol = svdCt.matrixV().col(2);

        Eigen::Vector3d t( t12sol(0)/t12sol(2), t12sol(1)/t12sol(2), t3 );

        double fsoln = 1/invf;
        Eigen::Matrix3d Rsoln = R;
        Eigen::Vector3d tsoln = t.normalized();

        Camera calib = Camera(SimplePinholeCameraModel::model_id, std::vector<double>{fsoln, 0.0, 0.0}, -1, -1);
        CameraPose pose(Rsoln,tsoln);

        if ( check_cheirality(pose, x1, x2 ) )
        {
            output->emplace_back(ImagePair(pose,calib,calib));
        }
    }

    // std::vector<double> f_solns;
    // fpoly.realRoots(f_solns);
    
    // int nsolns = f_solns.size();
    // output->clear();
    // output->reserve(nsolns);
    // for ( int i = 0; i < nsolns; i++ )
    // {
    //     if ( f_solns[i] == 0 ) continue;

    //     double f = f_solns[i];

    //     Eigen::Matrix4d C;
    //     C <<
    //     C1_1.eval(f), C1_2.eval(f), C1_3.eval(f), C1_4.eval(f),
    //     C2_1.eval(f), C2_2.eval(f), C2_3.eval(f), C2_4.eval(f),
    //     C3_1.eval(f), C3_2.eval(f), C3_3.eval(f), C3_4.eval(f),
    //     C4_1.eval(f), C4_2.eval(f), C4_3.eval(f), C4_4.eval(f);

    //     // compute null vector of C
    //     Eigen::JacobiSVD< Eigen::Matrix4d > svdC(C, Eigen::ComputeFullV);
    //     Eigen::Vector4d p = svdC.matrixV().col(3);

    //     double t3 = p(0)/p(3);
    //     double q = p(2)/p(3);
    //     double theta = 2*atan(q);

    //     Eigen::Matrix3d Rq = Eigen::Quaterniond(Eigen::AngleAxisd(theta,r)).toRotationMatrix();
    //     Eigen::Matrix3d R = Rq*Rn;

    //     Eigen::Matrix<double,2,3> Ct;
    //     Ct <<f, 0, f*(R(0,2)*f + R(0,0)*pref(0) + R(0,1)*pref(1)) - pquery(0)*(t3 + R(2,2)*f + R(2,0)*pref(0) + R(2,1)*pref(1)),
    //     0, f, f*(R(1,2)*f + R(1,0)*pref(0) + R(1,1)*pref(1)) - pquery(1)*(t3 + R(2,2)*f + R(2,0)*pref(0) + R(2,1)*pref(1));

    //     // compute null vector of Ct
    //     Eigen::JacobiSVD< Eigen::Matrix<double,2,3> > svdCt(Ct, Eigen::ComputeFullV);
    //     Eigen::Vector3d t12sol = svdCt.matrixV().col(2);

    //     Eigen::Vector3d t( t12sol(0)/t12sol(2), t12sol(1)/t12sol(2), t3 );

    //     Eigen::Matrix3d Rsoln = R;
    //     Eigen::Vector3d tsoln = t.normalized();

    //     Camera calib = Camera(SimplePinholeCameraModel::model_id, std::vector<double>{f, 0.0, 0.0}, -1, -1);
    //     CameraPose pose(Rsoln,tsoln);

    //     if ( check_cheirality(pose, x1, x2 ) )
    //     {
    //         output->emplace_back(ImagePair(pose,calib,calib));
    //     }
    // }

    return output->size();
}

} // namespace poselib