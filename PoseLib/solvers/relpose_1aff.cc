#include "relpose_5pt.h"
#include "relpose_common.h"

#include "p1ac.h"

#include "PoseLib/misc/essential.h"
#include "PoseLib/misc/sturm.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

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
    // M << 1, 0, -pquery(0), pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) - pquery(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2), 2*Rn(2,2)*r(1) - 2*Rn(1,2)*r(2) + pquery(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) - pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) - pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1)), Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1) - pquery(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 1, -pquery(1), pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) - pquery(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2), 2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pquery(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0)), Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1) - pquery(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 0, A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), nref(0)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) - pquery(0)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) + A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), (2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - pquery(0)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0)))) - A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), nref(0)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - pquery(0)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - Rn(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + A(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    // 0, 0, A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), nref(1)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) - pquery(0)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), (2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - pquery(0)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0)))) - A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), nref(1)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) - pquery(0)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - Rn(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + A(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1));

    const double sqrtdetA = sqrt(A.determinant());
    const double cref = sqrt(2)/2;
    const double sref = cref;
    const double cquery = (A(0,0)*cref + A(0,1)*sref)/sqrtdetA;
    const double squery = (A(1,0)*cref + A(1,1)*sref)/sqrtdetA;
    // version using sref, cref, cquery, squery, sqrtdetA
    M << 1, 0, -pquery(0), pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) - pquery(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2), 2*Rn(2,2)*r(1) - 2*Rn(1,2)*r(2) + pquery(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) - pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) - pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1)), Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1) - pquery(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 1, -pquery(1), pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) - pquery(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2), 2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pquery(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0)), Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1) - pquery(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 0, -cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), cref*(pquery(0)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - nref(0)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2))) + sref*(pquery(0)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - nref(1)*(pref(0)*(2*Rn(1,0)*r(0)*r(1) - Rn(0,0)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(0)*r(2)) - Rn(0,2)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2)) + 2*Rn(1,2)*r(0)*r(1) + 2*Rn(2,2)*r(0)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(1,1)*r(0)*r(1) - Rn(0,1)*(2*r(1)*r(1) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(0)*r(2))) - cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), cref*(nref(0)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - (2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + pquery(0)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sref*(nref(1)*(2*Rn(1,2)*r(2) - 2*Rn(2,2)*r(1) + pref(0)*(2*Rn(1,0)*r(2) - 2*Rn(2,0)*r(1)) + pref(1)*(2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))) - (2*Rn(1,1)*r(2) - 2*Rn(2,1)*r(1))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) + pquery(0)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), cref*(pquery(0)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(0)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) + Rn(0,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) + sref*(pquery(0)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(1)*(Rn(0,2) + Rn(0,0)*pref(0) + Rn(0,1)*pref(1)) + Rn(0,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - cquery*sqrtdetA*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)),
    0, 0, -sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)), cref*(pquery(1)*(nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2))) - nref(0)*(pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2))) + sref*(pquery(1)*(nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)) - (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2))) - nref(1)*(pref(0)*(2*Rn(0,0)*r(0)*r(1) - Rn(1,0)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,0)*r(1)*r(2)) - Rn(1,2)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(1) + 2*Rn(2,2)*r(1)*r(2)) + (nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,1)*r(0)*r(1) - Rn(1,1)*(2*r(0)*r(0) + 2*r(2)*r(2) - 1) + 2*Rn(2,1)*r(1)*r(2))) - sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(pref(0)*(2*Rn(0,0)*r(0)*r(2) - Rn(2,0)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,0)*r(1)*r(2)) - Rn(2,2)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + pref(1)*(2*Rn(0,1)*r(0)*r(2) - Rn(2,1)*(2*r(0)*r(0) + 2*r(1)*r(1) - 1) + 2*Rn(1,1)*r(1)*r(2)) + 2*Rn(0,2)*r(0)*r(2) + 2*Rn(1,2)*r(1)*r(2)), cref*((2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))) + pquery(1)*((2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(0)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sref*((2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(2) - 2*Rn(2,2)*r(0) + pref(0)*(2*Rn(0,0)*r(2) - 2*Rn(2,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(2) - 2*Rn(2,1)*r(0))) + pquery(1)*((2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1)) - nref(1)*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))))) + sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(2*Rn(0,2)*r(1) - 2*Rn(1,2)*r(0) + pref(0)*(2*Rn(0,0)*r(1) - 2*Rn(1,0)*r(0)) + pref(1)*(2*Rn(0,1)*r(1) - 2*Rn(1,1)*r(0))), cref*(pquery(1)*(nref(0)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(0)*(Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) + Rn(1,0)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) + sref*(pquery(1)*(nref(1)*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1)) - Rn(2,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - nref(1)*(Rn(1,2) + Rn(1,0)*pref(0) + Rn(1,1)*pref(1)) + Rn(1,1)*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))) - sqrtdetA*squery*(nref(2) + nref(0)*pref(0) + nref(1)*pref(1))*(Rn(2,2) + Rn(2,0)*pref(0) + Rn(2,1)*pref(1));
    
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

        Eigen::Vector3d mon( q*q, q, 1 );
        Eigen::Vector3d t = -C.block<3,3>(0,0)*mon;

        Eigen::Matrix3d Rsoln = Rq*Rn;
        Eigen::Vector3d tsoln = t.normalized();

        output->push_back(CameraPose(Rsoln,tsoln));
    }

    return nsolns;
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

} // namespace poselib