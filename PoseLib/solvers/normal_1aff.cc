// Author: Jonathan Ventura

#include "normal_1aff.h"

namespace poselib {

int normal_1aff(const Eigen::Vector3d &x, const Eigen::Vector3d &X, const Eigen::Matrix2d &A, const Eigen::Matrix3d &R, const Eigen::Vector3d &t, Eigen::Vector3d *output ) {
    if (output == nullptr) {
        return 0;
    }

    const Eigen::Vector2d x1(X(0)/X(2),X(1)/X(2));
    const double d = X(2);
    const Eigen::Vector2d x2(x(0)/x(2),x(1)/x(2));

    Eigen::Matrix<double,4,3> M;
    
    // terms: [ n1 n2 n3 ]
    M << d*(R(0,2) + R(0,1)*x1(1) - x2(0)*(R(2,2) + R(2,1)*x1(1))) + A(0,0)*x1(0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))), A(0,0)*x1(1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(0,0)*x1(1) - R(2,0)*x1(1)*x2(0)), A(0,0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(0,0) - R(2,0)*x2(0)),
    d*(R(1,2) + R(1,1)*x1(1) - x2(1)*(R(2,2) + R(2,1)*x1(1))) + A(1,0)*x1(0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))), A(1,0)*x1(1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(1,0)*x1(1) - R(2,0)*x1(1)*x2(1)), A(1,0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(1,0) - R(2,0)*x2(1)), A(0,1)*x1(0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(0,1)*x1(0) - R(2,1)*x1(0)*x2(0)), d*(R(0,2) + R(0,0)*x1(0) - x2(0)*(R(2,2) + R(2,0)*x1(0))) + A(0,1)*x1(1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))), A(0,1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(0,1) - R(2,1)*x2(0)),
     A(1,1)*x1(0)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(1,1)*x1(0) - R(2,1)*x1(0)*x2(1)), d*(R(1,2) + R(1,0)*x1(0) - x2(1)*(R(2,2) + R(2,0)*x1(0))) + A(1,1)*x1(1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))), A(1,1)*(t(2) + d*(R(2,2) + R(2,0)*x1(0) + R(2,1)*x1(1))) - d*(R(1,1) - R(2,1)*x2(1));

    // solve for normal using SVD
    Eigen::JacobiSVD<Eigen::Matrix<double,4,3>> svd(M, Eigen::ComputeFullV);
    *output = svd.matrixV().col(2);

    return 1;
}

} // namespace poselib
