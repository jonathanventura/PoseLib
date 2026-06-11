// Author: Jonathan Ventura

#include "PoseLib/misc/re3q3.h"
#include "p1ac.h"

namespace poselib {

void cayley(const Eigen::Vector3d &r, Eigen::Matrix3d &R, double &s ) {
    const double x = r(0);
    const double y = r(1);
    const double z = r(2);
    
    const double w = 1;
    const double xx = x*x;
    const double xy = x*y;
    const double xz = x*z;
    const double xw = x*w;
    const double yy = y*y;
    const double yz = y*z;
    const double yw = y*w;
    const double zz = z*z;
    const double zw = z*w;
    const double ww = w*w;

    s = w*w+x*x+y*y+z*z;
    R << ww+xx-yy-zz, 2*(xy-zw), 2*(yw+xz),
            2*(xy+zw), ww-xx+yy-zz, 2*(yz-xw),
            2*(xz-yw), 2*(xw+yz), ww-xx-yy+zz;
}

int p1ac(const Eigen::Vector3d &y, const Eigen::Vector3d &X, const Eigen::Vector3d &n, const Eigen::Matrix2d &A, std::vector<CameraPose> *output) {
    if (output == nullptr) {
        return 0;
    }
    output->clear();
    output->reserve(4);

    Eigen::Matrix<double,6,13> M;

    const Eigen::Vector2d x(X(0)/X(2),X(1)/X(2));
    const double d = X(2);
    
    // terms: [ t1, t2, t3, r1^2, r1*r2, r1*r3, r2^2, r2*r3, r3^2, r1, r2, r3, 1]
    M << -1, 0, y(0), - d*x(0) - d*y(0), -2*d*x(1), 2*d*x(0)*y(0) - 2*d, d*x(0) - d*y(0), 2*d*x(1)*y(0), d*x(0) + d*y(0), 2*d*x(1)*y(0), - 2*d - 2*d*x(0)*y(0), 2*d*x(1), d*y(0) - d*x(0),
        0, -1, y(1), d*x(1) - d*y(1), -2*d*x(0), 2*d*x(0)*y(1), - d*x(1) - d*y(1), 2*d*x(1)*y(1) - 2*d, d*x(1) + d*y(1), 2*d + 2*d*x(1)*y(1), -2*d*x(0)*y(1), -2*d*x(0), d*y(1) - d*x(1),
        0, 0, A(0,0)*n(2) + A(0,0)*n(0)*x(0) + A(0,0)*n(1)*x(1), d*n(0)*y(0) - A(0,0)*d*n(2) - d*n(1)*x(1) - d*n(2) - A(0,0)*d*n(0)*x(0) - A(0,0)*d*n(1)*x(1), 2*d*n(0)*x(1), 2*d*n(0) + 2*d*n(2)*y(0) + 2*d*n(1)*x(1)*y(0) + 2*A(0,0)*d*n(0)*x(0)*x(0) + 2*A(0,0)*d*n(2)*x(0) + 2*A(0,0)*d*n(1)*x(0)*x(1), d*n(2) - A(0,0)*d*n(2) + d*n(1)*x(1) + d*n(0)*y(0) - A(0,0)*d*n(0)*x(0) - A(0,0)*d*n(1)*x(1), 2*A(0,0)*d*n(1)*x(1)*x(1) - 2*d*n(0)*x(1)*y(0) + 2*A(0,0)*d*n(2)*x(1) + 2*A(0,0)*d*n(0)*x(0)*x(1), d*n(2) + A(0,0)*d*n(2) + d*n(1)*x(1) - d*n(0)*y(0) + A(0,0)*d*n(0)*x(0) + A(0,0)*d*n(1)*x(1), 2*A(0,0)*d*n(1)*x(1)*x(1) - 2*d*n(0)*x(1)*y(0) + 2*A(0,0)*d*n(2)*x(1) + 2*A(0,0)*d*n(0)*x(0)*x(1), 2*d*n(0) - 2*d*n(2)*y(0) - 2*d*n(1)*x(1)*y(0) - 2*A(0,0)*d*n(0)*x(0)*x(0) - 2*A(0,0)*d*n(2)*x(0) - 2*A(0,0)*d*n(1)*x(0)*x(1), -2*d*n(0)*x(1), A(0,0)*d*n(2) - d*n(2) - d*n(1)*x(1) - d*n(0)*y(0) + A(0,0)*d*n(0)*x(0) + A(0,0)*d*n(1)*x(1),
        0, 0, A(1,0)*n(2) + A(1,0)*n(0)*x(0) + A(1,0)*n(1)*x(1), d*n(0)*y(1) - d*n(0)*x(1) - A(1,0)*d*n(2) - A(1,0)*d*n(0)*x(0) - A(1,0)*d*n(1)*x(1), - 2*d*n(2) - 2*d*n(1)*x(1), 2*d*n(2)*y(1) + 2*d*n(1)*x(1)*y(1) + 2*A(1,0)*d*n(0)*x(0)*x(0) + 2*A(1,0)*d*n(2)*x(0) + 2*A(1,0)*d*n(1)*x(0)*x(1), d*n(0)*x(1) - A(1,0)*d*n(2) + d*n(0)*y(1) - A(1,0)*d*n(0)*x(0) - A(1,0)*d*n(1)*x(1), 2*d*n(0) - 2*d*n(0)*x(1)*y(1) + 2*A(1,0)*d*n(1)*x(1)*x(1) + 2*A(1,0)*d*n(2)*x(1) + 2*A(1,0)*d*n(0)*x(0)*x(1), A(1,0)*d*n(2) - d*n(0)*x(1) - d*n(0)*y(1) + A(1,0)*d*n(0)*x(0) + A(1,0)*d*n(1)*x(1), 2*A(1,0)*d*n(1)*x(1)*x(1) - 2*d*n(0)*x(1)*y(1) - 2*d*n(0) + 2*A(1,0)*d*n(2)*x(1) + 2*A(1,0)*d*n(0)*x(0)*x(1), - 2*d*n(2)*y(1) - 2*d*n(1)*x(1)*y(1) - 2*A(1,0)*d*n(0)*x(0)*x(0) - 2*A(1,0)*d*n(2)*x(0) - 2*A(1,0)*d*n(1)*x(0)*x(1), - 2*d*n(2) - 2*d*n(1)*x(1), A(1,0)*d*n(2) + d*n(0)*x(1) - d*n(0)*y(1) + A(1,0)*d*n(0)*x(0) + A(1,0)*d*n(1)*x(1),
        0, 0, A(0,1)*n(2) + A(0,1)*n(0)*x(0) + A(0,1)*n(1)*x(1), d*n(1)*x(0) - A(0,1)*d*n(2) + d*n(1)*y(0) - A(0,1)*d*n(0)*x(0) - A(0,1)*d*n(1)*x(1), - 2*d*n(2) - 2*d*n(0)*x(0), 2*d*n(1) - 2*d*n(1)*x(0)*y(0) + 2*A(0,1)*d*n(0)*x(0)*x(0) + 2*A(0,1)*d*n(2)*x(0) + 2*A(0,1)*d*n(1)*x(0)*x(1), d*n(1)*y(0) - d*n(1)*x(0) - A(0,1)*d*n(2) - A(0,1)*d*n(0)*x(0) - A(0,1)*d*n(1)*x(1), 2*d*n(2)*y(0) + 2*d*n(0)*x(0)*y(0) + 2*A(0,1)*d*n(1)*x(1)*x(1) + 2*A(0,1)*d*n(2)*x(1) + 2*A(0,1)*d*n(0)*x(0)*x(1), A(0,1)*d*n(2) - d*n(1)*x(0) - d*n(1)*y(0) + A(0,1)*d*n(0)*x(0) + A(0,1)*d*n(1)*x(1), 2*d*n(2)*y(0) + 2*d*n(0)*x(0)*y(0) + 2*A(0,1)*d*n(1)*x(1)*x(1) + 2*A(0,1)*d*n(2)*x(1) + 2*A(0,1)*d*n(0)*x(0)*x(1), 2*d*n(1) + 2*d*n(1)*x(0)*y(0) - 2*A(0,1)*d*n(0)*x(0)*x(0) - 2*A(0,1)*d*n(2)*x(0) - 2*A(0,1)*d*n(1)*x(0)*x(1), 2*d*n(2) + 2*d*n(0)*x(0), A(0,1)*d*n(2) + d*n(1)*x(0) - d*n(1)*y(0) + A(0,1)*d*n(0)*x(0) + A(0,1)*d*n(1)*x(1),
        0, 0, A(1,1)*n(2) + A(1,1)*n(0)*x(0) + A(1,1)*n(1)*x(1), d*n(2) - A(1,1)*d*n(2) + d*n(0)*x(0) + d*n(1)*y(1) - A(1,1)*d*n(0)*x(0) - A(1,1)*d*n(1)*x(1), 2*d*n(1)*x(0), 2*A(1,1)*d*n(0)*x(0)*x(0) - 2*d*n(1)*x(0)*y(1) + 2*A(1,1)*d*n(2)*x(0) + 2*A(1,1)*d*n(1)*x(0)*x(1), d*n(1)*y(1) - A(1,1)*d*n(2) - d*n(0)*x(0) - d*n(2) - A(1,1)*d*n(0)*x(0) - A(1,1)*d*n(1)*x(1), 2*d*n(1) + 2*d*n(2)*y(1) + 2*d*n(0)*x(0)*y(1) + 2*A(1,1)*d*n(1)*x(1)*x(1) + 2*A(1,1)*d*n(2)*x(1) + 2*A(1,1)*d*n(0)*x(0)*x(1), d*n(2) + A(1,1)*d*n(2) + d*n(0)*x(0) - d*n(1)*y(1) + A(1,1)*d*n(0)*x(0) + A(1,1)*d*n(1)*x(1), 2*d*n(2)*y(1) - 2*d*n(1) + 2*d*n(0)*x(0)*y(1) + 2*A(1,1)*d*n(1)*x(1)*x(1) + 2*A(1,1)*d*n(2)*x(1) + 2*A(1,1)*d*n(0)*x(0)*x(1), 2*d*n(1)*x(0)*y(1) - 2*A(1,1)*d*n(0)*x(0)*x(0) - 2*A(1,1)*d*n(2)*x(0) - 2*A(1,1)*d*n(1)*x(0)*x(1), 2*d*n(1)*x(0), A(1,1)*d*n(2) - d*n(2) - d*n(0)*x(0) - d*n(1)*y(1) + A(1,1)*d*n(0)*x(0) + A(1,1)*d*n(1)*x(1);

    // G-J elimination
    Eigen::Matrix<double,6,10> G = M.block<6,6>(0,0).partialPivLu().solve(M.block<6,10>(0,3));
    Eigen::Matrix<double,3,10> T = G.block<3,10>(0,0);
    Eigen::Matrix<double,3,10> C = G.block<3,10>(3,0);
    
    Eigen::Matrix<double,3,8> cayley_solutions;
    int nsolns = re3q3::re3q3(C, &cayley_solutions);
    
    for ( int i = 0; i < nsolns; i++ )
    {
        Eigen::Vector3d r = cayley_solutions.col(i);
        Eigen::Matrix3d R;
        double s;
        cayley(r,R,s);
        
        Eigen::Matrix<double,10,1> X;
        X << r(0)*r(0), r(0)*r(1), r(0)*r(2), r(1)*r(1), r(1)*r(2), r(2)*r(2), r(0), r(1), r(2), 1;
        Eigen::Vector3d t = -T*X;

        R /= s;
        t *= d/s;
        
        output->emplace_back(R,t);
    }

    return output->size();
}

} // namespace poselib
