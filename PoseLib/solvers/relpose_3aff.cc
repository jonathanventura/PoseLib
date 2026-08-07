#include "relpose_3aff.h"

#include "PoseLib/misc/univariate.h"

#include <Eigen/Dense>
#include <cmath>

namespace poselib {

int relpose_3aff(const std::vector<Eigen::Vector3d> &u, const std::vector<Eigen::Vector3d> &v, const std::vector<Eigen::Matrix2d> &A,
                std::vector<Eigen::Matrix3d> *fundamental_matrices) {

    // Compute epipolar constraints
    Eigen::Matrix<double, 9, 9> epipolar_constraints;
    for (int i = 0; i < 3; ++i) {
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

    Eigen::Matrix<double, 9, 9> Q = epipolar_constraints.fullPivHouseholderQr().matrixQ();
    Eigen::Matrix<double, 9, 2> N = Q.rightCols(2);

    // coefficients for det(F(x)) = 0
    const double c3 = N(0, 0) * N(4, 0) * N(8, 0) - N(0, 0) * N(5, 0) * N(7, 0) - N(1, 0) * N(3, 0) * N(8, 0) +
                      N(1, 0) * N(5, 0) * N(6, 0) + N(2, 0) * N(3, 0) * N(7, 0) - N(2, 0) * N(4, 0) * N(6, 0);
    const double c2 = N(0, 0) * N(4, 0) * N(8, 1) + N(0, 0) * N(4, 1) * N(8, 0) - N(0, 0) * N(5, 0) * N(7, 1) -
                      N(0, 0) * N(5, 1) * N(7, 0) + N(0, 1) * N(4, 0) * N(8, 0) - N(0, 1) * N(5, 0) * N(7, 0) -
                      N(1, 0) * N(3, 0) * N(8, 1) - N(1, 0) * N(3, 1) * N(8, 0) + N(1, 0) * N(5, 0) * N(6, 1) +
                      N(1, 0) * N(5, 1) * N(6, 0) - N(1, 1) * N(3, 0) * N(8, 0) + N(1, 1) * N(5, 0) * N(6, 0) +
                      N(2, 0) * N(3, 0) * N(7, 1) + N(2, 0) * N(3, 1) * N(7, 0) - N(2, 0) * N(4, 0) * N(6, 1) -
                      N(2, 0) * N(4, 1) * N(6, 0) + N(2, 1) * N(3, 0) * N(7, 0) - N(2, 1) * N(4, 0) * N(6, 0);
    const double c1 = N(0, 0) * N(4, 1) * N(8, 1) - N(0, 0) * N(5, 1) * N(7, 1) + N(0, 1) * N(4, 0) * N(8, 1) +
                      N(0, 1) * N(4, 1) * N(8, 0) - N(0, 1) * N(5, 0) * N(7, 1) - N(0, 1) * N(5, 1) * N(7, 0) -
                      N(1, 0) * N(3, 1) * N(8, 1) + N(1, 0) * N(5, 1) * N(6, 1) - N(1, 1) * N(3, 0) * N(8, 1) -
                      N(1, 1) * N(3, 1) * N(8, 0) + N(1, 1) * N(5, 0) * N(6, 1) + N(1, 1) * N(5, 1) * N(6, 0) +
                      N(2, 0) * N(3, 1) * N(7, 1) - N(2, 0) * N(4, 1) * N(6, 1) + N(2, 1) * N(3, 0) * N(7, 1) +
                      N(2, 1) * N(3, 1) * N(7, 0) - N(2, 1) * N(4, 0) * N(6, 1) - N(2, 1) * N(4, 1) * N(6, 0);
    const double c0 = N(0, 1) * N(4, 1) * N(8, 1) - N(0, 1) * N(5, 1) * N(7, 1) - N(1, 1) * N(3, 1) * N(8, 1) +
                      N(1, 1) * N(5, 1) * N(6, 1) + N(2, 1) * N(3, 1) * N(7, 1) - N(2, 1) * N(4, 1) * N(6, 1);

    // Solve the cubic (guarded against degenerate c3 ≈ 0 case)
    double roots[3];
    int n_roots;
    if (std::abs(c3) < 1e-14) {
        // Cubic degenerates to quadratic: c2*x^2 + c1*x + c0 = 0
        n_roots = univariate::solve_quadratic_real(c2, c1, c0, roots);
    } else {
        double inv_c3 = 1.0 / c3;
        n_roots = univariate::solve_cubic_real(c2 * inv_c3, c1 * inv_c3, c0 * inv_c3, roots);
    }

    // Reshape back into 3x3 matrices
    fundamental_matrices->clear();
    fundamental_matrices->reserve(n_roots);
    for (int i = 0; i < n_roots; ++i) {
        Eigen::Matrix<double, 9, 1> f = N.col(0) * roots[i] + N.col(1);
        f.normalize();
        fundamental_matrices->push_back(Eigen::Map<Eigen::Matrix3d>(f.data()));
    }

    return n_roots;
}

} // namespace poselib