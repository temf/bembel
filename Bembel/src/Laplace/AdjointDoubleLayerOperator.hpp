// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_SRC_LAPLACE_ADJOINTDOUBLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_LAPLACE_ADJOINTDOUBLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class LaplaceAdjointDoubleLayerOperator in order to
// define traits
class LaplaceAdjointDoubleLayerOperator;

template <>
struct LinearOperatorTraits<LaplaceAdjointDoubleLayerOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = 0,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup Laplace
 */
class LaplaceAdjointDoubleLayerOperator
    : public LinearOperatorBase<LaplaceAdjointDoubleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceAdjointDoubleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<typename LinearOperatorTraits<
                        LaplaceAdjointDoubleLayerOperator>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get quadrature weights
    auto ws = p1(2);
    auto wt = p2(2);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute surface measures from tangential derivatives
    auto x_n = x_f_dx.cross(x_f_dy);
    auto y_kappa = y_f_dx.cross(y_f_dy).norm();

    // integrand without basis functions
    auto integrand = evaluateKernelGrad(x_f, x_n, y_f) * y_kappa * ws * wt;

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    // (*intval) += super_space.BasisInteraction(s, t) * evaluateKernel(x_f,
    // y_f)
    // * x_kappa * y_kappa * ws * wt;
    super_space.addScaledBasisInteraction(intval, integrand, s, t);

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute surface measure in x and unnormalized normal in y from tangential
    // derivatives
    auto x_n = x_f_dx.cross(x_f_dy);
    auto y_kappa = y_f_dx.cross(y_f_dy).norm();

    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernelGrad(x_f, x_n, y_f) * y_kappa;

    return intval;
  }

  /**
   * \brief Gradient of fundamental solution of Laplace problem
   */
  double evaluateKernelGrad(const Eigen::Vector3d &x,
                            const Eigen::Vector3d &x_n,
                            const Eigen::Vector3d &y) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    return -c.dot(x_n) / 4. / BEMBEL_PI / r3;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LAPLACE_ADJOINTDOUBLELAYEROPERATOR_HPP_
