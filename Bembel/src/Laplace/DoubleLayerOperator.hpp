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
#ifndef BEMBEL_SRC_LAPLACE_DOUBLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_LAPLACE_DOUBLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class LaplaceDoubleLayerOperator in order to define
// traits
class LaplaceDoubleLayerOperator;

template <>
struct LinearOperatorTraits<LaplaceDoubleLayerOperator> {
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
class LaplaceDoubleLayerOperator
    : public LinearOperatorBase<LaplaceDoubleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceDoubleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<LaplaceDoubleLayerOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    // integrand without basis functions
    auto integrand =
        evaluateKernelGrad(p1.get_f(), p2.get_f(), p2.get_normal()) *
        p1.get_surface_measure() * p1.get_w() * p2.get_w();

    // multiply basis functions with integrand and add to intval
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_f(),
                                          p2.get_f());

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernelGrad(p1.get_f(), p2.get_f(), p2.get_normal()) *
                p1.get_surface_measure();
    return intval;
  }

  /**
   * \brief Gradient of fundamental solution of Laplace problem
   */
  double evaluateKernelGrad(const Eigen::Vector3d &x, const Eigen::Vector3d &y,
                            const Eigen::Vector3d &y_n) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    return c.dot(y_n) / 4. / BEMBEL_PI / r3;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LAPLACE_DOUBLELAYEROPERATOR_HPP_
