// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_SRC_LAPLACE_SINGLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_LAPLACE_SINGLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class LaplaceSingleLayerOperator in order to define
// traits
class LaplaceSingleLayerOperator;
/**
 * \brief Specification of the LinerOperatorTraits for Laplace.
 */
template <>
struct LinearOperatorTraits<LaplaceSingleLayerOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = -1,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup Laplace
 * \brief This class implements the specification of the integration for the
 * single layer operator for Laplace.
 */
class LaplaceSingleLayerOperator
    : public LinearOperatorBase<LaplaceSingleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceSingleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<LaplaceSingleLayerOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    // integrand without basis functions
    double integrand = evaluateKernel(p1.get_f(), p2.get_f()) *
                       p1.get_surface_measure() * p2.get_surface_measure() *
                       p1.get_w() * p2.get_w();

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    // (*intval) += super_space.basisInteraction(s, t) * evaluateKernel(x_f,
    // y_f)
    // * x_kappa * y_kappa * ws * wt;
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_xi(),
                                          p2.get_xi());

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernel(p1.get_f(), p2.get_f()) *
                p1.get_surface_measure() * p2.get_surface_measure();

    return intval;
  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
  double evaluateKernel(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y) const {
    return 1. / 4. / BEMBEL_PI / (x - y).norm();
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LAPLACE_SINGLELAYEROPERATOR_HPP_
