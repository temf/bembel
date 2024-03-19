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
#ifndef BEMBEL_SRC_LAPLACE_SINGLELAYERPOTENTIALGRADIENT_HPP_
#define BEMBEL_SRC_LAPLACE_SINGLELAYERPOTENTIALGRADIENT_HPP_

namespace Bembel {
// forward declaration of class LaplaceSingleLayerPotentialGradient in order
// to define traits
template <typename LinOp>
class LaplaceSingleLayerPotentialGradient;

template <typename LinOp>
struct PotentialTraits<LaplaceSingleLayerPotentialGradient<LinOp>> {
  typedef Eigen::VectorXd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 3;
};

/**
 * \ingroup Laplace
 */
template <typename LinOp>
class LaplaceSingleLayerPotentialGradient
    : public PotentialBase<LaplaceSingleLayerPotentialGradient<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceSingleLayerPotentialGradient() {}
  Eigen::Matrix<
      typename PotentialReturnScalar<
          typename LinearOperatorTraits<LinOp>::Scalar, double>::Scalar,
      3, 1>
  evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // evaluate kernel
    auto kernel = evaluateKernelGrad(point, p.get_f());

    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);

    // integrand without basis functions
    auto integrand =
        kernel * cauchy_value * p.get_surface_measure() * p.get_w();

    return integrand;
  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
  Eigen::VectorXd evaluateKernelGrad(const Eigen::Vector3d &x,
                                     const Eigen::Vector3d &y) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    return -c / r3 / 4. / BEMBEL_PI;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LAPLACE_SINGLELAYERPOTENTIALGRADIENT_HPP_
