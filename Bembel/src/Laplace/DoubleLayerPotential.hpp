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
#ifndef BEMBEL_SRC_LAPLACE_DOUBLELAYERPOTENTIAL_HPP_
#define BEMBEL_SRC_LAPLACE_DOUBLELAYERPOTENTIAL_HPP_

namespace Bembel {
// forward declaration of class LaplaceDoubleLayerPotential in order to define
// traits
template <typename LinOp>
class LaplaceDoubleLayerPotential;

template <typename LinOp>
struct PotentialTraits<LaplaceDoubleLayerPotential<LinOp>> {
  typedef Eigen::VectorXd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup Laplace
 */
template <typename LinOp>
class LaplaceDoubleLayerPotential
    : public PotentialBase<LaplaceDoubleLayerPotential<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceDoubleLayerPotential() {}
  Eigen::Matrix<
      typename PotentialReturnScalar<
          typename LinearOperatorTraits<LinOp>::Scalar, double>::Scalar,
      1, 1>
  evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);

    // integrand without basis functions
    auto integrand = evaluateKernelGrad(point, p.get_f()).dot(p.get_normal()) *
                     cauchy_value * p.get_w();

    return integrand;
  }

  /**
   * \brief Gradient of fundamental solution of Laplace problem
   */
  Eigen::Vector3d evaluateKernelGrad(const Eigen::Vector3d &x,
                                     const Eigen::Vector3d &y) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    return c / 4. / BEMBEL_PI / r3;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LAPLACE_DOUBLELAYERPOTENTIAL_HPP_
