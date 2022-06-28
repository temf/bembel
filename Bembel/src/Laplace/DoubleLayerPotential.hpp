// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYERPOTENTIAL_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYERPOTENTIAL_H_

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
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // compute unnormalized normal from tangential derivatives
    auto x_n = x_f_dx.cross(x_f_dy);

    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);

    // integrand without basis functions
    auto integrand =
        evaluateKernelGrad(point, x_f).dot(x_n) * cauchy_value * ws;

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
#endif
