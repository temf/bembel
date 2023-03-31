// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_HELMHOLTZ_HELMHOLTZDOUBLELAYERPOTENTIAL_H_
#define BEMBEL_LINEAROPERATOR_HELMHOLTZ_HELMHOLTZDOUBLELAYERPOTENTIAL_H_

namespace Bembel {
// forward declaration of class HelmholtzDoubleLayerPotential in order to define
// traits
template <typename LinOp>
class HelmholtzDoubleLayerPotential;

template <typename LinOp>
struct PotentialTraits<HelmholtzDoubleLayerPotential<LinOp>> {
  typedef Eigen::VectorXcd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup Helmholtz
 */
template <typename LinOp>
class HelmholtzDoubleLayerPotential
    : public PotentialBase<HelmholtzDoubleLayerPotential<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  HelmholtzDoubleLayerPotential() {}
  Eigen::Matrix<typename PotentialReturnScalar<
                    typename LinearOperatorTraits<LinOp>::Scalar,
                    std::complex<double>>::Scalar,
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
    // dot: adjoint in first variable
    auto integrand = evaluateKernelGrad(point, x_f, x_n) * cauchy_value * ws;

    return integrand;
  }

  /**
   * \brief Gradient of fundamental solution of Helmholtz problem
   */
  std::complex<double> evaluateKernelGrad(const Eigen::Vector3d &x,
                                          const Eigen::Vector3d &y,
                                          const Eigen::Vector3d &y_n) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    auto i = std::complex<double>(0., 1.);
    return c.dot(y_n) * std::exp(-i * wavenumber_ * r) *
           (1. + i * wavenumber_ * r) / 4. / BEMBEL_PI / r3;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setters
  //////////////////////////////////////////////////////////////////////////////
  void set_wavenumber(std::complex<double> wavenumber) {
    wavenumber_ = wavenumber;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  std::complex<double> get_wavenumber() { return wavenumber_; }

 private:
  std::complex<double> wavenumber_;
};

}  // namespace Bembel
#endif
