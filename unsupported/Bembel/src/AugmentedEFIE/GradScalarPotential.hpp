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

#ifndef BEMBEL_SRC_AUGMENTEDEFIE_GRADSCALARPOTENTIAL_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_GRADSCALARPOTENTIAL_HPP_

namespace Bembel {
// forward declaration of class GradScalarPotential in order to define
// traits
template <typename LinOp>
class GradScalarPotential;

template <typename LinOp>
struct PotentialTraits<GradScalarPotential<LinOp>> {
  typedef Eigen::VectorXcd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 3;
};

/**
 * \ingroup AugmentedEFIE
 */
template <typename LinOp>
class GradScalarPotential
    : public PotentialBase<GradScalarPotential<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  GradScalarPotential() {}
  Eigen::Matrix<typename PotentialReturnScalar<
                    typename LinearOperatorTraits<LinOp>::Scalar,
                    std::complex<double>>::Scalar,
                3, 1>
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

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();

    // evaluate gradient of kernel
    auto kernel_gradient = evaluateKernelGrad(point, x_f);

    // assemble Galerkin solution
    auto cauchy_data = fun_ev.evaluate(element, p);

    // integrand without basis functions, note that the surface measure
    // disappears for the divergence
    auto integrand = kernel_gradient * cauchy_data * ws * x_kappa;

    return integrand;
  }
  /**
   * \brief Gradient of fundamental solution of Helmholtz problem
   */
  Eigen::VectorXcd evaluateKernelGrad(const Eigen::Vector3d &x,
                                      const Eigen::Vector3d &y) const {
    auto c = x - y;
    auto r = c.norm();
    auto r3 = r * r * r;
    auto i = std::complex<double>(0., 1.);
    return (std::exp(-i * wavenumber_ * r) * (-1. - i * wavenumber_ * r) / 4. /
            BEMBEL_PI / r3) *
           c;
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
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_GRADSCALARPOTENTIAL_HPP_
