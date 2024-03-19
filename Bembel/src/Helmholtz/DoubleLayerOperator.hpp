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
#ifndef BEMBEL_SRC_HELMHOLTZ_DOUBLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_HELMHOLTZ_DOUBLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class HelmholtzDoubleLayerOperator in order to define
// traits
class HelmholtzDoubleLayerOperator;

template <>
struct LinearOperatorTraits<HelmholtzDoubleLayerOperator> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
  enum {
    OperatorOrder = 0,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup Helmholtz
 */
class HelmholtzDoubleLayerOperator
    : public LinearOperatorBase<HelmholtzDoubleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  HelmholtzDoubleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<HelmholtzDoubleLayerOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    // integrand without basis functions
    auto integrand =
        evaluateKernelGrad(p1.get_f(), p1.get_normal(), p2.get_f()) *
        p1.get_surface_measure() * p1.get_w() * p2.get_w();

    // multiply basis functions with integrand and add to intval
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_xi(),
                                          p2.get_xi());

    return;
  }

  Eigen::Matrix<std::complex<double>, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // interpolation
    Eigen::Matrix<std::complex<double>, 1, 1> intval;
    intval(0) = evaluateKernelGrad(p1.get_f(), p1.get_normal(), p2.get_f()) *
                p2.get_surface_measure();
    return intval;
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
#endif  // BEMBEL_SRC_HELMHOLTZ_DOUBLELAYEROPERATOR_HPP_
