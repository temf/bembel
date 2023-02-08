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
#ifndef BEMBEL_SRC_HELMHOLTZ_SINGLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_HELMHOLTZ_SINGLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class HelmholtzSingleLayerOperator in order to define
// traits
class HelmholtzSingleLayerOperator;

template <>
struct LinearOperatorTraits<HelmholtzSingleLayerOperator> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
  enum {
    OperatorOrder = -1,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup Helmholtz
 */
class HelmholtzSingleLayerOperator
    : public LinearOperatorBase<HelmholtzSingleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  HelmholtzSingleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<HelmholtzSingleLayerOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    // integrand without basis functions
    std::complex<double> integrand =
        evaluateKernel(p1.get_f(), p2.get_f()) * p1.get_surface_measure() *
        p2.get_surface_measure() * p1.get_w() * p2.get_w();

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    // (*intval) += super_space.BasisInteraction(s, t) * evaluateKernel(x_f,
    // y_f)
    // * x_kappa * y_kappa * ws * wt;
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_xi(),
                                          p2.get_xi());

    return;
  }

  Eigen::Matrix<std::complex<double>, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // interpolation
    Eigen::Matrix<std::complex<double>, 1, 1> intval;
    intval(0) = evaluateKernel(p1.get_f(), p2.get_f()) *
                p1.get_surface_measure() * p2.get_surface_measure();

    return intval;
  }

  /**
   * \brief Fundamental solution of Helmholtz problem
   */
  std::complex<double> evaluateKernel(const Eigen::Vector3d &x,
                                      const Eigen::Vector3d &y) const {
    auto r = (x - y).norm();
    return std::exp(-std::complex<double>(0., 1.) * wavenumber_ * r) / 4. /
           BEMBEL_PI / r;
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
#endif  // BEMBEL_SRC_HELMHOLTZ_SINGLELAYEROPERATOR_HPP_
