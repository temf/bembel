// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_LINEAROPERATOR_MAXWELL_MAXWELLSINGLELAYEROPERATOR_H_
#define BEMBEL_LINEAROPERATOR_MAXWELL_MAXWELLSINGLELAYEROPERATOR_H_

namespace Bembel {
// forward declaration of class MaxwellSingleLayerOperator in order to define
// traits
class MaxwellSingleLayerOperator;

template <>
struct LinearOperatorTraits<MaxwellSingleLayerOperator> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
  enum {
    OperatorOrder = -1,
    Form = DifferentialForm::DivConforming,
    NumberOfFMMComponents = 2
  };
};

class MaxwellSingleLayerOperator
    : public LinearOperatorBase<MaxwellSingleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  MaxwellSingleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<MaxwellSingleLayerOperator>::Scalar,
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
    auto h = 1. / (1 << super_space.get_refinement_level());  // h = 1 ./ (2^M)

    // integrand without basis functions
    auto kernel_evaluation = evaluateKernel(x_f, y_f) * ws * wt;
    auto integrand_vector = kernel_evaluation;
    auto integrand_divergence = -kernel_evaluation / wavenumber2_ / h / h;

    // vector part: mulitply shape functions with integrand and add to buffer
    super_space.addScaledVectorBasisInteraction(intval, integrand_vector, s, t,
                                                x_f_dx, x_f_dy, y_f_dx, y_f_dy);

    // divergence part: multiply shape functions with integrand and add to
    // buffer
    super_space.addScaledVectorBasisDivergenceInteraction(
        intval, integrand_divergence, s, t);

    return;
  }

  /**
   * \brief Fundamental solution of Helmholtz/Maxwell problem
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
    wavenumber2_ = wavenumber_ * wavenumber_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  std::complex<double> get_wavenumber() { return wavenumber_; }

 private:
  std::complex<double> wavenumber_;
  std::complex<double> wavenumber2_;
};

}  // namespace Bembel
#endif
