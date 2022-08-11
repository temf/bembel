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
#ifndef BEMBEL_SRC_LINEARFORM_DIRICHLETTRACE_HPP_
#define BEMBEL_SRC_LINEARFORM_DIRICHLETTRACE_HPP_

namespace Bembel {

template <typename Scalar>
class DirichletTrace;

template <typename ScalarT>
struct LinearFormTraits<DirichletTrace<ScalarT>> {
  typedef ScalarT Scalar;
};

/**
 *  \ingroup LinearForm
 *  \brief This class provides an implementation of the Dirichlet trace operator
 * and a corresponding method to evaluate the linear form corresponding to the
 * right hand side of the system via quadrature.
 */
template <typename Scalar>
class DirichletTrace : public LinearFormBase<DirichletTrace<Scalar>, Scalar> {
 public:
  DirichletTrace() {}
  void set_function(const std::function<Scalar(Eigen::Vector3d)> &function) {
    function_ = function;
  }
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

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

    // integrand without basis functions
    auto integrand = function_(x_f) * x_kappa * ws;

    // multiply basis functions with integrand
    super_space.addScaledBasis(intval, integrand, s);

    return;
  }

 private:
  std::function<Scalar(Eigen::Vector3d)> function_;
};
}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEARFORM_DIRICHLETTRACE_HPP_
