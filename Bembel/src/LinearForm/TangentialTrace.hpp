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
#ifndef BEMBEL_SRC_LINEARFORM_TANGENTIALTRACE_HPP_
#define BEMBEL_SRC_LINEARFORM_TANGENTIALTRACE_HPP_

namespace Bembel {

template <typename Scalar>
class TangentialTrace;

template <typename ScalarT>
struct LinearFormTraits<TangentialTrace<ScalarT>> {
  typedef ScalarT Scalar;
};

/**
 *    \ingroup LinearForm
 *    \brief This class provides a specialization of the linear form required
 *for the solution of the electric field integral equation.
 **/
template <typename Scalar>
class TangentialTrace : public LinearFormBase<TangentialTrace<Scalar>, Scalar> {
 public:
  TangentialTrace() {}
  void set_function(
      const std::function<Eigen::Matrix<Scalar, 3, 1>(Eigen::Vector3d)>
          &function) {
    function_ = function;
  }
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> *intval) const {
    int polynomial_degree = super_space.get_polynomial_degree();
    int polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    // get evaluation points on unit square
    Eigen::Vector2d s = p.segment<2>(0);

    // get quadrature weights
    double ws = p(2);

    // get points on geometry and tangential derivatives
    Eigen::Vector3d x_f = p.segment<3>(3);
    Eigen::Vector3d x_f_dx = p.segment<3>(6);
    Eigen::Vector3d x_f_dy = p.segment<3>(9);

    // compute surface measures from tangential derivatives
    Eigen::Vector3d x_n = x_f_dx.cross(x_f_dy).normalized();

    // tangential component + quadrature weights
    Eigen::Matrix<Scalar, 3, 1> fun_x_f = function_(x_f);
    Eigen::Matrix<Scalar, 3, 1> tangential_component = fun_x_f.cross(x_n) * ws;

    // extract tangential component
    Scalar component_x = x_f_dx.dot(tangential_component);
    Scalar component_y = x_f_dy.dot(tangential_component);

    // evaluate shape functions
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> phiPhiVec =
        super_space.basis(s);

    // multiply basis functions with integrand
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2> phiPhiMat(
        polynomial_degree_plus_one_squared, 2);
    phiPhiMat.col(0) = component_x * phiPhiVec;
    phiPhiMat.col(1) = component_y * phiPhiVec;

    // compute integrals
    (*intval) += phiPhiMat;
    return;
  }

 private:
  std::function<Eigen::Matrix<Scalar, 3, 1>(Eigen::Vector3d)> function_;
};
}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEARFORM_TANGENTIALTRACE_HPP_
