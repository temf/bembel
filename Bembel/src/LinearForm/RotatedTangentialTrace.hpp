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
#ifndef BEMBEL_SRC_LINEARFORM_ROTATEDTANGENTIALTRACE_HPP_
#define BEMBEL_SRC_LINEARFORM_ROTATEDTANGENTIALTRACE_HPP_

namespace Bembel {

template <typename Scalar>
class RotatedTangentialTrace;

template <typename ScalarT>
struct LinearFormTraits<RotatedTangentialTrace<ScalarT>> {
  typedef ScalarT Scalar;
};

/**
 *    \ingroup LinearForm
 *    \brief This class provides a specialization of the linear form required
 *for the solution of the electric field integral equation.
 **/
template <typename Scalar>
class RotatedTangentialTrace
    : public LinearFormBase<RotatedTangentialTrace<Scalar>, Scalar> {
 public:
  RotatedTangentialTrace() {}
  void set_function(
      const std::function<Eigen::Matrix<Scalar, 3, 1>(Eigen::Vector3d)>
          &function) {
    function_ = function;
  }
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 2> *intval) const {
    // tangential component + quadrature weights
    // use n x f x n = f-<f,n>n to avoid troubles with -flto flag in combination
    // of .cross()
    Eigen::Matrix<Scalar, 3, 1> fun_x_f = function_(p.get_f());
    Eigen::Vector3d x_n = p.get_unit_normal();
    Eigen::Matrix<Scalar, 3, 1> tangential_component =
        (fun_x_f - fun_x_f.dot(x_n) * x_n) * p.get_w();

    // extract tangential component
    Eigen::Matrix<Scalar, 2, 1> components =
        p.get_jacobian().transpose() * tangential_component;

    // evaluate shape functions
    auto phiPhiVec = super_space.basis(p.get_xi());

    // compute integrals
    (*intval) += phiPhiVec * components.transpose();
    return;
  }

 private:
  std::function<Eigen::Matrix<Scalar, 3, 1>(Eigen::Vector3d)> function_;
};
}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEARFORM_ROTATEDTANGENTIALTRACE_HPP_
