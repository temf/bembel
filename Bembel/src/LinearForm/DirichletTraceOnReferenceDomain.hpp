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
#ifndef BEMBEL_SRC_LINEARFORM_DIRICHLETTRACEONREFERENCEDOMAIN_H_
#define BEMBEL_SRC_LINEARFORM_DIRICHLETTRACEONREFERENCEDOMAIN_H_

namespace Bembel {

template <typename Scalar>
class DirichletTraceOnReferenceDomain;

template <typename ScalarT>
struct LinearFormTraits<DirichletTraceOnReferenceDomain<ScalarT>> {
  typedef ScalarT Scalar;
};

/**
 *  \ingroup LinearForm
 *  \brief This class provides an implementation of the Dirichlet trace operator
 * and a corresponding method to evaluate the linear form corresponding to the
 * right hand side of the system via quadrature.
 */
template <typename Scalar>
class DirichletTraceOnReferenceDomain
    : public LinearFormBase<DirichletTraceOnReferenceDomain<Scalar>, Scalar> {
 public:
  DirichletTraceOnReferenceDomain() {}
  void set_function(
      const std::function<Scalar(int, Eigen::Vector2d)> &function) {
    function_ = function;
  }
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *intval) const {
    // integrand without basis functions
    //std::cout << p.get_xi() << std::endl;
    Scalar integrand = function_(p.get_patch(), p.get_xi()) *
                       p.get_surface_measure() * p.get_w();
    // multiply basis functions with integrand
    super_space.addScaledBasis(intval, integrand, p.get_xi());
    return;
  }

 private:
  std::function<Scalar(int, Eigen::Vector2d)> function_;
};
}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEARFORM_DIRICHLETTRACEONREFERENCEDOMAIN_H_
