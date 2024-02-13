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

#ifndef BEMBEL_SRC_UTIL_SURFACEL2ERROR_HPP_
#define BEMBEL_SRC_UTIL_SURFACEL2ERROR_HPP_

namespace Bembel {

template <typename Op, typename Derived, typename Functor>
double surfaceL2error(const AnsatzSpace<Op> &ansatz_space,
                      const Eigen::MatrixBase<Derived> &vec,
                      const Functor &functor, int deg = 4) {
  typedef typename Derived::Scalar Scalar;
  FunctionEvaluator<Op> fun_val(ansatz_space);
  fun_val.set_function(vec);
  Scalar retval = 0;
  GaussSquare<Constants::maximum_quadrature_degree> GS;
  auto Q = GS[deg];
  SurfacePoint qp;
  const auto &super_space = ansatz_space.get_superspace();
  const ElementTree &et = super_space.get_mesh().get_element_tree();
  for (auto element = et.cpbegin(); element != et.cpend(); ++element) {
    for (auto i = 0; i < Q.w_.size(); ++i) {
      super_space.map2surface(*element, Q.xi_.col(i), Q.w_(i), &qp);
      // get points on geometry and tangential derivatives
      const auto &x_f = qp.segment<3>(3);
      const auto &x_f_dx = qp.segment<3>(6);
      const auto &x_f_dy = qp.segment<3>(9);
      const auto &normal = x_f_dx.cross(x_f_dy);
      // compute surface measures from tangential derivatives
      Scalar x_kappa = normal.norm();
      const Scalar val = fun_val.evaluate(*element, qp)(0);
      retval += x_kappa * Q.w_(i) * element->get_h() * element->get_h() *
                (functor(x_f) - val) * (functor(x_f) - val);
    }
  }
  return sqrt(retval);
}

}  // namespace Bembel
#endif  // BEMBEL_SRC_UTIL_SURFACEL2ERROR_HPP_
