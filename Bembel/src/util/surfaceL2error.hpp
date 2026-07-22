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
  const auto longvec = (ansatz_space.get_transformation_matrix() * vec).eval();
  const auto &super_space = ansatz_space.get_superspace();
  const ElementTree &et = super_space.get_mesh().get_element_tree();
  for (auto element = et.cpbegin(); element != et.cpend(); ++element) {
    for (auto i = 0; i < Q.w_.size(); ++i) {
      super_space.map2surface(*element, Q.xi_.col(i), Q.w_(i), &qp);
      // integrand without basis functions
      const Scalar val = fun_val.evaluate(*element, qp)(0);
      retval += qp.get_surface_measure() * Q.w_(i) * element->get_h() *
                element->get_h() *
                (functor(qp.get_f()) - val) *
                (functor(qp.get_f()) - val);
    }
  }
  return sqrt(retval);
}

}  // namespace Bembel
#endif  // BEMBEL_SRC_UTIL_SURFACEL2ERROR_HPP_
