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
#ifndef BEMBEL_SRC_DUFFYTRICK_INTEGRATE1_HPP_
#define BEMBEL_SRC_DUFFYTRICK_INTEGRATE1_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *    \brief no-problem quadrature routine, elements are sufficiently far
 *           away from each other, but not far-field yet (this is just an
 *           isotropic tensor product quadrature!)
 *    \todo  be sure that map2element computes the weight h*Q.w(i) such that
 *           the integrand may then be scaled by qp1.weight * qp2.weight
 **/
template <typename Derived, class T>
void integrate1(const LinearOperatorBase<Derived> &LinOp, const T &super_space,
                const ElementTreeNode &e1, int rot1, const ElementTreeNode &e2,
                int rot2, const Eigen::MatrixXd &ffield_qnodes1,
                const Eigen::MatrixXd &ffield_qnodes2, const Cubature &Q,
                Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                              Eigen::Dynamic, Eigen::Dynamic> *intval) {
  double h = e1.get_h();
  intval->setZero();
  SurfacePoint qp1, qp2;
  for (auto i = 0; i < Q.w_.size(); ++i) {
    super_space.map2surface(e1, Q.xi_.col(i), h * Q.w_(i), &qp1);
    for (auto j = 0; j < Q.w_.size(); ++j) {
      super_space.map2surface(e2, Q.xi_.col(j), h * Q.w_(j), &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
    }
  }

  BEMBEL_UNUSED_(rot1);
  BEMBEL_UNUSED_(rot2);
  BEMBEL_UNUSED_(ffield_qnodes1);
  BEMBEL_UNUSED_(ffield_qnodes2);
  return;
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif  // BEMBEL_SRC_DUFFYTRICK_INTEGRATE1_HPP_
