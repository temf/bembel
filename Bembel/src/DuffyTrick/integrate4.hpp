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
#ifndef BEMBEL_SRC_DUFFYTRICK_INTEGRATE4_HPP_
#define BEMBEL_SRC_DUFFYTRICK_INTEGRATE4_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *    \brief quadrature routine for identical elements
 *    \todo  be sure that map2element computes the weight h*Q.w(i) such that
 *           the integrand may then be scaled by qp1.weight * qp2.weight
 *           here we just set one weight to the actual weight, while the
 *other one will be set to 1. This is to remain conforming to the structure
 *           of integrate0/1
 *           Information that map2element has to provide:
 *           xi; w; Chi(xi); dChidx(xi); dChidy(xi);
 **/
template <typename Derived, class T>
void integrate4(const LinearOperatorBase<Derived> &LinOp, const T &super_space,
                const ElementTreeNode &e1, int rot1, const ElementTreeNode &e2,
                int rot2, const ElementSurfacePoints &ffield_qnodes1,
                const ElementSurfacePoints &ffield_qnodes2, const Cubature &Q,
                Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                              Eigen::Dynamic, Eigen::Dynamic> *intval) {
  intval->setZero();
  double h = e1.get_h();
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  double t4 = 0;
  SurfacePoint qp1, qp2, qp3, qp4, qp5, qp6;
  Eigen::Vector2d pt1;
  // llc of the element wrt [0,1]^2
  for (auto i = 0; i < Q.w_.size(); ++i) {
    Eigen::Vector2d xi = Q.xi_.col(i);
    double w = h * h * Q.w_(i) * std::pow(xi(0), 3.);
    xi(1) *= xi(0);
    super_space.map2surface(e1, tau(xi(0), xi(1), rot1), w, &qp1);
    super_space.map2surface(e1, tau(xi(1), xi(0), rot1), w, &qp2);
    super_space.map2surface(e2, tau(xi(0), xi(1), rot2), w, &qp3);
    super_space.map2surface(e2, tau(xi(1), xi(0), rot2), w, &qp4);
    for (auto j = 0; j < Q.w_.size(); ++j) {
      auto eta = xi(0) * Q.xi_.col(j);
      super_space.map2surface(e2, tau(eta(0), eta(1), rot2), Q.w_(j), &qp5);
      super_space.map2surface(e1, tau(eta(0), eta(1), rot1), Q.w_(j), &qp6);
      LinOp.evaluateIntegrand(super_space, qp1, qp5, intval);
      LinOp.evaluateIntegrand(super_space, qp2, qp5, intval);
      LinOp.evaluateIntegrand(super_space, qp6, qp3, intval);
      LinOp.evaluateIntegrand(super_space, qp6, qp4, intval);
    }
  }
  BEMBEL_UNUSED_(ffield_qnodes1);
  BEMBEL_UNUSED_(ffield_qnodes2);
  return;
}
}  // namespace DuffyTrick
}  // namespace Bembel

#endif  // BEMBEL_SRC_DUFFYTRICK_INTEGRATE4_HPP_
