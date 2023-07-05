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
#ifndef BEMBEL_SRC_DUFFYTRICK_INTEGRATE3_HPP_
#define BEMBEL_SRC_DUFFYTRICK_INTEGRATE3_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *    \brief quadrature routine for the common edge case
 *    \todo  be sure that map2element computes the weight h*Q.w(i) such that
 *           the integrand may then be scaled by qp1.weight * qp2.weight
 *           here we just set one weight to the actual weight, while the
 *           other one will be set to 1. This is to remain conforming to the
 *           structure of integrate0/1
 *           Information that map2surface has to provide:
 *           xi; w; Chi(xi); dChidx(xi); dChidy(xi);
 **/
template <typename Derived, class T>
void integrate3(const LinearOperatorBase<Derived> &LinOp, const T &super_space,
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
  SurfacePoint qp1, qp2;
  // llc of the element wrt [0,1]^2
  for (auto i = 0; i < Q.w_.size(); ++i) {
    auto xi = Q.xi_.col(i);
    double w = h * h * Q.w_(i) * xi(1) * xi(1);
    t1 = xi(0) * (1 - xi(1));
    t2 = (1 - xi(0)) * (1 - xi(1));
    for (auto j = 0; j < Q.w_.size(); ++j) {
      auto eta = xi(1) * Q.xi_.col(j);
      t3 = xi(0) * (1 - eta(0));
      t4 = (1 - xi(0)) * (1 - eta(0));
      super_space.map2surface(e1, tau(t1, eta(0), rot1), w, &qp1);
      super_space.map2surface(e2, tau(t2, eta(1), rot2), Q.w_(j) * (1 - xi(1)),
                              &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
      super_space.map2surface(e1, tau(1 - t1, eta(0), rot1), w, &qp1);
      super_space.map2surface(e2, tau(1 - t2, eta(1), rot2),
                              Q.w_(j) * (1 - xi(1)), &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
      super_space.map2surface(e1, tau(t3, xi(1), rot1), w, &qp1);
      super_space.map2surface(e2, tau(t4, eta(1), rot2), Q.w_(j) * (1 - eta(0)),
                              &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
      super_space.map2surface(e1, tau(1 - t3, xi(1), rot1), w, &qp1);
      super_space.map2surface(e2, tau(1 - t4, eta(1), rot2),
                              Q.w_(j) * (1 - eta(0)), &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
      super_space.map2surface(e1, tau(t4, eta(1), rot1), w, &qp1);
      super_space.map2surface(e2, tau(t3, xi(1), rot2), Q.w_(j) * (1 - eta(0)),
                              &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
      super_space.map2surface(e1, tau(1 - t4, eta(1), rot1), w, &qp1);
      super_space.map2surface(e2, tau(1 - t3, xi(1), rot2),
                              Q.w_(j) * (1 - eta(0)), &qp2);
      LinOp.evaluateIntegrand(super_space, qp1, qp2, intval);
    }
  }
  BEMBEL_UNUSED_(ffield_qnodes1);
  BEMBEL_UNUSED_(ffield_qnodes2);
  return;
}
}  // namespace DuffyTrick
}  // namespace Bembel

#endif  // BEMBEL_SRC_DUFFYTRICK_INTEGRATE3_HPP_
