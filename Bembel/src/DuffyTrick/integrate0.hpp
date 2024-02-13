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
#ifndef BEMBEL_SRC_DUFFYTRICK_INTEGRATE0_HPP_
#define BEMBEL_SRC_DUFFYTRICK_INTEGRATE0_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *    \brief far-field quadrature routine, which is based on precomputed values
 *           in order to quickly evaluate the integrand in the case that the
 *           far-field quadrature degree can be used
 **/
template <typename Derived, class T>
void integrate0(const LinearOperatorBase<Derived> &LinOp, const T &super_space,
                const ElementTreeNode &e1, int rot1, const ElementTreeNode &e2,
                int rot2, const ElementSurfacePoints &ffield_qnodes1,
                const ElementSurfacePoints &ffield_qnodes2, const Cubature &Q,
                Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                              Eigen::Dynamic, Eigen::Dynamic> *intval) {
  intval->setZero();
  for (auto i = 0; i < Q.w_.size(); ++i)
    for (auto j = 0; j < Q.w_.size(); ++j)
      LinOp.evaluateIntegrand(super_space, ffield_qnodes1[i], ffield_qnodes2[j],
                              intval);
  BEMBEL_UNUSED_(rot1);
  BEMBEL_UNUSED_(rot2);
  BEMBEL_UNUSED_(Q);
  return;
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif  // BEMBEL_SRC_DUFFYTRICK_INTEGRATE0_HPP_
