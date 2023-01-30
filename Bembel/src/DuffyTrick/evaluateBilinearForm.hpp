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
#ifndef BEMBEL_SRC_DUFFYTRICK_EVALUATEBILINEARFORM_HPP_
#define BEMBEL_SRC_DUFFYTRICK_EVALUATEBILINEARFORM_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *  \brief  This function wraps the quadrature routines for the duffy trick
 *          and returns all integrals for the given pair of elements
 */
template <typename Derived, class T, class CubatureVector>
void evaluateBilinearForm(
    const LinearOperatorBase<Derived>& linOp, const T& super_space,
    const ElementTreeNode& e1, const ElementTreeNode& e2,
    const CubatureVector& GS, const Eigen::MatrixXd& ffield_qnodes1,
    const Eigen::MatrixXd& ffield_qnodes2,
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, Eigen::Dynamic>* intval) {
  //////////////////////////////////////////////////////////////////////////////
  double dist = 0;
  int ffield_deg =
      linOp.get_FarfieldQuadratureDegree(super_space.get_polynomial_degree());
  int nfield_deg = 0;
  auto cp = compareElements(e1, e2, &dist);
  nfield_deg = linOp.getNearfieldQuadratureDegree(
      super_space.get_polynomial_degree(), dist, e1.level_);
  // make sure that the quadratur degree is at least the far field degree
  nfield_deg = nfield_deg >= ffield_deg ? nfield_deg : ffield_deg;
  assert(nfield_deg < Constants::maximum_quadrature_degree &&
         "nfield_deg too large, increase maximum_quadrature_degree");
  auto Q = GS[nfield_deg];
  switch (cp(2)) {
    case 0:
      if (nfield_deg == ffield_deg) {
        integrate0(linOp, super_space, e1, 0, e2, 0, ffield_qnodes1,
                   ffield_qnodes2, Q, intval);
        return;
      } else {
        integrate1(linOp, super_space, e1, 0, e2, 0, ffield_qnodes1,
                   ffield_qnodes2, Q, intval);
        return;
      }
    case 1:
      assert(!"you should not have ended up here!");
    case 2:
      integrate2(linOp, super_space, e1, 0, e2, 0, ffield_qnodes1,
                 ffield_qnodes2, Q, intval);
      return;
    case 3:
      integrate3(linOp, super_space, e1, cp(0), e2, cp(1), ffield_qnodes1,
                 ffield_qnodes2, Q, intval);
      return;
    case 4:
      integrate4(linOp, super_space, e1, cp(0), e2, cp(1), ffield_qnodes1,
                 ffield_qnodes2, Q, intval);
      return;
    default:
      assert(!"you should not have ended up here!");
  }
  return;
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif  // BEMBEL_SRC_DUFFYTRICK_EVALUATEBILINEARFORM_HPP_
