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

#ifndef BEMBEL_SRC_DUFFYTRICK_FARFIELDQUADRATURENODES_HPP_
#define BEMBEL_SRC_DUFFYTRICK_FARFIELDQUADRATURENODES_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *  \brief evaluates a given quadrature on all surface panels storage format
 *         is qNodes.col(k) = [xi, w, Chi(xi); dsChi(xi); dtChi(xi)]\in\Rbb^12
 **/
template <class T>
std::vector<ElementSurfacePoints> computeFfieldQnodes(const T &super_space,
                                                      const Cubature &Q) {
  std::vector<ElementSurfacePoints> ffield_qnodes;
  int next = 0;
  // assume isotropic mesh width h!
  double h = (super_space.get_mesh().get_element_tree().cpbegin())->get_h();
  auto nE = super_space.get_mesh().get_number_of_elements();
  auto pbegin = super_space.get_mesh().get_element_tree().cpbegin();
  auto pend = super_space.get_mesh().get_element_tree().cpend();
  ffield_qnodes.reserve(nE);
  SurfacePoint surfpt;

  for (auto it = pbegin; it != pend; ++it) {
    ffield_qnodes.emplace_back(Q.xi_.cols());
    for (auto k = 0; k < Q.xi_.cols(); ++k) {
      // the quadrature weight is scaled by mesh width
      // this corresponds to a scaling of the basis functions
      // with respect to the L^2 norm!
      super_space.map2surface(*it, Q.xi_.col(k), h * Q.w_(k),
                              &ffield_qnodes[it->id_][k]);
    }
  }
  return ffield_qnodes;
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif  // BEMBEL_SRC_DUFFYTRICK_FARFIELDQUADRATURENODES_HPP_
