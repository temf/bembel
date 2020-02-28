// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_DUFFYTRICK_FARFIELDQUADRATURENODES_H_
#define BEMBEL_DUFFYTRICK_FARFIELDQUADRATURENODES_H_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *  \brief evaluates a given quadrature on all surface panels storage format
 *         is qNodes.col(k) = [xi, w, Chi(xi); dsChi(xi); dtChi(xi)]\in\Rbb^12
 **/
template <class T>
Eigen::Matrix<double, 12, Eigen::Dynamic> computeFfieldQnodes(
    const T &super_space, const Cubature &Q) {
  Eigen::Matrix<double, 12, Eigen::Dynamic> ffield_qnodes;
  int next = 0;
  // assume isotropic mesh width h!
  double h = (super_space.get_mesh().get_element_tree().cpbegin())->get_h();
  auto nE = super_space.get_mesh().get_number_of_elements();
  auto pbegin = super_space.get_mesh().get_element_tree().cpbegin();
  auto pend = super_space.get_mesh().get_element_tree().cpend();
  ffield_qnodes.resize(12, nE * Q.xi_.cols());
  SurfacePoint surfpt;

  for (auto it = pbegin; it != pend; ++it)
    for (auto k = 0; k < Q.xi_.cols(); ++k) {
      // the quadrature weight is scaled by mesh width
      // this corresponds to a scaling of the basis functions
      // with respect to the L^2 norm!
      super_space.map2surface(*it, Q.xi_.col(k), h * Q.w_(k), &surfpt);
      ffield_qnodes.col(it->id_ * Q.xi_.cols() + k) = surfpt;
    }
  return ffield_qnodes;
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif
