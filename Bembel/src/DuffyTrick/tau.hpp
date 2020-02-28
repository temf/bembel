// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_DUFFYTRICK_TAU_H_
#define BEMBEL_DUFFYTRICK_TAU_H_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *  \brief computes rotations for the Duffy trick
 **/
Eigen::Vector2d tau(double x, double y, int thecase) {
  Eigen::Vector2d retval;
  switch (thecase) {
    case 1:
      retval << 1 - y, x;
      return retval;
    case 2:
      retval << 1 - x, 1 - y;
      return retval;
    case 3:
      retval << y, 1 - x;
      return retval;
    default:
      retval << x, y;
      return retval;
  }
}
}  // namespace DuffyTrick
}  // namespace Bembel
#endif
