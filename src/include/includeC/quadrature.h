// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_QUADRATURE_
#define __BEMBEL_C_QUADRATURE_

namespace Bembel {
/**
 *  \brief         Defines the quadrature formulas on the interval [0,1].
 */
typedef struct {
  int nop; /**< number of integration points and     */
  /*
   * weights
   */

  double *xi; /**< integration points                   */

  double *w; /**< integration weights                  */
} quadrature;
}  // namespace Bembel
#endif
