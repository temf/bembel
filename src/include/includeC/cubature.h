// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_CUBATURE
#define __BEMBEL_C_CUBATURE

#include "myvector.h"
namespace Bembel {
/**
 *  \brief         Defines the quadrature formulas on the square.
 */
typedef struct {
  int nop; /**< number of integration points and     */
  /*
   * weights
   */

  Bembel::vector2 *xi; /**< integration points                   */

  double *w; /**< integration weights                  */
} cubature;
}  // namespace Bembel
#endif
