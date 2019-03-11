// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_hmatrixfactory__
#define __BEMBEL_C_hmatrixfactory__

#include "discretization.h"
#include "hmatrixsettings.h"
namespace Bembel {
typedef struct {
  discretization *disc;
  hmatrixsettings *hmatset;

  /*
   * temporary variables initialized and destroyed during matrix assembly
   */
  randwerte **RW;
  cubature *Q;

  /* These variables are only for statistical reasons */
  int assemfmats;
  int assemsfmats;
  int assemrkmats;

} hmatrixfactory;
}  // namespace Bembel
#endif
