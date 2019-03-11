// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "pdeproblem.h"

namespace Bembel {
void pdeproblem_print(pdeproblem *pde) {
  printf("\n");
  printf("PDE is %s\n", pde->name.c_str());
  printf("%d H-matrices of type %s\n", pde->nct, pde->symflags.c_str());
  printf("Wavenumber is %g", pde->kappa[0]);
  if (pde->kappa[1]) printf("+%gi\n", pde->kappa[1]);
  printf("\n");
}
}  // namespace Bembel
