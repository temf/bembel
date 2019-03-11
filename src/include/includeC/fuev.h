// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_FUEV_
#define __BEMBEL_C_FUEV_

#include "element_tree.h"
#include "gamma.h"

namespace Bembel {
int fuev(et_node *E, double *rhoc, double *rhon, int nf, int M);
}
#endif
