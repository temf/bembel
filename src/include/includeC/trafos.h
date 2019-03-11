// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_TRAFOS__
#define __BEMBEL_C_TRAFOS__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "myvector.h"

namespace Bembel {
vector2 Kappa(vector2 a, vector2 b, double h);
vector2 Tau(double b_x, double b_y, int CASE);
int compare(vector3 *P, int *F1, int *F2, int *ind1, int *ind2);
}  // namespace Bembel
#endif
