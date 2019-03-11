// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
/*=======================================*
 *  Stuetzstellen Xi und Gewichte G der  *
 *  Gauss-Quadraturformeln auf [0,1].    *
 *=======================================*/

#ifndef __BEMBEL_C_GAUSS_LEGENDRE_
#define __BEMBEL_C_GAUSS_LEGENDRE_

#include <stdio.h>
#include <stdlib.h>
#include "quadrature.h"

namespace Bembel {

void init_Gauss_Legendre(quadrature **Q, int g);
void free_Gauss_Legendre(quadrature **Q, int g);

extern const int g_max;
}  // namespace Bembel
#endif
