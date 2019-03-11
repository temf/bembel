// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_GAUSS_SQUARE_
#define __BEMBEL_C_GAUSS_SQUARE_
/********************
 *  Gauss_Square.h  *
 ********************/
#include "cubature.h"
/*============================================*
 *  Definiert die Gauss-Quadraturformeln auf  *
 *  dem Referenzviereck [0,1]^2.              *
 *============================================*/
namespace Bembel {
void init_Gauss_Square(cubature **Q, int g);

void free_Gauss_Square(cubature **Q, int g);
}  // namespace Bembel
#endif
