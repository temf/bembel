// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_int_core__
#define __BEMBEL_C_int_core__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "constants.h"
#include "cubature.h"
#include "element_tree.h"
#include "gamma.h"
#include "gauss_square.h"
#include "hmatrixfactory.h"
#include "myvector.h"
#include "pdeproblem.h"
#include "trafos.h" /* Aehnlichkeitstransformationen */

namespace Bembel {
void get_stiff(hmatrixfactory *, vector3 *, double *, et_node *, int, int, int);
}
#endif
