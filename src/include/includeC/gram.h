// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_GRAM_
#define __BEMBEL_C_GRAM_
#include <math.h>
#include <string.h>
#include "cluster_tree.h"
#include "constants.h"
#include "cubature.h"
#include "discretization.h"
#include "element_tree.h"
#include "gamma.h"
#include "gauss_square.h"
#include "myvector.h"
#include "sparse.h"
namespace Bembel {
void gram(sparse *A, discretization *disc);

void gramTangentialL2(sparse *A, discretization *disc);
void gramHdiv(sparse *A, discretization *disc);

int gram_timesV(sparse *G, double *x, double *y, discretization *disc);

void graminv_timesV(sparse *G, double *b, double *x, double epsi,
                    discretization *disc);
}  // namespace Bembel
#endif
