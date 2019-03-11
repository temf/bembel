// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_POT__
#define __BEMBEL_C_POT__

#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "cubature.h"
#include "data.h"
#include "discretization.h"
#include "element_tree.h"
#include "gamma.h"
#include "gauss_square.h"
#include "mycblas.h"
#include "myvector.h"
#include "pdeproblem.h"

namespace Bembel {
void pot(double *, double **, vector3 *, int, discretization *, pdeproblem *);
void potcomplex(double *, double **, vector3 *, int, discretization *,
                pdeproblem *);

void potMaxwell(double *rho, double **Pot, vector3 *R, int nr,
                discretization *disc, pdeproblem *pde);

std::vector<double> maxwellTpPotential(std::vector<std::vector<double>> &tps,
                                       double *rho, discretization *disc,
                                       pdeproblem *pde);
}  // namespace Bembel
#endif
