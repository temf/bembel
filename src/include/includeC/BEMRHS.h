// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_BEMRHS_
#define __BEMBEL_C_BEMRHS_

#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "cubature.h"
#include "data.h"
#include "discretization.h"
#include "gamma.h"
#include "gauss_square.h"
#include "mycblas.h"
#include "myvector.h"

namespace Bembel {
void BEMRHS_Dirichlet(double **, discretization *, double (*)(vector3));
void BEMRHS_HelmholtzDirichlet(double **, discretization *);
void BEMRHS_MaxwellDirichletTrace(double **rhs, discretization *disc,
                                  void (*f)(vector3 field[2], vector3 X,
                                            double kappa[2]));

// This one has been replaced by a wrapped version get_rhs
// void BEMRHS_MaxwellDirichletTangentialTrace(double **rhs, discretization
// *disc,
//                                             void (*f)(vector3 field[2],
//                                                       vector3 X,
//                                                       double kappa[2]));

void BEMRHS_TangentialL2(double **rhs, discretization *disc,
                         void (*f)(vector3 field[2], vector3 X,
                                   double kappa[2]));
void BEMRHS_Hdiv(double **rhs, discretization *disc,
                 void (*f)(vector3 field[2], vector3 X, double kappa[2]),
                 void (*fcurl)(vector3 field[2], vector3 X, double kappa[2]));

void BEMRHS_Neumann(double **, discretization *, vector3 (*)(vector3));
}  // namespace Bembel
#endif
