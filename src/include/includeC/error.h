// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_ERROR_
#define __BEMBEL_C_ERROR_

#include <math.h>
#include <stdio.h>
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

namespace Bembel {
double error(double *Pot, vector3 *Q, int nq);
double errorcomplex(double *, vector3 *Q, int nq, double kappa[2]);
double errorMaxwell(double *, vector3 *Q, int nq, double kappa[2],
                    void (*f)(vector3 d[2], vector3 a, double kappa[2]));
double errorTangentialL2TangentialTrace(double *rho, discretization *disc,
                                        pdeproblem *pde,
                                        void (*mf)(vector3[2], vector3,
                                                   double[2]));
double errorTangentialHdiv(double *rho, discretization *disc, pdeproblem *pde,
                           void (*H)(vector3[2], vector3, double[2]),
                           void (*E)(vector3[2], vector3, double[2]));
}  // namespace Bembel

#endif
