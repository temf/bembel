// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_MULTIPOLE__
#define __BEMBEL_C_MULTIPOLE__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cluster_tree.h"
#include "constants.h"
#include "element_tree.h"
#include "gamma.h"
#include "gauss_legendre.h"
#include "hmatrixsettings.h"

namespace Bembel {

int cheb_roots(double *, int);
int interpol_kernel(hmatrixfactory *, double ***, double ***, et_node *,
                    et_node *, int);
int init_polynomials(double ***, int);
double eval_polynomial(double *, int, double);
int init_moments(double ***, int, int, int, discretization *);
int assemble_moments(int M, int lvl, double ***Tmom, double ****Ttr, int np,
                     discretization *);
int init_momentmatrices(discretization *, hmatrixsettings *, double ***,
                        void (**)(double *, double, double));
int init_transfermatrices(discretization *, hmatrixsettings *, double ****);
int free_momentmatrices(discretization *, hmatrixsettings *, double **);
int free_transfermatrices(discretization *, hmatrixsettings *, double ***);

}  // namespace Bembel
#endif
