// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_H_print_geometry__
#define __BEMBEL_C_H_print_geometry__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#include "constants.h"
#include "data.h"
#include "discretization.h"
#include "element_tree.h"
#include "fuev.h"
#include "meshdata.h"
#include "myvector.h"
namespace Bembel {
int print_geometry_Maxwell(discretization *disc, double *rho, char *dname,
                           double kappa[2]);
int print_potential_Maxwell(double *Pot, vector3 *Q, int nq, char *dname,
                            double kappa[2]);

int print_geometry_level(et_node *E, vector3 *P, int nf, double *rho, int l,
                         char *dname);
}  // namespace Bembel
#endif
