// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_init_gridpoints__
#define __BEMBEL_C_init_gridpoints__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "gamma.h"
#include "myvector.h"

namespace Bembel {
typedef struct _gridcell gridcell;
struct _gridcell {
  int c1, c2, c3, c4, c5, c6, c7, c8;
};

int init_gridpoints(vector3 **Q, int *nq, gridcell **Cells, int *nc, double h,
                    double h_disc, const geometry &Chi);
int print_cells(vector3 *P, int nq, double *Pot, gridcell *Cells, int nc,
                char *dname);
}  // namespace Bembel
#endif
