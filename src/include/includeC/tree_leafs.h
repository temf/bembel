// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_TREE_LEAFS_
#define __BEMBEL_C_TREE_LEAFS_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cluster_tree.h"
#include "constants.h"
#include "element_tree.h"
#include "int_core.h"
#include "multipoleH2.h"
#include "myvector.h"

namespace Bembel {

int init_rkmat(hmatrixfactory *, vector3 *, ct_rkmat *[], ct_rkmat *[],
               et_node *, et_node *, int);

int init_fmat(hmatrixfactory *, vector3 *, ct_fmat *[], ct_fmat *[], et_node *,
              et_node *, int);

int init_sym_fmat(hmatrixfactory *, vector3 *, ct_fmat *[], et_node *,
                  et_node *, int);

int free_fmat(ct_fmat *);

int free_rkmat(ct_rkmat *);

int init_parametrix(int);
int free_parametrix(int);

}  // namespace Bembel
#endif
