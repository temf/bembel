// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_H2TRANSFORM_
#define __BEMBEL_C_H2TRANSFORM_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cluster_tree.h"
#include "constants.h"
#include "mycblas.h"

namespace Bembel {

typedef struct _H2tree H2tree;

struct _H2tree {
  int depth;
  int nsons;
  H2tree *sons;
  double *data;
};

void new_H2tree(H2tree *t, int depth, int nsons, int np2);
void free_H2tree(H2tree *t);
void add_H2tree(H2tree *t, H2tree *a, int np2);

int ForwardTransform(ct_root *, double *x, H2tree *transx);

int BackwardTransform(ct_root *, double *x, H2tree *transx);
}  // namespace Bembel
#endif
