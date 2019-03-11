// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_CLUSTER_TREE_
#define __BEMBEL_C_CLUSTER_TREE_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "element_tree.h"
#include "gauss_legendre.h"
#include "hmatrixfactory.h"
#include "int_core.h"
#include "trafos.h"

namespace Bembel {

typedef struct _ct_node ct_node;
typedef struct _ct_root ct_root;
typedef struct _ct_fmat ct_fmat;
typedef struct _ct_ker ct_ker;
typedef struct _ct_rkmat ct_rkmat;

struct _ct_node {
  ct_node *pdad;
  ct_node **pch;
  ct_fmat *fmat;

  ct_rkmat *rkmat;
  int cluster1;
  int cluster2;
  int level;
  int nch;
};

struct _ct_root {
  ct_node *root;
  discretization *disc;
  hmatrixsettings *hmatset;
  char sym;
  int p;
  int M;
  int rank;
  double **Tmom_left;
  double **Tmom_right;
  double ***Ttr;
};

struct _ct_leaf {
  ct_node *father;
};

struct _ct_fmat {
  int bs;
  double *A;
  double fnorm;
};

struct _ct_ker {
  double *l;
  double *r;
  double *A;
  int rk;
};

struct _ct_rkmat {
  ct_ker *ker;
  double **l;
  double **r;
  double fnorm;
  int k;
  int bs;
  int lind1;
  int cind1;
};
}  // namespace Bembel
#include "tree_leafs.h"
/*
 * tree_leafs has to be included here, in order to make structs available
 */

namespace Bembel {

int ct_get_nsons(ct_node *);
int compare_cluster(et_node *, et_node *, int, hmatrixsettings *);
void init_cluster_tree(hmatrixfactory *, ct_root **);
void fast_print_cluster_tree(ct_node *, char *);
int print_cluster_tree(double, double, ct_node *, FILE *, FILE *);
int free_cluster_tree(ct_root *);

}  // namespace Bembel
#endif
