// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "element_tree.h"

namespace Bembel {
/**
 *  \brief   Initializes the hierarchical element tree.
 *
 *  \param[in,out] E        Root of the hierarchical element tree
 *  \param[in]     P        Point list of the single scale basis
 *  \param[in]     F        Element list of the single scale basis
 *  \param[in]     p        Number of patches
 *  \param[in]     M        Level, i.e. 2^M*2^M elements per patch
 *  \param[in]     np       Number of points used for the element list
 *
 *  \return        Number of ansatz functions.
 *
 */
int init_element_tree(et_root *E, vector3 *P, int **F, int p, int M, int np) {
  int i = 0;
  int l = 0;
  int m = 0;
  int n = 0;      /* p*n*n Elemente auf Level m */
  int N = 1 << M; /* p * N * N Panels auf Level M */
  // int nf = p * N * N;
  int na;
  // int *dictv; /* list with index-pairs (point,ansatzf.) */
  // sparsei edges;
  vector3 d1; /* cicrumcenter of element 0 */
  vector3 d2; /* circumcenter of element 2 */
  et_node *pF = NULL;
  double r1; /* cicrumradius of element 1 */
  double r2; /* cicrumradius of element 3 */

  /*
   * initialize memory for each level
   */
  E->nop = p;
  E->patch = (et_node **)calloc(1, p * sizeof(et_node *));
  pF = E->patch[0] = (et_node *)calloc(1, p * sizeof(et_node));
  for (l = 0; l < p; ++l) {
    E->patch[l] = pF + l;
    pF[l].level = 0;
    pF[l].number = l;
    pF[l].patch = l;
    pF[l].index_s = 0;
    pF[l].index_t = 0;
    pF[l].father = NULL;
  }
  n = p;
  for (m = 0; m < M; ++m) {
    pF[0].son[0] = (et_node *)calloc(1, 4 * n * sizeof(et_node));
    for (l = 0; l < n; ++l) {
      /*
       * assigning memory to sons
       */
      pF[l].son[0] = 4 * l + pF[0].son[0];
      pF[l].son[1] = 4 * l + pF[0].son[0] + 1;
      pF[l].son[2] = 4 * l + pF[0].son[0] + 2;
      pF[l].son[3] = 4 * l + pF[0].son[0] + 3;
      /*
       * assigning numbers to sons
       */
      pF[l].son[0]->number = 4 * l;
      pF[l].son[1]->number = 4 * l + 1;
      pF[l].son[2]->number = 4 * l + 2;
      pF[l].son[3]->number = 4 * l + 3;
      /*
       * assigning level to sons
       */
      pF[l].son[0]->level = m + 1;
      pF[l].son[1]->level = m + 1;
      pF[l].son[2]->level = m + 1;
      pF[l].son[3]->level = m + 1;
      /*
       * assigning patch to sons
       */
      pF[l].son[0]->patch = pF[l].patch;
      pF[l].son[1]->patch = pF[l].patch;
      pF[l].son[2]->patch = pF[l].patch;
      pF[l].son[3]->patch = pF[l].patch;
      /*
       * assigning indices to sons
       */
      pF[l].son[0]->index_s = 2 * pF[l].index_s;
      pF[l].son[1]->index_s = 2 * pF[l].index_s + 1;
      pF[l].son[2]->index_s = 2 * pF[l].index_s + 1;
      pF[l].son[3]->index_s = 2 * pF[l].index_s;
      pF[l].son[0]->index_t = 2 * pF[l].index_t;
      pF[l].son[1]->index_t = 2 * pF[l].index_t;
      pF[l].son[2]->index_t = 2 * pF[l].index_t + 1;
      pF[l].son[3]->index_t = 2 * pF[l].index_t + 1;
      /*
       * assigning father to sons
       */
      pF[l].son[0]->father = pF + l;
      pF[l].son[1]->father = pF + l;
      pF[l].son[2]->father = pF + l;
      pF[l].son[3]->father = pF + l;
    }
    pF = pF[0].son[0];
    n *= 4;
  }
  /*
   * assigning sons, vertices, balls to leafs
   */
  for (l = 0; l < p * N * N; ++l) {
    pF[l].son[0] = NULL;
    pF[l].son[1] = NULL;
    pF[l].son[2] = NULL;
    pF[l].son[3] = NULL;
    n = pF[l].patch * N * N + pF[l].index_t * N + pF[l].index_s;
    pF[l].vertex[0] = F[n][0];
    pF[l].vertex[1] = F[n][1];
    pF[l].vertex[2] = F[n][2];
    pF[l].vertex[3] = F[n][3];
    unify(&d1, &r1, P[F[n][0]], 0, P[F[n][2]], 0);
    unify(&d2, &r2, P[F[n][1]], 0, P[F[n][3]], 0);
    unify(&(pF[l].midpoint), &(pF[l].radius), d1, r1, d2, r2);
  }
  /*
   * assigning recursively vertices, balls to father elements
   */
  n = p * N * N;
  for (i = M - 1; i >= 0; --i) {
    pF = pF[0].father;
    n /= 4;
    for (l = 0; l < n; ++l) {
      pF[l].vertex[0] = pF[l].son[0]->vertex[0];
      pF[l].vertex[1] = pF[l].son[1]->vertex[1];
      pF[l].vertex[2] = pF[l].son[2]->vertex[2];
      pF[l].vertex[3] = pF[l].son[3]->vertex[3];
      unify(&d1, &r1, pF[l].son[0]->midpoint, pF[l].son[0]->radius,
            pF[l].son[2]->midpoint, pF[l].son[2]->radius);
      unify(&d2, &r2, pF[l].son[1]->midpoint, pF[l].son[1]->radius,
            pF[l].son[3]->midpoint, pF[l].son[3]->radius);
      unify(&(pF[l].midpoint), &(pF[l].radius), d1, r1, d2, r2);
    }
  }
  return na;
}

/**
 *  \brief   Frees memory of hierarchical element tree.
 *
 *  \param[in,out] E        Root of the hierarchical element tree
 *
 */
int free_element_tree(et_root *E) {
  et_node *pF = E->patch[0];
  et_node *pS = NULL;

  while (pF) {
    pS = pF[0].son[0];
    free(pF);
    pF = pS;
  }

  free(E->patch);

  return 0;
}
}  // namespace Bembel