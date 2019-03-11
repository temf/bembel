// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "discretization.h"

namespace Bembel {
/**
 *  \brief         Distributes a vector represented in the nodal basis to a
 *                 vector represented in the basis of the big matrix. Result
 *                 will be added to 'ul'.
 *
 *  \param[in]     T              sparse projector matrix
 *  \param[in]     us             Vector to distribute
 *  \param[out]    ul             Result plus input content of ul.
 *  \param[in]     E              Element list
 *  \param[in]     M              Level
 *
 */
void proj_distr_et(discretization *disc, double *us, double *ul) {
  int i;
  int j;
  et_node *E = disc->mesh->E.patch[0];
  sparse *T = disc->T;
  /*
   * get pointer to the leafs of the element tree
   */
  while (E->son[0]) E = E->son[0];

  assert(T->n == disc->na);
  /*
   * apply the projector
   */
  for (i = 0; i < T->m; ++i)
    for (j = 0; j < T->row_number[i]; ++j)
      ul[i] += T->value[i][j] * us[T->index[i][j]];

  return;
}

/**
 *  \brief         Restricts a vector represented in the basis of the big matrix
 *                 to a vector in the nodal basis. Result will be added to 'us'.
 *
 *  \param[in]     T              sparse projector matrix
 *  \param[in]     ul             Vector to restrict
 *  \param[in,out] us             Result plus input content of us.
 *  \param[in]     E              Element list
 *  \param[in]     M              Level
 *
 */
void proj_restr_et(discretization *disc, double *ul, double *us) {
  int i;
  int j;
  et_node *E = disc->mesh->E.patch[0];
  sparse *T = disc->T;
  /*
   * get pointer to the leafs of the element tree
   */
  while (E->son[0]) E = E->son[0];

  assert(T->n == disc->na);
  /*
   * apply the projector
   */
  for (i = 0; i < T->m; ++i)
    for (j = 0; j < T->row_number[i]; ++j)
      us[T->index[i][j]] += T->value[i][j] * ul[i];

  return;
}

void free_discretization(discretization *disc) {
  free_sparse(disc->T);
  free(disc->T);

  return;
}
}  // namespace Bembel