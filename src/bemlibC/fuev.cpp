// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "fuev.h"
namespace Bembel {
/**
 *  \brief         Evaluates the density function at the vertices of the
 *                 elements for the canonical scalar product and adds the
 *                 values to the values in rhon.
 *
 *  \param[in]     E              Element list
 *  \param[in]     rhoc           Coefficients of rho
 *  \param[in,out] rhon           Values of the function evaluation
 *  \param[in]     nf             Length of E
 *  \param[in]     M              Number of levels
 *
 */
int fuev(et_node *E, double *rhoc, double *rhon, int nf, int M) {
  int i;
  int j;
  int n = 1 << M;
  int np;

  /*
   * get pointer to the leafs of the element tree
   */
  while (E->son[0]) E = E->son[0];

  /*
   * get np
   */
  np = 0;
  for (i = 0; i < nf; ++i)
    for (j = 0; j < 4; ++j)
      if (E[i].vertex[j] > np) np = E[i].vertex[j];
  ++np;

  /*
   * scale nodal representation
   */
  for (i = 0; i < np; ++i) rhon[i] = n * rhoc[i];

  return 0;
}
}  // namespace Bembel