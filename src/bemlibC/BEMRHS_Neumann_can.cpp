// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "BEMRHS.h"

namespace Bembel {
namespace Rhs {
/**
 *  \brief         Initializes the right hand side of the Galerkin method for
 *                 the canonical scalar product.
 *
 *  \param[in,out] rhs            Place where to store the rhs.
 *  \param[in]     discretization Discretization properties.
 *  \param[in]     mf             Function for right-hand side.
 *
 */
void BEMRHS_Neumann(double **rhs, discretization *disc,
                    vector3 (*mf)(vector3)) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double c;                   /* output of quadrature routines */
  double h = 1. / n;          /* step size */
  double *d;
  double *rhsl;
  vector2 s;   /* left lower corner on square for zi */
  vector2 t;   /* quadrature point on square */
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */
  /*
   * lowest level of the element tree
   */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhsl = (double *)calloc(nf * disc->a_bs, sizeof(double));

  /*
   * find first node in element list
   */
  for (i = 0; i < mesh->M; ++i) pE = pE[0].son[0];

  d = (double *)calloc(disc->a_bs, sizeof(double));
  for (zi = 0; zi < nf; ++zi) {
    /*
     * find place on the unit square
     */
    s.x = pE[zi].index_s * h;
    s.y = pE[zi].index_t * h;

    for (k = 0; k < Q[g].nop; ++k) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
      c = h * Q[g].w[k] *
          vector3_skalp(Chi[pE[zi].patch].n_f(t), mf(Chi[pE[zi].patch].f(t)));

      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c, d, rhsl + disc->a_bs * zi);
    }
  }

  (*rhs) = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhsl, *rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhsl);
  free(d);
  return;
}
}  // namespace Rhs
}  // namespace Bembel