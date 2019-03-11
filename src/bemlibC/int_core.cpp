// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include "int_core.h"

namespace Bembel {
/**
 * \brief computes for the index pair (zi,si) the entry A(zi,si) of the
 *        stiffness matrix
 *
 *  \param[in]     P        Point list
 *  \param[in]     pE       pointer to the panels of the element tree
 *  \param[in]     M        discretization level, necessary to determine
 *                          precision
 *  \param[in]     zi       row index of stiffness matrix entry
 *  \param[in]     si       column index of stiffness matrix entry
 *
 */
void get_stiff(hmatrixfactory *hmatfac, vector3 *P, double *c, et_node *pE,
               int M, int zi, int si) {
  int n = 1 << M; /* number of elements per patch */
  int i1 = 0;     /* index of the first patch */
  int i2 = 0;     /* index of the element zi in Harbrecht order */
  int j1 = 0;     /* index of the second patch */
  int j2 = 0;     /* index of the element si in Harbrecht order */
  int CASE = 0;   /* distinction of cases for the integration */
  int k = 0;      /* Indices for the rotations in the Duffy
                   * trick */
  int l = 0;
  int g = 0;         /* degree of cubature */
  double h = 1. / n; /* stepsize on the reference domain */
  double dx = 0;     /* x/y/z-Distance of the zi and si element */
  double dy = 0;
  double dz = 0;
  double dist = 0; /* total distance of s and t element */
  vector2 s;       /* lower left corners of the element zi */
  vector2 t;       /* lower left corners of the element zi */
  pdeproblem *pde = hmatfac->disc->pde;
  const geometry &geom = *hmatfac->disc->mesh->geom;
  discretization *disc = hmatfac->disc;
  int g_far = disc->g_far;

  /*
   * Initialization
   */
  s.x = pE[zi].index_s * h;
  s.y = pE[zi].index_t * h;
  t.x = pE[si].index_s * h;
  t.y = pE[si].index_t * h;
  i1 = pE[zi].patch;
  j1 = pE[si].patch;
  i2 = pE[zi].patch * n * n + pE[zi].index_t * n + pE[zi].index_s;
  j2 = pE[si].patch * n * n + pE[si].index_t * n + pE[si].index_s;

  /*
   * Computation of cubature degree
   */
  dx = pE[zi].midpoint.x - pE[si].midpoint.x;
  dy = pE[zi].midpoint.y - pE[si].midpoint.y;
  dz = pE[zi].midpoint.z - pE[si].midpoint.z;
  dist = sqrt(dx * dx + dy * dy + dz * dz) - pE[zi].radius - pE[si].radius;
  /*
   * distance relative to element size
   */

  g = pde->quadrature_order(dist, disc->quadrature_accuracy, M, disc->a_o);
  if (g < g_far) g = g_far;

  /*
   * distinction of cases for the integration CASE = 0 -> far-field (g =
   * g_far) CASE = 1 -> nothing to do, but not far-field CASE = 2 -> identical
   * elements CASE = 3 -> common edge CASE = 4 -> common vertex
   */

  /*
   * determine topoligical situation of the two elements
   */
  if (dist > 1e-6) /* distinct elements */
    CASE = 0;
  else if (zi == si) /* identical elements */
    CASE = 2;
  else /* are there common vertices / edges ? */
    CASE = compare(P, pE[zi].vertex, pE[si].vertex, &k, &l);

  /*
   * check if precomputed far field cubature is applicable
   */
  if ((CASE == 0) && (g > g_far)) CASE = 1;

  /*
   * determine cubature routine according to CASE
   */
  switch (CASE) {
    case 0:
      pde->IntPhi0(c, i2, j2, &hmatfac->Q[g_far], hmatfac->RW, disc);
      break;
    case 1:
      pde->IntPhi1(c, s, t, h, geom[i1], geom[j1], &hmatfac->Q[g], disc);
      break;
    case 2:
      pde->IntPhi2(c, s, h, geom[i1], &hmatfac->Q[g], disc);
      break;
    case 3:
      pde->IntPhi3(c, s, t, h, k, l, geom[i1], geom[j1], &hmatfac->Q[g], disc);
      break;
    case 4:
      pde->IntPhi4(c, s, t, h, k, l, geom[i1], geom[j1], &hmatfac->Q[g], disc);
      break;
    default:
      assert(CASE >= 0 && CASE <= 4);
  }

  return;
}
}  // namespace Bembel