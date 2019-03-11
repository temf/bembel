// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;

double divcheck(double *rho, discretization *disc) {
  meshdata *mesh = disc->mesh;
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int a_bs = disc->a_bs;
  int N = a_bs * nf;
  int zi;            /* row index of the element */
  double h = 1. / n; /* step size */
  double w;
  double value[2];
  double w_div[2];
  double *d_dx;
  double *d_dy;
  double *rhol;
  vector2 s;
  vector2 t;
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1); /* Kubatur-Formeln */

  /*
   * allocate memory
   */
  rhol = (double *)calloc(4 * a_bs * nf, sizeof(double));
  proj_distr_et(disc, rho, rhol);

  /*
   * find first node in element list
   */
  for (int i = 0; i < mesh->M; ++i) pE = pE[0].son[0];

  value[0] = 0.;
  ;
  value[1] = 0.;
  ;
  d_dx = (double *)calloc(disc->a_bs, sizeof(double));
  d_dy = (double *)calloc(disc->a_bs, sizeof(double));
  for (zi = 0; zi < nf; ++zi) {
    /*
     * find place on the unit square
     */
    s.x = pE[zi].index_s * h;
    s.y = pE[zi].index_t * h;

    for (int k = 0; k < Q[g].nop; ++k) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
      /*
       * evaluate divergence of Galerkin solution on unit square, scale
       * with quadrature weight
       */
      memset(d_dx, 0, a_bs * sizeof(double));
      disc->phiphi_dx(d_dx, Q[g].xi[k]);
      w_div[0] = myddot(a_bs, rhol + a_bs * zi, d_dx);
      w_div[1] = myddot(a_bs, rhol + a_bs * zi + N, d_dx);
      memset(d_dy, 0, a_bs * sizeof(double));
      disc->phiphi_dy(d_dy, Q[g].xi[k]);
      w_div[0] += myddot(a_bs, rhol + a_bs * zi + 2 * N, d_dy);
      w_div[1] += myddot(a_bs, rhol + a_bs * zi + 3 * N, d_dy);
      value[0] += w_div[0] * Q[g].w[k];
      value[1] += w_div[1] * Q[g].w[k];
    }
  }

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhol);
  free(d_dx);
  free(d_dy);

  return fabs(value[0] + value[1]);
}

/**
 *  \author        Juergen Doelz, Mods. by Felix Wolf
 */
int Test::test_divconf() {
  int count = 0;

  auto bool_assert = [&count](bool in) {
    if (in) {
      return;
    } else
      count++;
    return;
  };

  int M = 2;
  int iorder = 3;
  double *result;
  double *unit;
  MaxwellSingle maxwell;
  Geometry geom(Bembel::Test::mkSphere());
  Discretization<MaxwellSingle> modern_disc;
  modern_disc.init_Discretization(geom, maxwell, iorder + 1, 1, M);

  discretization disc = modern_disc.get_disc();

  // check for conformity
  unit = (double *)calloc(disc.na, sizeof(double));
  result = (double *)calloc(disc.na, sizeof(double));
  for (int i = 0; i < disc.na; ++i) {
    memset(unit, 0, disc.na * sizeof(double));
    unit[i] = 1.;
    result[i] = divcheck(unit, &disc);
    // std::cout << result[i] << std::endl;
    bool_assert(fabs(result[i]) < 1e-14);
  }

  free(unit);
  free(result);

  return count;
}
