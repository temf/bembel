// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "gram.h"

namespace Bembel {
static const int maxiter = 600; /* for GMRES */

/**
 *  \brief      Matrix-vector multiplication V*x for sparse matrices
 *
 *  \param[in]     G              Sparse matrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *  \param[in]     p              Number of patches of H.
 *  \param[in]     M              Number of levels.
 *  \param[in]     E              Element list.
 */
int gram_timesV(sparse *G, double *x, double *y, discretization *disc) {
  int i;
  int j;
  double *lx; /* source big vector */
  double *ly; /* target big vector */

  /*
   * get temporary memory
   */
  lx = (double *)calloc(G->n, sizeof(double));
  ly = (double *)calloc(G->m, sizeof(double));

  /*
   * distribute vectors
   */
  proj_distr_et(disc, x, lx);

  /*
   * apply big matrix to big vectors
   */
  for (i = 0; i < G->m; ++i)
    for (j = 0; j < G->row_number[i]; ++j)
      ly[i] += G->value[i][j] * lx[G->index[i][j]];

  /*
   * restrict destination vector
   */
  proj_restr_et(disc, ly, y);

  free(lx);
  free(ly);

  return 0;
}

/**
 *  \brief         Preconditioned conjugate gradient method for sparse matrices.
 *
 *  Preconditioned conjugate gradient method using a Jacobi preconditioner to
 *  solve A*x = b.
 *
 *  \param[in]     G              Gram matrix.
 *  \param[in]     b              Right hand side.
 *  \param[in,out] x              Starting value and result.
 *  \param[in]     epsi           Tolerance for the solution.
 *  \param[in]     na             Length of the vectors x and b.
 *  \param[in]     E              Element list.
 *  \param[in]     M              Discretization level.
 *
 */
void graminv_timesV(sparse *G, double *b, double *x, double epsi,
                    discretization *disc) {
  int i;
  int k;
  int na = disc->na;
  double P;
  double Q;
  double omg;
  double *r;
  double *d;
  double *Ad;
  double *D;
  double *DH;

  /*
   * Initialisierung
   */
  D = (double *)malloc(na * sizeof(double));
  DH = (double *)malloc(disc->T->m * sizeof(double));
  r = (double *)malloc(na * sizeof(double));
  d = (double *)malloc(na * sizeof(double));
  Ad = (double *)malloc(na * sizeof(double));

  /*
   * Preconditioning: D ~= diag(A)^(-1), but not exactly
   */
  for (i = 0; i < disc->T->m; ++i) DH[i] = get_sparse(G, i, i);
  memset(D, 0, na * sizeof(double));
  proj_restr_et(disc, DH, D);
  for (i = 0; i < na; ++i)
    if (fabs(D[i]) < 1e-8)
      D[i] = 1.;
    else
      D[i] = 1. / D[i];

  /*
   * r = b - A*x
   */
  memset(r, 0, na * sizeof(double));
  gram_timesV(G, x, r, disc);
  for (i = 0; i < na; ++i) r[i] = b[i] - r[i];

  /*
   * d = D*r und Q = (r,D*r)
   */
  Q = 0;
  for (i = 0; i < na; i++) {
    d[i] = D[i] * r[i];
    Q += d[i] * r[i];
  }

  /*
   * Iteration
   */
  for (k = 0; sqrt(Q) > epsi; k++) {
    /*
     * Ad = A*d
     */
    memset(Ad, 0, na * sizeof(double));
    gram_timesV(G, d, Ad, disc);

    /*
     * omg = Q / (d,Ad)
     */
    omg = 0;
    for (i = 0; i < na; i++) omg += d[i] * Ad[i];
    omg = Q / omg;
    P = Q;

    /*
     * x = x + omg * d, r = r - omg * Ad und Q = (r,D*r)
     */
    Q = 0;
    for (i = 0; i < na; i++) {
      x[i] += omg * d[i];
      r[i] -= omg * Ad[i];
      Q += r[i] * D[i] * r[i];
    }

    /*
     * d = D * r + Q/P * d
     */
    omg = Q / P;
    for (i = 0; i < na; i++) d[i] = D[i] * r[i] + omg * d[i];
  }

  printf("                           %d iterations \n", k);

  free(D);
  free(DH);
  free(r);
  free(d);
  free(Ad);

  return;
}
}  // namespace Bembel