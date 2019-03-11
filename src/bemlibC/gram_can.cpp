// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "gram.h"

namespace Bembel {
/**
 * \brief          Assembles mass matrix
 *
 *  \param[in,out] G              Sparse matrix
 *  \param[in]     discretization Discretization struct
 *  \param[in]     mf             Function which shall be multiplied
 *
 */
void gram(sparse *G, discretization *disc) {
  int i;                                  /* increment variable */
  int k;                                  /* increment variable */
  int g = disc->g_far;                    /* quadrature degree */
  const int p = disc->mesh->geom->size(); /* dumber of parameter domains */
  int M = disc->mesh->M;
  int n = 1 << M;     /* n*n patches per parameter domain */
  int nf = p * n * n; /* number of elements */
  int zi;             /* row index of the element */
  int a_bs = disc->a_bs;
  int a_bs2 = disc->a_bs2;
  double d;
  double h = 1. / n; /* step size */
  double *c;         /* output of quadrature routines */
  vector2 s;         /* left lower corner on square for zi */
  vector2 t;         /* quadrature point on square */
  cubature *Q;       /* quadrature formulas */
  const geometry &Chi = *disc->mesh->geom;
  et_node *pE =
      disc->mesh->E.patch[0]; /* pointer to the first element on the */

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);
  init_sparse(G, a_bs * nf, a_bs * nf, a_bs);

  /*
   * find first node in element list
   */
  for (i = 0; i < M; ++i) pE = pE[0].son[0];

  c = (double *)calloc(a_bs2, sizeof(double));

  for (zi = 0; zi < nf; ++zi) {
    /*
     * find place on the unit square
     */
    s.x = pE[zi].index_s * h;
    s.y = pE[zi].index_t * h;

    memset(c, 0, a_bs2 * sizeof(double));
    for (k = 0; k < Q[g].nop; ++k) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
      d = Q[g].w[k] * vector3_norm(Chi[pE[zi].patch].n_f(t));
      disc->Phi_times_Phi(c, d, Q[g].xi[k], Q[g].xi[k]);
    }

    for (k = 0; k < a_bs2; ++k)
      add_sparse(G, a_bs * zi + k % a_bs, a_bs * zi + k / a_bs,
                 c[k]); /* L^2-normiert!
                         */
  }

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(c);

  return;
}

/**
 * \brief          Assembles mass matrix w.r.t. L^2 innerproduct of tangential
 * space
 *
 *  \param[in,out] G              Sparse matrix
 *  \param[in]     discretization Discretization struct
 *  \param[in]     mf             Function which shall be multiplied
 *
 */
void gramTangentialL2(sparse *A, discretization *disc) {
  int i;                                  /* increment variable */
  int k;                                  /* increment variable */
  int g = disc->g_far;                    /* quadrature degree */
  const int p = disc->mesh->geom->size(); /* dumber of parameter domains */
  int M = disc->mesh->M;
  int n = 1 << M;     /* n*n patches per parameter domain */
  int nf = p * n * n; /* number of elements */
  int zi;             /* row index of the element */
  int a_bs = disc->a_bs;
  int a_bs2 = disc->a_bs2;
  double h = 1. / n;
  double *c;    /* output of quadrature routines */
  vector2 s;    /* quadrature point on square */
  vector2 t;    /* quadrature point on square */
  vector3 t_dx; /* quadrature point on square */
  vector3 t_dy; /* quadrature point on square */
  cubature *Q;  /* quadrature formulas */
  // parametrix *Chi;                      /* parametrizations */
  et_node *pE = disc->mesh->E.patch[0]; /* pointer to the first
                                         * element on the */

  const geometry &Chi = *disc->mesh->geom;

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);
  init_sparse(A, 4 * a_bs * nf, 4 * a_bs * nf, 2 * a_bs);

  /*
   * find first node in element list
   */
  for (i = 0; i < M; ++i) pE = pE[0].son[0];

  c = (double *)calloc(4 * a_bs2, sizeof(double));

  for (zi = 0; zi < nf; ++zi) {
    /*
     * find place on the unit square
     */
    s.x = pE[zi].index_s * h;
    s.y = pE[zi].index_t * h;

    memset(c, 0, 4 * a_bs2 * sizeof(double));

    for (k = 0; k < Q[g].nop; ++k) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
      t_dx = Chi[pE[zi].patch].df_dx(t);
      t_dy = Chi[pE[zi].patch].df_dy(t);
      disc->VPhi_scal_VPhi(
          c, Q[g].w[k] * h / vector3_norm(vector3_mul(t_dx, t_dy)), Q[g].xi[k],
          Q[g].xi[k], t_dx, t_dy, t_dx, t_dy);
    }

    for (k = 0; k < a_bs2; ++k) {
      // real
      add_sparse(A, a_bs * zi + k / a_bs, a_bs * zi + k % a_bs, c[k]);
      add_sparse(A, a_bs * zi + k / a_bs, a_bs * (zi + 2 * nf) + k % a_bs,
                 c[k + a_bs2]);
      add_sparse(A, a_bs * (zi + 2 * nf) + k / a_bs, a_bs * zi + k % a_bs,
                 c[k + 2 * a_bs2]);
      add_sparse(A, a_bs * (zi + 2 * nf) + k / a_bs,
                 a_bs * (zi + 2 * nf) + k % a_bs, c[k + 3 * a_bs2]);
      // imaginary
      add_sparse(A, a_bs * (zi + nf) + k / a_bs, a_bs * (zi + nf) + k % a_bs,
                 c[k]);
      add_sparse(A, a_bs * (zi + nf) + k / a_bs,
                 a_bs * (zi + 3 * nf) + k % a_bs, c[k + a_bs2]);
      add_sparse(A, a_bs * (zi + 3 * nf) + k / a_bs,
                 a_bs * (zi + nf) + k % a_bs, c[k + 2 * a_bs2]);
      add_sparse(A, a_bs * (zi + 3 * nf) + k / a_bs,
                 a_bs * (zi + 3 * nf) + k % a_bs, c[k + 3 * a_bs2]);
    }
  }

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(c);

  return;
}

/**
 * \brief          Assembles mass matrix w.r.t. Hdiv innerproduct of tangential
 *                 space.
 *
 *  \param[in,out] G              Sparse matrix
 *  \param[in]     discretization Discretization struct
 *  \param[in]     mf             Function which shall be multiplied
 *
 */
void gramHdiv(sparse *A, discretization *disc) {
  int i;                                  /* increment variable */
  int k;                                  /* increment variable */
  int g = disc->g_far;                    /* quadrature degree */
  const int p = disc->mesh->geom->size(); /* dumber of parameter domains */
  int M = disc->mesh->M;
  int n = 1 << M;     /* n*n patches per parameter domain */
  int nf = p * n * n; /* number of elements */
  int zi;             /* row index of the element */
  int a_bs = disc->a_bs;
  int a_bs2 = disc->a_bs2;
  double normal;
  double h = 1. / n;
  double *c;    /* output of quadrature routines */
  vector2 s;    /* quadrature point on square */
  vector2 t;    /* quadrature point on square */
  vector3 t_dx; /* quadrature point on square */
  vector3 t_dy; /* quadrature point on square */
  cubature *Q;  /* quadrature formulas */
  // parametrix *Chi;                      /* parametrizations */
  et_node *pE = disc->mesh->E.patch[0]; /* pointer to the first
                                         * element on the */

  const geometry &Chi = *disc->mesh->geom;

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);
  init_sparse(A, 4 * a_bs * nf, 4 * a_bs * nf, 2 * a_bs);

  /*
   * find first node in element list
   */
  for (i = 0; i < M; ++i) pE = pE[0].son[0];

  c = (double *)calloc(4 * a_bs2, sizeof(double));

  for (zi = 0; zi < nf; ++zi) {
    /*
     * find place on the unit square
     */
    s.x = pE[zi].index_s * h;
    s.y = pE[zi].index_t * h;

    memset(c, 0, 4 * a_bs2 * sizeof(double));

    for (k = 0; k < Q[g].nop; ++k) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
      t_dx = Chi[pE[zi].patch].df_dx(t);
      t_dy = Chi[pE[zi].patch].df_dy(t);
      normal = vector3_norm(Chi[pE[zi].patch].n_f(t));
      disc->VPhi_scal_VPhi(c, Q[g].w[k] * h / normal, Q[g].xi[k], Q[g].xi[k],
                           t_dx, t_dy, t_dx, t_dy);
      disc->Div_Phi_times_Div_Phi(c, Q[g].w[k] / h / normal, Q[g].xi[k],
                                  Q[g].xi[k]);
    }

    for (k = 0; k < a_bs2; ++k) {
      // real
      add_sparse(A, a_bs * zi + k / a_bs, a_bs * zi + k % a_bs, c[k]);
      add_sparse(A, a_bs * zi + k / a_bs, a_bs * (zi + 2 * nf) + k % a_bs,
                 c[k + a_bs2]);
      add_sparse(A, a_bs * (zi + 2 * nf) + k / a_bs, a_bs * zi + k % a_bs,
                 c[k + 2 * a_bs2]);
      add_sparse(A, a_bs * (zi + 2 * nf) + k / a_bs,
                 a_bs * (zi + 2 * nf) + k % a_bs, c[k + 3 * a_bs2]);
      // imaginary
      add_sparse(A, a_bs * (zi + nf) + k / a_bs, a_bs * (zi + nf) + k % a_bs,
                 c[k]);
      add_sparse(A, a_bs * (zi + nf) + k / a_bs,
                 a_bs * (zi + 3 * nf) + k % a_bs, c[k + a_bs2]);
      add_sparse(A, a_bs * (zi + 3 * nf) + k / a_bs,
                 a_bs * (zi + nf) + k % a_bs, c[k + 2 * a_bs2]);
      add_sparse(A, a_bs * (zi + 3 * nf) + k / a_bs,
                 a_bs * (zi + 3 * nf) + k % a_bs, c[k + 3 * a_bs2]);
    }
  }

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(c);

  return;
}
}  // namespace Bembel