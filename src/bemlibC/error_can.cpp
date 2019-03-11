// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "error.h"

namespace Bembel {
/**
 *  \brief         Computes the error of the potential in the maximum norm.
 *
 *  \param[in]     Pot            Values of the potential.
 *  \param[in]     Q              Points where the values have been computed.
 *  \param[in]     nq             Length of Pot and Q.
 *
 *  \returns       Error in the maximum norm.
 *
 */
double error(double *Pot, vector3 *Q, int nq) {
  int k;
  double max;
  double min;
  double c;

  printf("                           Maximum potential error is ");

  /*
   * extract constant
   */
  min = max = Pot[0] - f(Q[0]);
  for (k = 1; k < nq; k++) {
    c = Pot[k] - f(Q[k]);
    if (max < c) max = c;
    if (min > c) min = c;
  }

  /*
   * substract constant
   */
  for (k = 0; k < nq; k++) Pot[k] -= 0.5 * (min + max);

  /*
   * measure error
   */
  max = 0;
  for (k = 0; k < nq; k++) {
    c = fabs(Pot[k] - f(Q[k]));
    if (max < c) max = c;
  }

  printf("%g\n", max);

  return max;
}

/**
 *  \brief         Computes the error of the potential in the maximum norm.
 *
 *  \param[in]     Potreal        Real values of the potential.
 *  \param[in]     Potimat        Imaginary values of the potential.
 *  \param[in]     Q              Points where the values have been computed.
 *  \param[in]     nq             Length of Pot and Q.
 *
 *  \returns       Error in the maximum norm.
 *
 */
double errorcomplex(double *Pot, vector3 *Q, int nq, double kappa[3]) {
  int k;
  double max;
  double c;
  double d[2];
  double *Potreal = Pot;
  double *Potimag = Pot + nq;

  printf("                           Maximum potential error is ");

  max = 0;
  for (k = 0; k < nq; k++) {
    Helmholtzf(d, Q[k], kappa);
    d[0] -= Potreal[k];
    d[1] -= Potimag[k];
    c = sqrt(d[0] * d[0] + d[1] * d[1]);
    if (max < c) max = c;
  }

  printf("%g\n", max);

  return max;
}

/**
 *  \brief         Computes the error of the potential in the maximum norm.
 *
 *  \param[in]     Potreal        Real values of the potential.
 *  \param[in]     Potimat        Imaginary values of the potential.
 *  \param[in]     Q              Points where the values have been computed.
 *  \param[in]     nq             Length of Pot and Q.
 *
 *  \returns       Error in the maximum norm.
 *
 */
double errorMaxwell(double *Pot, vector3 *Q, int nq, double kappa[2],
                    void (*Sol)(vector3 d[2], vector3 a, double kappa[2])) {
  int k;
  double max;
  double c;
  vector3 d[2];
  vector3 d2[2];

  printf("                           Maximum potential error is ");

  max = 0;
  for (k = 0; k < nq; k++) {
    Sol(d2, Q[k], kappa);
    d[0].x = d2[0].x - Pot[k];
    d[1].x = d2[1].x - Pot[k + nq];
    d[0].y = d2[0].y - Pot[k + 2 * nq];
    d[1].y = d2[1].y - Pot[k + 3 * nq];
    d[0].z = d2[0].z - Pot[k + 4 * nq];
    d[1].z = d2[1].z - Pot[k + 5 * nq];
    c = sqrt(vector3_skalp(d[0], d[0]) + vector3_skalp(d[1], d[1]));
    if (max < c) max = c;
  }

  printf("%g\n", max);

  return max;
}

/**
 *  \brief         Computes the L^2-error of a function in the tangential space
 *                 obtained by taking the tangential trace.
 *
 *  \param[in]     Coefficient vector
 *  \param[in]     Discretization
 *  \param[in]     Maxwell pde
 *  \param[in]     Function to compare to
 *
 *  \return        L^2 error with inner product from tangential space
 *
 *  Quadrature degree is chosen higher than necessary in order to avoid
 *  superconvergence effects.
 *
 */
double errorTangentialL2TangentialTrace(double *rho, discretization *disc,
                                        pdeproblem *pde,
                                        void (*mf)(vector3[2], vector3,
                                                   double[2])) {
  int i;
  int j;
  int g = disc->g_far;
  int M = disc->mesh->M;
  int n = 1 << M;
  int p = disc->mesh->geom->size();
  int nf = p * n * n;
  int a_bs = disc->a_bs;
  int N = a_bs * nf;
  double *kappa = pde->kappa;
  double e;
  double diff;
  double norm;
  double h = 1. / n;
  double A[3];
  double b[2];
  double w[4];
  double n_y;
  double det;
  double *d;
  double *rhol;
  vector2 s;
  vector2 t;
  vector3 rhoe[2];
  const geometry &Chi = *disc->mesh->geom; /* parametrizations */
  et_node *E = disc->mesh->E.patch[0];     /* pointer to the first element
                                            * on the */
  cubature *Q;
  std::array<vector3, 2> Chi_df;

  /* initialize geometry and quadrature                                      */
  init_Gauss_Square(&Q, g + 1);

  printf("                           L^2-error of density is ");

  /* get pointer to the leafs of the element tree                            */
  while (E->son[0]) E = E->son[0];

  /* project density to local functions                                      */
  rhol = (double *)calloc(6 * N, sizeof(double));
  proj_distr_et(disc, rho, rhol);

  diff = 0;
  norm = 0;
  d = (double *)calloc(a_bs, sizeof(double));
  for (i = 0; i < nf; ++i) {
    s.x = E[i].index_s * h;
    s.y = E[i].index_t * h;
    for (j = 0; j < Q[g].nop; ++j) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[j]));
      n_y = vector3_norm(Chi[E[i].patch].n_f(t));
      Chi_df = Chi[E[i].patch].jacobian(t);
      Chi_df[0] = vector3_Smul(1. / n_y, Chi_df[0]);
      Chi_df[1] = vector3_Smul(1. / n_y, Chi_df[1]);

      /* find coefficents for surface current in tangential space */
      vector3 Chi_n = vector3_Smul(1. / n_y, Chi[E[i].patch].n_f(t));
      mf(rhoe, Chi[E[i].patch].f(t), kappa);
      rhoe[0] = vector3_mul(rhoe[0], Chi_n);
      rhoe[1] = vector3_mul(rhoe[1], Chi_n);
      A[0] = vector3_skalp(Chi_df[0], Chi_df[0]);
      A[1] = vector3_skalp(Chi_df[0], Chi_df[1]);
      A[2] = vector3_skalp(Chi_df[1], Chi_df[1]);
      det = A[0] * A[2] - A[1] * A[1];
      b[0] = vector3_skalp(rhoe[0], Chi_df[0]);
      b[1] = vector3_skalp(rhoe[0], Chi_df[1]);
      w[0] = (A[2] * b[0] - A[1] * b[1]) / det;
      w[2] = (-A[1] * b[0] + A[0] * b[1]) / det;
      b[0] = vector3_skalp(rhoe[1], Chi_df[0]);
      b[1] = vector3_skalp(rhoe[1], Chi_df[1]);
      w[1] = (A[2] * b[0] - A[1] * b[1]) / det;
      w[3] = (-A[1] * b[0] + A[0] * b[1]) / det;

      norm += h * Q[g].w[j] *
              (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);

      /* compute coefficients of numerical approximation in tangential space */
      memset(d, 0, a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[j]);
      w[0] -= myddot(a_bs, rhol + a_bs * i, d);
      w[1] -= myddot(a_bs, rhol + a_bs * i + N, d);
      w[2] -= myddot(a_bs, rhol + a_bs * i + 2 * N, d);
      w[3] -= myddot(a_bs, rhol + a_bs * i + 3 * N, d);

      diff += h * Q[g].w[j] *
              (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
    }
  }

  e = sqrt(diff / norm);
  // e = sqrt(diff);

  printf("%g\n", e);

  free_Gauss_Square(&Q, g + 1);
  free(rhol);
  free(d);

  return e;
}

/**
 *  \brief         Computes the Hdiv-error of a function in the tangential space
 *                 obtained by taking the tangential trace of H and the normal
 *                 component of E.
 *
 *  \param[in]     Coefficient vector
 *  \param[in]     Discretization
 *  \param[in]     Maxwell pde
 *  \param[in]     Hfield
 *  \param[in]     Efield
 *
 *  \return        L^2 error with inner product from tangential space
 *
 *  Quadrature degree is chosen higher than necessary in order to avoid
 *  superconvergence effects.
 *
 */
double errorTangentialHdiv(double *rho, discretization *disc, pdeproblem *pde,
                           void (*Hf)(vector3[2], vector3, double[2]),
                           void (*Ef)(vector3[2], vector3, double[2])) {
  int i;
  int j;
  int g = disc->g_far;
  int M = disc->mesh->M;
  int n = 1 << M;
  int p = disc->mesh->geom->size();
  int nf = p * n * n;
  int a_bs = disc->a_bs;
  int N = a_bs * nf;
  double *kappa = pde->kappa;
  double e;
  double eL2;
  double ediv;
  double diff;
  double norm;
  double diffL2;
  double normL2;
  double diffdiv;
  double normdiv;
  double h = 1. / n;
  double A[3];
  double b[2];
  double w[4];
  double w_div[2];
  double n_y;
  double det;
  double *d;
  double *d_dx;
  double *d_dy;
  double *rhol;
  vector2 s;
  vector2 t;
  vector3 rhoe[2];
  const geometry &Chi = *disc->mesh->geom; /* parametrizations */
  et_node *E = disc->mesh->E.patch[0];     /* pointer to the first element
                                            * on the */
  cubature *Q;
  vector3 Chi_n;
  std::array<vector3, 2> Chi_df;

  /* initialize geometry and quadrature                                      */
  init_Gauss_Square(&Q, g + 1);

  printf("                           Hdiv-error of surface current is ");

  /* get pointer to the leafs of the element tree                            */
  while (E->son[0]) E = E->son[0];

  /* project density to local functions                                      */
  rhol = (double *)calloc(4 * N, sizeof(double));
  proj_distr_et(disc, rho, rhol);
  // for (i = 0 ; i < 4 * N; ++i)
  // rhol[i] *= 4 * M_PI;

  diff = 0;
  norm = 0;
  diffL2 = 0;
  normL2 = 0;
  diffdiv = 0;
  normdiv = 0;
  d = (double *)calloc(a_bs, sizeof(double));
  d_dx = (double *)calloc(a_bs, sizeof(double));
  d_dy = (double *)calloc(a_bs, sizeof(double));
  for (i = 0; i < nf; ++i) {
    s.x = E[i].index_s * h;
    s.y = E[i].index_t * h;
    for (j = 0; j < Q[g].nop; ++j) {
      t = vector2_add(s, vector2_Smul(h, Q[g].xi[j]));
      Chi_n = Chi[E[i].patch].n_f(t);
      n_y = vector3_norm(Chi[E[i].patch].n_f(t));
      Chi_n = vector3_Smul(1. / n_y, Chi_n);
      Chi_df = Chi[E[i].patch].jacobian(t);
      Chi_df[0] = vector3_Smul(1. / n_y, Chi_df[0]);
      Chi_df[1] = vector3_Smul(1. / n_y, Chi_df[1]);

      /* find coefficents for surface current in tangential space */
      Hf(rhoe, Chi[E[i].patch].f(t), kappa);
      rhoe[0] = vector3_mul(rhoe[0], Chi_n);
      rhoe[1] = vector3_mul(rhoe[1], Chi_n);
      A[0] = vector3_skalp(Chi_df[0], Chi_df[0]);
      A[1] = vector3_skalp(Chi_df[0], Chi_df[1]);
      A[2] = vector3_skalp(Chi_df[1], Chi_df[1]);
      det = A[0] * A[2] - A[1] * A[1];
      b[0] = vector3_skalp(rhoe[0], Chi_df[0]);
      b[1] = vector3_skalp(rhoe[0], Chi_df[1]);
      w[0] = (A[2] * b[0] - A[1] * b[1]) / det;
      w[2] = (-A[1] * b[0] + A[0] * b[1]) / det;
      b[0] = vector3_skalp(rhoe[1], Chi_df[0]);
      b[1] = vector3_skalp(rhoe[1], Chi_df[1]);
      w[1] = (A[2] * b[0] - A[1] * b[1]) / det;
      w[3] = (-A[1] * b[0] + A[0] * b[1]) / det;

      norm += h * Q[g].w[j] *
              (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
      normL2 += h * Q[g].w[j] *
                (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);

      /* compute coefficients of numerical approximation in tangential space */
      memset(d, 0, a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[j]);
      w[0] -= myddot(a_bs, rhol + a_bs * i, d);
      w[1] -= myddot(a_bs, rhol + a_bs * i + N, d);
      w[2] -= myddot(a_bs, rhol + a_bs * i + 2 * N, d);
      w[3] -= myddot(a_bs, rhol + a_bs * i + 3 * N, d);

      diff += h * Q[g].w[j] *
              (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
      diffL2 += h * Q[g].w[j] *
                (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);

#if 1
      /* compute normal component of Efield */
      Ef(rhoe, Chi[E[i].patch].f(t), kappa);
      w_div[0] = vector3_skalp(rhoe[0], Chi_n);
      w_div[1] = vector3_skalp(rhoe[1], Chi_n);

      norm += h * n_y * Q[g].w[j] * (w_div[0] * w_div[0] + w_div[1] * w_div[1]);
      normdiv +=
          h * n_y * Q[g].w[j] * (w_div[0] * w_div[0] + w_div[1] * w_div[1]);

      /* compute surface divergence of numerical approximation */
      memset(d_dx, 0, a_bs * sizeof(double));
      disc->phiphi_dx(d_dx, Q[g].xi[j]);
      w_div[0] -= myddot(a_bs, rhol + a_bs * i, d_dx);
      w_div[1] -= myddot(a_bs, rhol + a_bs * i + N, d_dx);
      memset(d_dy, 0, a_bs * sizeof(double));
      disc->phiphi_dy(d_dy, Q[g].xi[j]);
      w_div[0] -= myddot(a_bs, rhol + a_bs * i + 2 * N, d_dy);
      w_div[1] -= myddot(a_bs, rhol + a_bs * i + 3 * N, d_dy);

      diff += Q[g].w[j] * (w_div[0] * w_div[0] + w_div[1] * w_div[1]) / n_y;
      diffdiv += Q[g].w[j] * (w_div[0] * w_div[0] + w_div[1] * w_div[1]) / n_y;
#endif
    }
  }

  e = sqrt(diff / norm);
  eL2 = sqrt(diffL2 / norm);
  ediv = sqrt(diffdiv / norm);

  printf("%g, (L2: %g, div: %g)\n", e, eL2, ediv);

  free_Gauss_Square(&Q, g + 1);
  free(rhol);
  free(d);

  return e;
}
}  // namespace Bembel