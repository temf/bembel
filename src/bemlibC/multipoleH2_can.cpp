// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "multipoleH2.h"

namespace Bembel {
/**
 *  \brief         Computes np Chebychev points.
 *
 *  \param[out]    x              Array of length np with the Chebyshef points
 *  \param[in]     np             Number of Chebyshev points
 *
 */
int cheb_roots(double *x, int np) {
  int i = 0;

  for (i = 1; i <= np; ++i)
    x[np - i] = .5 * (cos(M_PI * (2. * i - 1.) / (2. * np)) + 1.);

  return 0;
}

/**
 *  \brief         Evaluates the kernel function in the Chebyshev points.
 *
 *  \param[out]    ker1           Lower output matrix, row-wise
 *  \param[out]    ker2           Upper output matrix, row-wise
 *  \param[in]     pC1            Element tree node in x direction
 *  \param[in]     pC2            Element tree node in y direction
 *  \param[in]     np             Number of Chebyshev points
 *
 */
int interpol_kernel(hmatrixfactory *hmatfac, double ***ker1, double ***ker2,
                    et_node *pC1, et_node *pC2, int np) {
  int i = 0;
  int j = 0;
  int ell = 0;
  double h = 1. / (1 << pC1->level);
  double *x = NULL;
  vector2 p1;
  vector2 p2;
  pdeproblem *pde = hmatfac->disc->pde;
  const geometry &geom = *hmatfac->disc->mesh->geom;
  double *kappa = pde->kappa;
  int np_fac = pde->np_max_fac;
  int np2 = np * np;
  int k = np2 * np_fac;

  /*
   * init interpolation points
   */
  x = (double *)malloc(np * sizeof(double));
  cheb_roots(x, np);

  for (ell = 0; ell < pde->nct; ++ell) {
    ker1[ell][0] = (double *)malloc(k * k * sizeof(double));
    if (pde->symflags[ell] == 'N')
      ker2[ell][0] = (double *)malloc(k * k * sizeof(double));
  }

  for (i = 0; i < np2; ++i)
    for (j = 0; j < np2; ++j) {
      /*
       * interpolation points are enumerated by i/np and i%np to access
       * all tensor product points
       */
      p1.x = h * (x[i / np] + pC1->index_s);
      p1.y = h * (x[i % np] + pC1->index_t);
      p2.x = h * (x[j / np] + pC2->index_s);
      p2.y = h * (x[j % np] + pC2->index_t);
      pde->interpol_kernel(ker1, ker2, geom[pC1->patch], geom[pC2->patch], p1,
                           p2, i, j, k, kappa);
    }

  free(x);
  return 0;
}

/**
 *  \brief         Computes the coefficients of the Lagrangian polynomials
 *                 corresponding to np Chebychev interpolation points.
 *
 *  \param[out]    L              Coefficients of the np Lagrangian polyonomials
 *  \param[in]     np             Number of Chebychev interpolation points.
 *
 *  The memory for L will be allocated in this routine.
 *
 */
int init_polynomials(double ***L, int np) {
  int i = 0;
  int j = 0;
  int k = 0;
  double *x = NULL;

  /*
   * allocate memory for Chebychev nodes and all Lagrangian polynomials
   */
  x = (double *)malloc(np * sizeof(double));
  (*L) = (double **)calloc(1, np * (sizeof(double *) + np * sizeof(double)));

  /*
   * compute Chebychev nodes
   */
  cheb_roots(x, np);

  for (i = 0; i < np; ++i) {
    /*
     * "allocate" memory for i.th Lagranian polynomial
     */
    (*L)[i] = (double *)(*L + np) + i * np;

    /*
     * compute coefficients for i.th Lagrangian polynomial
     */
    (*L)[i][i] = 1;
    for (j = 1; j < np; ++j)
      for (k = np - 1; k >= j; --k)
        (*L)[i][k] = ((*L)[i][k] - (*L)[i][k - 1]) / (x[k] - x[k - j]);
  }

  /*
   * free memory
   */
  free(x);

  return 0;
}

/**
 *  \brief         Evaluates the polynomial with the np coefficients in L on the
 *                 point xi using the Horner scheme.
 *
 *  \param[in]     L              Coefficients of the polyonomial
 *  \param[in]     np             Number of Coefficients
 *  \param[in]     xi             Evaluation point
 *
 *  \return        Value of the polynomial on xi
 *
 */
double eval_polynomial(double *L, int np, double xi) {
  int i = 0;
  double *x = NULL;
  double val = 0;

  x = (double *)malloc(np * sizeof(double));

  cheb_roots(x, np);
  val = L[np - 1];

  for (i = np - 2; i >= 0; --i) val = val * (xi - x[i]) + L[i];

  free(x);

  return val;
}

/**
 *  \brief         Integrates the one-dimensional polynomials on every
 *                 subinterval of length h. Moments are then obtained by
 *                 tensorizing the 1d moments.
 *
 *  \param[out]    T              Moments of level M
 *  \param[in]     np             Number of Chebychev points
 *  \param[in]     M              Refinement level
 *
 *  Computes only the moments on level M.
 *
 */
int init_moments(double ***T, int np, int l, int M, int a_o,
                 void (*phi)(double *, double, double)) {
  int i = 0;
  int j = 0;
  int k = 0;
  int g = 0;
  int n = 1 << (M - l);
  int N = 1 << l;
  double h = 1. / n;
  double H = 1. / N;
  double w;
  double **L = NULL;
  quadrature *Q;

  /*
   * allocate memory
   */
  (*T) =
      (double **)calloc(1, np * (sizeof(double *) + a_o * n * sizeof(double)));

  /*
   * compute quadrature degree and initialize quadrature routines
   */
  g = .5 * (np + a_o - 2);
  init_Gauss_Legendre(&Q, g + 1);

  /*
   * initialize polynomials
   */
  init_polynomials(&L, np);

  for (i = 0; i < np; ++i) {
    /*
     * "allocate" memory for integrals of i.th polynomial
     */
    (*T)[i] = (double *)(*T + np) + a_o * i * n;

    /*
     * integrate i.th polynomial on j.th subinterval of length h
     */
    for (j = 0; j < n; ++j)
      for (k = 0; k < Q[g].nop; ++k) {
        w = sqrt(H * h) * Q[g].w[k] *
            eval_polynomial(L[i], np, h * (j + Q[g].xi[k]));
        phi(&((*T)[i][a_o * j]), w, Q[g].xi[k]);
      }
  }

  /*
   * free memory
   */
  free(L);
  free_Gauss_Legendre(&Q, g + 1);

  return 0;
}

/**
 *  \brief         Integrates the one-dimensional polynomials on every
 *                 subinterval of length h. Moments are then obtained by
 *                 tensorizing the 1d moments.
 *
 *  \param[out]    T              Moments of all levels
 *  \param[in]     np             Number of Chebychev points
 *  \param[in]     M              Refinement level
 *
 *  Computes the moments on all levels in a recursive way.
 *
 *
 */
int init_all_moments(double ****T, int np, int M, int a_o,
                     void (*phi)(double *, double, double)) {
  int k = 0;

  /*
   * allocate memory for all moments on all levels
   */
  (*T) = (double ***)malloc((M + 1) * sizeof(double **));

  /*
   * compute one-dimensional moments on all levels
   */
  for (k = 0; k <= M; ++k) init_moments((*T) + k, np, k, M, a_o, phi);

  return 0;
}

/**
 *  \brief         Assemble tensor moments for the respective cluster
 *
 *  \param[in]     M              Refinemenet level
 *  \param[in]     lvl            Level up to which to compute the moments
 *  \param[out]    Tmom           Tensor product moments on all levels
 *  \param[in]     np             Number of Chebychev points
 *
 */
int assemble_moments(int M, int lvl, double ***Tmom, int np,
                     discretization *disc,
                     void (**philist)(double *, double, double)) {
  int i = 0;
  int j = 0;
  int ii = 0;
  int k = 0;
  int m1;
  int m2;
  int inc = 1;
  int *index_s = NULL;
  int *index_t = NULL;
  int bs = 1 << (M - lvl);
  int bs2 = bs * bs;
  int np2 = np * np;
  int a_o = disc->a_o;
  int a_bs = disc->a_bs;
  int fac = disc->pde->np_max_fac;
  double ****T1;
  double ****T2;

  if (lvl < 0) {
    Tmom[0] = NULL;
    return 0;
  }

  /*
   * initialize all moment matrices
   */
  T1 = (double ****)calloc(fac, sizeof(double ***));
  T2 = (double ****)calloc(fac, sizeof(double ***));
  for (i = 0; i < fac; ++i) {
    init_all_moments(&T1[i], np, M, a_o, philist[2 * i + 0]);
    init_all_moments(&T2[i], np, M, a_o, philist[2 * i + 1]);
  }

  /*
   * this fancy nice piece of code is able to reconstruct the indices of the
   * patch edge points on the unit square (edge point indices are only
   * level dependent and can thus be recomputed every time it is necessary)
   * thus, we do not depend on the particular patch anymore, but can directly
   * determine the moment for each particular level
   */
  index_s = (int *)calloc(1, bs2 * sizeof(int));
  index_t = (int *)calloc(1, bs2 * sizeof(int));
  inc = bs;
  for (j = 0; j < M - lvl; ++j) {
    inc = 1 << 2 * j;
    for (i = 0; i < inc; ++i) {
      index_s[4 * i] = 2 * index_s[i];
      index_s[4 * i + 1] = 2 * index_s[i] + 1;
      index_s[4 * i + 2] = 2 * index_s[i] + 1;
      index_s[4 * i + 3] = 2 * index_s[i];
      index_t[4 * i] = 2 * index_t[i];
      index_t[4 * i + 1] = 2 * index_t[i];
      index_t[4 * i + 2] = 2 * index_t[i] + 1;
      index_t[4 * i + 3] = 2 * index_t[i] + 1;
    }
  }

  /*
   * assemble tensor-product moments on level lvl
   */
  (*Tmom) = (double **)malloc(np2 * fac * sizeof(double *));
  for (ii = 0; ii < fac; ++ii) {
    for (i = 0; i < np; ++i)
      for (j = 0; j < np; ++j) {
        (*Tmom)[ii * np2 + i * np + j] =
            (double *)malloc(a_bs * bs2 * sizeof(double));
        for (k = 0; k < bs2; ++k)
          for (m1 = 0; m1 < a_o; ++m1)
            for (m2 = 0; m2 < a_o; ++m2)
              (*Tmom)[ii * np2 + i * np + j][a_bs * k + m1 * a_o + m2] =
                  T1[ii][lvl][i][index_s[k] * a_o + m2] *
                  T2[ii][lvl][j][index_t[k] * a_o + m1];
      }
  }

  free(index_s);
  free(index_t);
  for (i = 0; i < fac; ++i) {
    for (j = 0; j <= M; ++j) {
      free(T1[i][j]);
      free(T2[i][j]);
    }
    free(T1[i]);
    free(T2[i]);
  }
  free(T1);
  free(T2);

  return 0;
}

/**
 *  \brief         Initialize one and two dimensional moments.
 *
 *  \param[in]     M              Discretization level.
 *
 */
int init_momentmatrices(discretization *disc, hmatrixsettings *hmatset,
                        double ***Tmom,
                        void (**philist)(double *, double, double)) {
  assemble_moments(disc->mesh->M, disc->mesh->M - hmatset->min_bsize, Tmom,
                   hmatset->np_max, disc, philist);

  return 0;
}

/**
 *  \brief         Assemble transfer matrices
 *
 *  \param[in]     discretization
 *  \param[in]     hmatse[]
 *  \param[out]    Ttr            Transfer matrices, valid for all levels
 *
 */
int init_transfermatrices(discretization *disc, hmatrixsettings *hmatset,
                          double ****Ttr) {
  int i = 0;
  int j = 0;
  int ii = 0;
  int jj = 0;
  int k = 0;
  int np = hmatset->np_max;
  int np2 = np * np;
  double ***E = NULL;
  double *x = NULL;
  double **L = NULL;

  /*
   * initialize Lagrange polynomials and Chebychev roots
   */
  init_polynomials(&L, np);
  x = (double *)malloc(np * sizeof(double));
  cheb_roots(x, np);

  /*
   * initialize values of Lagrange functions on Chebychev points as follows:
   *   * E[0] containes the values of the Lagrange polynomials on the Chebychev
   *          points scaled to [0, 0.5].
   *   * E[1] containes the values of the Lagrange polynomials on the Chebychev
   *          points scaled to [0.5, 1].
   *   * E[*][i][j] containes the value of the j.th polynomial on the i.th point
   */
  E = (double ***)calloc(
      1, 2 * sizeof(double **) +
             2 * np * (sizeof(double *) + 2 * np * np * sizeof(double)));
  E[0] = (double **)E + 2;
  E[1] = E[0] + np;
  for (i = 0; i < np; ++i)
    for (k = 0; k < 2; ++k) {
      E[k][i] = (double *)(E[1] + np) + i * np + k * np * np;
      for (j = 0; j < np; ++j)
        E[k][i][j] = eval_polynomial(L[j], np, 0.5 * (x[i] + k));
    }

  /*
   * allocate space for transfer matrices
   */
  *Ttr = (double ***)calloc(1, 4 * sizeof(double **));
  for (i = 0; i < 4; ++i) {
    (*Ttr)[i] = (double **)calloc(1, np2 * sizeof(double *));
    for (j = 0; j < np2; ++j)
      (*Ttr)[i][j] = (double *)calloc(1, np2 * sizeof(double));
  }

  /*
   * This construction results in the Order R0 R2 R3 R1 to apply to the
   * moments furhermore, the transfer matrices are stored in row major
   * format
   */
  for (k = 0; k < 4; ++k)
    for (i = 0; i < np; ++i)
      for (j = 0; j < np; ++j)
        for (ii = 0; ii < np; ++ii)
          for (jj = 0; jj < np; ++jj)
            (*Ttr)[k][j * np + jj][i * np + ii] =
                E[k / 2][i][j] * E[k % 2][ii][jj];

  /*
   * free memory
   */
  free(E);
  free(L);
  free(x);

  return 0;
}

/**
 *  \brief         Free one and two dimensional moments.
 *
 *  \param[in]     M              Discretization level.
 *
 */
int free_momentmatrices(discretization *disc, hmatrixsettings *hmatset,
                        double **Tmom) {
  int i = 0;

  if (!Tmom) return 0;

  /*
   * free tensor product moments
   */
  for (i = 0; i < hmatset->np_max * hmatset->np_max; ++i) free(Tmom[i]);
  free(Tmom);

  return 0;
}

/**
 *  \brief         Free one and two dimensional moments.
 *
 *  \param[in]     M              Discretization level.
 *
 */
int free_transfermatrices(discretization *disc, hmatrixsettings *hmatset,
                          double ***Ttr) {
  int i = 0;
  int j = 0;

  /*
   * free transfer matrices
   */
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < hmatset->np_max * hmatset->np_max; ++j) free(Ttr[i][j]);
    free(Ttr[i]);
  }
  free(Ttr);

  return 0;
}
}  // namespace Bembel