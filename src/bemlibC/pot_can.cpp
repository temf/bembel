// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "pot.h"

namespace Bembel {
/**
 *  \brief         Evaluates the potential in the points in R.
 *
 *  \param[in]     rho              Density of the integral equation.
 *  \param[out]    Pot            Result.
 *  \param[in]     R              Points where the potential has to be
 *                                evaluated.
 *  \param[in]     nr             Length of R.
 *  \param[in]     disc           Discretization.
 *
 */
void pot(double *rho, double **Pot, vector3 *R, int nr, discretization *disc,
         pdeproblem *pde) {
  meshdata *mesh = disc->mesh;
  int i;                            /* increment variable */
  int k;                            /* rotation and increment variable */
  int g = disc->g_pot;              /* quadrature degree */
  const int p = mesh->geom->size(); /* number of parameter domains */
  int n = 1 << mesh->M;             /* n*n element per parameter domain */
  int nf = p * n * n;               /* total number of element */
  int zi;                           /* row index of the element */
  double h = 1. / n;                /* step size */
  double *d;
  double *rhol;
  double *mypot;
  vector2 s;                         /* left lower corner on square for zi */
  vector2 t;                         /* quadrature point on square */
  vector3 y;                         /* point on surface */
  vector3 n_y;                       /* point on surface */
  cubature *Q;                       /* quadrature formulas */
  const geometry &Chi = *mesh->geom; /* parametrizations */
  et_node *pE = mesh->E.patch[0];    /* pointer to the first element
                                      * on the */
  /*
   * lowest level of the element tree
   */

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1); /* Kubatur-Formeln */

  /*
   * allocate memory
   */
  (*Pot) = (double *)calloc(nr, sizeof(double));
  rhol = (double *)calloc(disc->a_bs * nf, sizeof(double));

  /*
   * find first node in element list
   */
  for (i = 0; i < mesh->M; ++i) pE = pE[0].son[0];

  proj_distr_et(disc, rho, rhol);

/*
 * Do a parallelization with subsequent reduction. This can be done
 * wonderful in OpenMP 4.0 and later, but to be backwards compatible, we do
 * it the old school way
 */
#pragma omp parallel default(shared) private(d, zi, s, k, t, y, n_y, mypot)
  {
    mypot = (double *)calloc(nr, sizeof(double));
    d = (double *)calloc(disc->a_bs, sizeof(double));
#pragma omp for
    for (zi = 0; zi < nf; ++zi) {
      /*
       * find place on the unit square
       */
      s.x = pE[zi].index_s * h;
      s.y = pE[zi].index_t * h;

      /*
       * Der Doppelschicht hat das Minuszeichen!
       */
      for (k = 0; k < Q[g].nop; ++k) {
        t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
        y = Chi[pE[zi].patch].f(t);
        n_y = Chi[pE[zi].patch].n_f(t);
        disc->phiphi(d, Q[g].xi[k]);
        pde->pot_eval(mypot, R, nr, Q[g].w[k], y, n_y, rhol + disc->a_bs * zi,
                      h, d, disc);
      }
    }
    free(d);
#pragma omp critical
    {
      for (zi = 0; zi < nr; ++zi) (*Pot)[zi] += mypot[zi];
    }
    free(mypot);
  }

  /*
   * normalize solution
   */
  pde->pot_normalize(*Pot, nr, h);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhol);

  return;
}

/**
 *  \brief         Evaluates the potential in the points in R.
 *
 *  \param[in]     rho              Density of the integral equation.
 *  \param[out]    Pot            Result.
 *  \param[in]     R              Points where the potential has to be
 *                                evaluated.
 *  \param[in]     nr             Length of R.
 *  \param[in]     disc           Discretization.
 *
 */
void potcomplex(double *rho, double **Pot, vector3 *R, int nr,
                discretization *disc, pdeproblem *pde) {
  meshdata *mesh = disc->mesh;
  int i;                            /* increment variable */
  int k;                            /* rotation and increment variable */
  int g = disc->g_pot;              /* quadrature degree */
  const int p = mesh->geom->size(); /* number of parameter domains */
  int n = 1 << mesh->M;             /* n*n element per parameter domain */
  int nf = p * n * n;               /* total number of element */
  int zi;                           /* row index of the element */
  double h = 1. / n;                /* step size */
  double *d;
  double *rhol;
  double *mypot;
  vector2 s;   /* left lower corner on square for zi */
  vector2 t;   /* quadrature point on square */
  vector3 y;   /* point on surface */
  vector3 n_y; /* point on surface */
  cubature *Q; /* quadrature formulas */
  // parametrix *Chi = mesh->geom->Chi; /* parametrizations */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the lowest level of the
                                   * element tree */

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1); /* Kubatur-Formeln */

  /*
   * allocate memory
   */
  (*Pot) = (double *)calloc(2 * nr, sizeof(double));
  rhol = (double *)calloc(2 * disc->a_bs * nf, sizeof(double));

  /*
   * find first node in element list
   */
  for (i = 0; i < mesh->M; ++i) pE = pE[0].son[0];

  proj_distr_et(disc, rho, rhol);

/*
 * Do a parallelization with subsequent reduction. This can be done
 * wonderful in OpenMP 4.0 and later, but to be backwards compatible, we do
 * it the old school way
 */
#pragma omp parallel default(shared) private(d, zi, s, k, t, y, n_y, mypot)
  {
    mypot = (double *)calloc(2 * nr, sizeof(double));
    d = (double *)calloc(disc->a_bs, sizeof(double));
#pragma omp for
    for (zi = 0; zi < nf; ++zi) {
      /*
       * find place on the unit square
       */
      s.x = pE[zi].index_s * h;
      s.y = pE[zi].index_t * h;

      /*
       * iterate over quadrature points
       */
      for (k = 0; k < Q[g].nop; ++k) {
        t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
        y = Chi[pE[zi].patch].f(t);
        n_y = Chi[pE[zi].patch].n_f(t);
        disc->phiphi(d, Q[g].xi[k]);
        pde->pot_eval(mypot, R, nr, Q[g].w[k], y, n_y, rhol + disc->a_bs * zi,
                      h, d, disc);
      }
    }
    free(d);
#pragma omp critical
    {
      for (zi = 0; zi < 2 * nr; ++zi) (*Pot)[zi] += mypot[zi];
    }
    free(mypot);
  }

  /*
   * normalize solution
   */
  pde->pot_normalize(*Pot, nr, h);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhol);

  return;
}

std::vector<double> maxwellTpPotential(std::vector<std::vector<double>> &tps,
                                       double *rho, discretization *disc,
                                       pdeproblem *pde) {
  assert(tps.size() == 3);
  const int nx = tps[0].size();
  const int ny = tps[1].size();
  const int nz = tps[2].size();
  const int sz = nx * nz * ny;

  double *Pot;

  vector3 *Q;
  Q = (vector3 *)malloc(sz * sizeof(vector3));

  int count = 0;

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        Q[count].x = tps[0][ix];
        Q[count].y = tps[1][iy];
        Q[count].z = tps[2][iz];
        count++;
      }
  potMaxwell(rho, &Pot, Q, count, disc, pde);

  std::vector<double> pots(count);
  for (int in = 0; in < count; in++)
    pots[in] =
        std::sqrt(Pot[in] * Pot[in] + Pot[in + 2 * sz] * Pot[in + 2 * sz] +
                  Pot[in + 4 * sz] * Pot[in + 4 * sz]);

  free(Q);
  return pots;
}

std::vector<vector3> maxwellTpVPotential(std::vector<std::vector<double>> &tps,
                                         double *rho, discretization *disc,
                                         pdeproblem *pde) {
  assert(tps.size() == 3);
  const int nx = tps[0].size();
  const int ny = tps[1].size();
  const int nz = tps[2].size();
  const int sz = nx * nz * ny;

  double *Pot;

  vector3 *Q;
  Q = (vector3 *)malloc(sz * sizeof(vector3));

  int count = 0;

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        Q[count].x = tps[0][ix];
        Q[count].y = tps[1][iy];
        Q[count].z = tps[2][iz];
        count++;
      }
  potMaxwell(rho, &Pot, Q, count, disc, pde);

  std::vector<vector3> pots;
  pots.reserve(count);
  for (int in = 0; in < count; in++)
    pots.push_back(vector3_make(Pot[in], Pot[in + 2 * nx * ny * nz],
                                Pot[in + 4 * nx * ny * nz]));
  free(Q);
  free(Pot);
  return pots;
}
}  // namespace Bembel