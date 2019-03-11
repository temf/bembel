// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "Rhs.hpp"

namespace Bembel {
namespace Rhs {
/**
 *  \brief        A simple wrapper for the initialization of the right hand
 * side.
 */
Eigen::VectorXcd computeRhs(
    Discretization<MaxwellSingle> &ddisc,
    std::function<Eigen::Vector3cd(Eigen::Vector3d, std::complex<double>)>
        fun) {
  discretization *disc = &(ddisc.get_disc());
  auto local_data = [fun](vector3 field[2], vector3 &&X, double kapl[2]) {
    Eigen::Vector3cd tmp = fun(Eigen::Vector3d(X.x, X.y, X.z),
                               std::complex<double>(kapl[0], kapl[1]));
    field[0].x = tmp(0).real();
    field[0].y = tmp(1).real();
    field[0].z = tmp(2).real();
    field[1].x = tmp(0).imag();
    field[1].y = tmp(1).imag();
    field[1].z = tmp(2).imag();
    return;
  };

  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double h = 1. / n;          /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double c[4];
  double *d;
  double *rhs1l;
  double *rhs1lc;
  double *rhs2l;
  double *rhs2lc;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  std::array<vector3, 2> Chi_df;
  vector3 Chi_n;
  vector3 field[2];
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhs1l = (double *)calloc(nf * 4 * disc->a_bs, sizeof(double));
  rhs1lc = rhs1l + nf * disc->a_bs;
  rhs2l = rhs1lc + nf * disc->a_bs;
  rhs2lc = rhs2l + nf * disc->a_bs;

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
      local_data(field, Chi[pE[zi].patch].f(t), kappa);
      w = h * Q[g].w[k];
      Chi_n = Chi[pE[zi].patch].n_f(t);
      Chi_n = vector3_Smul(1. / vector3_norm(Chi_n), Chi_n);
      field[0] = vector3_mul(Chi_n, vector3_mul(field[0], Chi_n));
      field[1] = vector3_mul(Chi_n, vector3_mul(field[1], Chi_n));
      Chi_df = Chi[pE[zi].patch].jacobian(t);
      c[0] = w * vector3_skalp(Chi_df[0], field[0]);
      c[1] = w * vector3_skalp(Chi_df[0], field[1]);
      c[2] = w * vector3_skalp(Chi_df[1], field[0]);
      c[3] = w * vector3_skalp(Chi_df[1], field[1]);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[2], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[3], d, rhs2lc + disc->a_bs * zi);
    }
  }

  double *rhs = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhs1l, rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhs1l);
  free(d);

  Eigen::VectorXcd out = cmplxptr2eigen(rhs, disc->na);

  free(rhs);
  return out;
}

/**
 *  \brief         Initializes the right hand side of the Galerkin method for
 *                 the canonical scalar product.
 *
 *  \param[in,out] rhs            Place where to store the rhs.
 *  \param[in]     E              Element tree.
 *  \param[in]     M              Discretization level.
 *  \param[in]     na             Number of ansatz functions.
 */
void BEMRHS_MaxwellDirichletTangentialTrace(double **rhs, discretization *disc,
                                            void (*f)(vector3 field[2],
                                                      vector3 X,
                                                      double kappa[2])) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double h = 1. / n;          /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double c[4];
  double *d;
  double *rhs1l;
  double *rhs1lc;
  double *rhs2l;
  double *rhs2lc;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  std::array<vector3, 2> Chi_df;
  vector3 Chi_n;
  vector3 field[2];
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhs1l = (double *)calloc(nf * 4 * disc->a_bs, sizeof(double));
  rhs1lc = rhs1l + nf * disc->a_bs;
  rhs2l = rhs1lc + nf * disc->a_bs;
  rhs2lc = rhs2l + nf * disc->a_bs;

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
      f(field, Chi[pE[zi].patch].f(t), kappa);
      w = h * Q[g].w[k];
      Chi_n = Chi[pE[zi].patch].n_f(t);
      Chi_n = vector3_Smul(1. / vector3_norm(Chi_n), Chi_n);
      field[0] = vector3_mul(vector3_mul(field[0], Chi_n), Chi_n);
      field[1] = vector3_mul(vector3_mul(field[1], Chi_n), Chi_n);
      Chi_df = Chi[pE[zi].patch].jacobian(t);
      c[0] = w * vector3_skalp(Chi_df[0], field[0]);
      c[1] = w * vector3_skalp(Chi_df[0], field[1]);
      c[2] = w * vector3_skalp(Chi_df[1], field[0]);
      c[3] = w * vector3_skalp(Chi_df[1], field[1]);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[2], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[3], d, rhs2lc + disc->a_bs * zi);
    }
  }

  (*rhs) = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhs1l, *rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhs1l);
  free(d);
  return;
}

/**
 *  \brief         Initializes the right hand side of the Galerkin method for
 *                 the canonical scalar product.
 *
 *  \param[in,out] rhs            Place where to store the rhs.
 *  \param[in]     E              Element tree.
 *  \param[in]     M              Discretization level.
 *  \param[in]     na             Number of ansatz functions.
 */
void BEMRHS_MaxwellDirichletTrace(double **rhs, discretization *disc,
                                  void (*f)(vector3 field[2], vector3 X,
                                            double kappa[2])) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double h = 1. / n;          /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double c[4];
  double *d;
  double *rhs1l;
  double *rhs1lc;
  double *rhs2l;
  double *rhs2lc;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  std::array<vector3, 2> Chi_df;
  vector3 Chi_n;
  vector3 field[2];
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhs1l = (double *)calloc(nf * 4 * disc->a_bs, sizeof(double));
  rhs1lc = rhs1l + nf * disc->a_bs;
  rhs2l = rhs1lc + nf * disc->a_bs;
  rhs2lc = rhs2l + nf * disc->a_bs;

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
      Chi_n = Chi[pE[zi].patch].n_f(t);
      Chi_n = vector3_Smul(1. / vector3_norm(Chi_n), Chi_n);
      Chi_df = Chi[pE[zi].patch].jacobian(t);
      f(field, Chi[pE[zi].patch].f(t), kappa);
      w = h * Q[g].w[k];
      field[0] = vector3_mul(field[0], Chi_n);
      field[1] = vector3_mul(field[1], Chi_n);
      c[0] = w * vector3_skalp(Chi_df[0], field[0]);
      c[1] = w * vector3_skalp(Chi_df[0], field[1]);
      c[2] = w * vector3_skalp(Chi_df[1], field[0]);
      c[3] = w * vector3_skalp(Chi_df[1], field[1]);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[2], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[3], d, rhs2lc + disc->a_bs * zi);
    }
  }

  (*rhs) = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhs1l, *rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhs1l);
  free(d);
  return;
}

/**
 *  \brief         Initializes the right hand side of the Galerkin method for
 *                 the canonical scalar product.
 *
 *  \param[in,out] rhs            Place where to store the rhs.
 *  \param[in]     E              Element tree.
 *  \param[in]     M              Discretization level.
 *  \param[in]     na             Number of ansatz functions.
 */
void BEMRHS_TangentialL2(double **rhs, discretization *disc,
                         void (*f)(vector3 field[2], vector3 X,
                                   double kappa[2])) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double h = 1. / n;          /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double c[4];
  double *d;
  double *rhs1l;
  double *rhs1lc;
  double *rhs2l;
  double *rhs2lc;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  std::array<vector3, 2> Chi_df;
  vector3 Chi_n;
  vector3 field[2];
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhs1l = (double *)calloc(nf * 4 * disc->a_bs, sizeof(double));
  rhs1lc = rhs1l + nf * disc->a_bs;
  rhs2l = rhs1lc + nf * disc->a_bs;
  rhs2lc = rhs2l + nf * disc->a_bs;

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
      f(field, Chi[pE[zi].patch].f(t), kappa);
      w = h * Q[g].w[k];
      Chi_n = Chi[pE[zi].patch].n_f(t);
      Chi_n = vector3_Smul(1. / vector3_norm(Chi_n), Chi_n);
      field[0] = vector3_mul(field[0], Chi_n);
      field[1] = vector3_mul(field[1], Chi_n);
      Chi_df = Chi[pE[zi].patch].jacobian(t);
      c[0] = w * vector3_skalp(Chi_df[0], field[0]);
      c[1] = w * vector3_skalp(Chi_df[0], field[1]);
      c[2] = w * vector3_skalp(Chi_df[1], field[0]);
      c[3] = w * vector3_skalp(Chi_df[1], field[1]);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[2], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[3], d, rhs2lc + disc->a_bs * zi);
    }
  }

  (*rhs) = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhs1l, *rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhs1l);
  free(d);
  return;
}

/**
 *  \brief         Initializes the right hand side of the Galerkin method for
 *                 the surface-Hdiv scalar product. When solved, the projection
 *                 is an approximation to the tangential trace of f.
 *
 *  \param[in,out] rhs            Place where to store the rhs.
 *  \param[in]     E              Element tree.
 *  \param[in]     M              Discretization level.
 *  \param[in]     na             Number of ansatz functions.
 */
void BEMRHS_Hdiv(double **rhs, discretization *disc,
                 void (*f)(vector3 field[2], vector3 X, double kappa[2]),
                 void (*fcurl)(vector3 field[2], vector3 X, double kappa[2])) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int k;                      /* increment variable */
  int g = disc->g_far;        /* quadrature degree */
  int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M; /* n*n patches per parameter
                               * domain */
  int nf = p * n * n;         /* number of elements */
  int zi;                     /* row index of the element */
  double h = 1. / n;          /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double c[4];
  double *d;
  double *rhs1l;
  double *rhs1lc;
  double *rhs2l;
  double *rhs2lc;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  std::array<vector3, 2> Chi_df;
  vector3 Chi_n;
  vector3 Chi_f;
  vector3 field[2];
  vector3 fieldCurl[2];
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhs1l = (double *)calloc(nf * 4 * disc->a_bs, sizeof(double));
  rhs1lc = rhs1l + nf * disc->a_bs;
  rhs2l = rhs1lc + nf * disc->a_bs;
  rhs2lc = rhs2l + nf * disc->a_bs;

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
      Chi_f = Chi[pE[zi].patch].f(t);
      Chi_n = Chi[pE[zi].patch].n_f(t);
      Chi_n = vector3_Smul(1. / vector3_norm(Chi_n), Chi_n);
      Chi_df = Chi[pE[zi].patch].jacobian(t);
      f(field, Chi_f, kappa);
      w = h * Q[g].w[k];
      field[0] = vector3_mul(field[0], Chi_n);
      field[1] = vector3_mul(field[1], Chi_n);
      c[0] = w * vector3_skalp(Chi_df[0], field[0]);
      c[1] = w * vector3_skalp(Chi_df[0], field[1]);
      c[2] = w * vector3_skalp(Chi_df[1], field[0]);
      c[3] = w * vector3_skalp(Chi_df[1], field[1]);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[2], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[3], d, rhs2lc + disc->a_bs * zi);
#if 1
      fcurl(fieldCurl, Chi_f, kappa);
      w = Q[g].w[k];
      c[0] = w * vector3_skalp(fieldCurl[0], Chi_n);
      c[1] = w * vector3_skalp(fieldCurl[1], Chi_n);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi_dx(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs1l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs1lc + disc->a_bs * zi);
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi_dy(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhs2l + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhs2lc + disc->a_bs * zi);
#endif
    }
  }

  (*rhs) = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhs1l, *rhs);

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhs1l);
  free(d);
  return;
}
}  // namespace Rhs
}  // namespace Bembel
