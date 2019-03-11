// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "H2_level2.h"
#include "cluster_tree.h"
#include "discretization.h"
#include "gram.h"
#include "pdeproblem.h"

namespace Bembel {

static inline double LaplaceSingleLayerKernel(vector3 x, vector3 y) {
  return (1. / (4 * M_PI *
                sqrt((x.x - y.x) * (x.x - y.x) + (x.y - y.y) * (x.y - y.y) +
                     (x.z - y.z) * (x.z - y.z))));
}

static inline int LaplaceSingle_quadrature_accuracy(int a_o) {
  return 3 + 2 * (a_o - 1);
}

static inline int LaplaceSingle_g_far(int a_o) { return a_o + 1; }

static inline int LaplaceSingle_g_pot(int a_o) { return a_o + 1; }

static void LaplaceSingle_init_randwerte(
    randwerte ***RW,     /* zu berechnende Randwerte *
                          * fuer Fernfeld-Quadratur */
    cubature *Q,         /* Kubatur-Formel */
    const geometry &Chi, /* definiert die
                          * Parametrisierung */
    int n)               /* n*n Patches pro Parametergebiet */
{
  const int p = Chi.size(); /* Anzahl der Parametergebiete */
  int i1, i2, i3;           /* Laufindizes fuer das Patch */
  int zi;                   /* Zeilenindex hieraus: zi = i1*(m*m)+i2*m+i3 */
  int k;                    /* Laufindex fuer Gauss-Quadratur */
  vector2 s;                /* Linker, unterer Echpunkt von Patch zi */
  vector2 t;                /* Stuetzpunkt der Tensor-Gauss-Quadratur auf
                             * Q */
  double h = 1. / n;        /* Schrittweite */

  (*RW) = (randwerte **)malloc(p * n * n * sizeof(randwerte *));
  zi = 0;
  for (i1 = 0; i1 < p; i1++) {
    s.y = 0;
    for (i2 = 0; i2 < n; i2++) {
      s.x = 0;
      for (i3 = 0; i3 < n; i3++) {
        (*RW)[zi] = (randwerte *)malloc(Q->nop * sizeof(randwerte));
        for (k = 0; k < Q->nop; k++) {
          t = vector2_add(s, vector2_Smul(h, Q->xi[k]));
          (*RW)[zi][k].Chi = Chi[i1].f(t);
          (*RW)[zi][k].det_dChi = h * Q->w[k] * vector3_norm(Chi[i1].n_f(t));
        }
        s.x += h;
        zi++;
      }
      s.y += h;
    }
  }
  return;
}

/**
 * Gibt den Speicherplatz fuer die Randwerte frei
 *
 */
static void LaplaceSingle_free_randwerte(
    randwerte ***RW,    /* Randwerte */
    int n, const int p) /* n*n Patches pro Parametergebiet */
{
  int k; /* Laufindex */
  for (k = 0; k < p * n * n; k++) free((*RW)[k]);
  free(*RW);

  return;
}

/**
 *  Fernfeld-Quadratur-Routine
 *
 */
static void LaplaceSingle_IntPhi0(double *c, int no1, int no2, cubature *Q,
                                  randwerte **RW, void *disc) {
  int i, j; /* Laufindizes */
  double d; /* temporaere Groesse */
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    for (j = 0; j < Q->nop; j++) {
      d = RW[no1][i].det_dChi * RW[no2][j].det_dChi *
          LaplaceSingleLayerKernel(RW[no1][i].Chi, RW[no2][j].Chi);
      Phi_times_Phi(c, d, Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * No-Problem-Quadrature-Routine -> kanonisches Skalarprodukt
 *
 */
static void LaplaceSingle_IntPhi1(double *c, vector2 s, vector2 t, double h,
                                  const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double d, w;
  vector2 xi, eta;
  vector3 x;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = vector2_add(s, vector2_Smul(h, Q->xi[i]));
    x = Chi_s.f(xi);
    w = h * h * Q->w[i] * vector3_norm(Chi_s.n_f(xi));
    for (j = 0; j < Q->nop; j++) {
      eta = vector2_add(t, vector2_Smul(h, Q->xi[j]));
      d = w * Q->w[j] * vector3_norm(Chi_t.n_f(eta)) *
          LaplaceSingleLayerKernel(x, Chi_t.f(eta));
      Phi_times_Phi(c, d, Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * GLEICHE PATCHES [0,1]^2 -> kanonisches Skalarprodukt
 *
 */
static void LaplaceSingle_IntPhi2(double *c, vector2 s, double h,
                                  const Spl::Patch &Chi_s, cubature *Q,
                                  void *disc) {
  int i, j;
  double d, w, t1, t2, t3, t4;
  vector2 xi, eta, a, b;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    w = h * h * Q->w[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
    for (j = 0; j < Q->nop; j++) {
      eta = Q->xi[j];
      t1 = eta.x * (1 - xi.x);
      t2 = eta.y * (1 - xi.x * xi.y);
      t3 = t1 + xi.x;
      t4 = t2 + xi.x * xi.y;

      a.x = s.x + h * t1;
      a.y = s.y + h * t2;
      b.x = s.x + h * t3;
      b.y = s.y + h * t4;
      d = w * Q->w[j] * LaplaceSingleLayerKernel(Chi_s.f(a), Chi_s.f(b)) *
          vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      Phi_times_Phi(c, d, vector2_make(t1, t2), vector2_make(t3, t4));
      Phi_times_Phi(c, d, vector2_make(t3, t4), vector2_make(t1, t2));

      a.y = s.y + h * t4;
      b.y = s.y + h * t2;
      d = w * Q->w[j] * LaplaceSingleLayerKernel(Chi_s.f(a), Chi_s.f(b)) *
          vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      Phi_times_Phi(c, d, vector2_make(t1, t4), vector2_make(t3, t2));
      Phi_times_Phi(c, d, vector2_make(t3, t2), vector2_make(t1, t4));

      a.x = s.x + h * t2;
      a.y = s.y + h * t1;
      b.x = s.x + h * t4;
      b.y = s.y + h * t3;
      d = w * Q->w[j] * LaplaceSingleLayerKernel(Chi_s.f(a), Chi_s.f(b)) *
          vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      Phi_times_Phi(c, d, vector2_make(t2, t1), vector2_make(t4, t3));
      Phi_times_Phi(c, d, vector2_make(t4, t3), vector2_make(t2, t1));

      a.y = s.y + h * t3;
      b.y = s.y + h * t1;
      d = w * Q->w[j] * LaplaceSingleLayerKernel(Chi_s.f(a), Chi_s.f(b)) *
          vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      Phi_times_Phi(c, d, vector2_make(t2, t3), vector2_make(t4, t1));
      Phi_times_Phi(c, d, vector2_make(t4, t1), vector2_make(t2, t3));
    }
  }

  return;
}

/**
 * GEMEINSAME KANTE [0,1] -> kanonisches Skalarprodukt
 *
 */
static void LaplaceSingle_IntPhi3(double *c, vector2 s, vector2 t, double h,
                                  int ind_s, int ind_t, const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double d, w, t1, t2, t3, t4;
  vector2 xi, eta, a, b, u, v;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;
  memset(c, 0, a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    w = h * h * xi.y * xi.y * Q->w[i];
    t1 = xi.x * (1 - xi.y);
    t2 = (1 - xi.x) * (1 - xi.y);

    for (j = 0; j < Q->nop; j++) {
      eta = vector2_Smul(xi.y, Q->xi[j]);
      t3 = xi.x * (1 - eta.x);
      t4 = (1 - xi.x) * (1 - eta.x);

      a = Tau(t1, eta.x, ind_s);
      b = Tau(t2, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - xi.y) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);

      a = Tau(1 - t1, eta.x, ind_s);
      b = Tau(1 - t2, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - xi.y) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);

      a = Tau(t3, xi.y, ind_s);
      b = Tau(t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - eta.x) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);

      a = Tau(1 - t3, xi.y, ind_s);
      b = Tau(1 - t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - eta.x) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);

      a = Tau(t4, eta.y, ind_s);
      b = Tau(t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - eta.x) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);

      a = Tau(1 - t4, eta.y, ind_s);
      b = Tau(1 - t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      d = w * Q->w[j] * (1 - eta.x) *
          LaplaceSingleLayerKernel(Chi_s.f(u), Chi_t.f(v)) *
          vector3_norm(Chi_s.n_f(u)) * vector3_norm(Chi_t.n_f(v));
      Phi_times_Phi(c, d, a, b);
    }
  }

  return;
}

/**
 * GEMEINSAME ECKE IM NULLPUNKT -> kanonisches Skalarprodukt
 *
 */
static void LaplaceSingle_IntPhi4(double *c, vector2 s, vector2 t, double h,
                                  int ind_s, int ind_t, const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double w, n_x1, n_x2, n_y1, n_y2, n_z;
  vector2 xi, eta, a, u, a1, a2, b1, b2;
  vector3 x1, x2, y1, y2, z;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    w = h * h * pow(xi.x, 3) * Q->w[i];
    xi.y *= xi.x;
    a1 = Tau(xi.x, xi.y, ind_s);
    a2 = Tau(xi.y, xi.x, ind_s);
    b1 = Tau(xi.x, xi.y, ind_t);
    b2 = Tau(xi.y, xi.x, ind_t);

    u = Kappa(s, a1, h);
    x1 = Chi_s.f(u);
    n_x1 = vector3_norm(Chi_s.n_f(u));

    u = Kappa(s, a2, h);
    x2 = Chi_s.f(u);
    n_x2 = vector3_norm(Chi_s.n_f(u));

    u = Kappa(t, b1, h);
    y1 = Chi_t.f(u);
    n_y1 = vector3_norm(Chi_t.n_f(u));

    u = Kappa(t, b2, h);
    y2 = Chi_t.f(u);
    n_y2 = vector3_norm(Chi_t.n_f(u));

    for (j = 0; j < Q->nop; j++) {
      eta = vector2_Smul(xi.x, Q->xi[j]);

      a = Tau(eta.x, eta.y, ind_t);
      u = Kappa(t, a, h);
      z = Chi_t.f(u);
      n_z = vector3_norm(Chi_t.n_f(u));
      Phi_times_Phi(
          c, w * Q->w[j] * n_z * n_x1 * LaplaceSingleLayerKernel(x1, z), a1, a);
      Phi_times_Phi(
          c, w * Q->w[j] * n_z * n_x2 * LaplaceSingleLayerKernel(x2, z), a2, a);

      a = Tau(eta.x, eta.y, ind_s);
      u = Kappa(s, a, h);
      z = Chi_s.f(u);
      n_z = vector3_norm(Chi_s.n_f(u));
      Phi_times_Phi(
          c, w * Q->w[j] * n_z * n_y1 * LaplaceSingleLayerKernel(z, y1), a, b1);
      Phi_times_Phi(
          c, w * Q->w[j] * n_z * n_y2 * LaplaceSingleLayerKernel(z, y2), a, b2);
    }
  }
  return;
}

/**
 *  \brief         Determines the quadrature grade for two elements with
 *                 distance dist to reach accuracy \f$h=2^{-\alpha M}\f$.
 *
 *  \param[in]     dist           Distance of the two elements.
 *  \param[in]     alpha          Measure for the quadrature accuracy.
 *  \param[in]     M              Discretization level.
 *
 *  \return        Quadrature grade.
 *
 *  Some theory about the quadrature grade: If we have polynomial ansatz
 *  functions of degree p and we are working with the double layer potential
 *  ansatz, then the quadrature degree g to reach accuracy h is
 *     \f[
 *     g\geq\frac{\log(h^{-2p-1}\mathrm{dist}^{p-2})}{2\log(2\rho(p)\chi)},
 *     \f]
 *  where \f$ \chi=\max(\mathrm{dist}/h,1) \f$.
 *  If \f$\chi\f$ behaves like 1, then we have almost singular integrals which
 *  have to be evaluated with a high quadrature grade. In the case where
 *  \f$\chi\f$ behaves like \f$h^{-1}\f$ we are in the far-field case and we
 *  can choose a constant quadrature grade with precomputed boundary values
 *  as initialised in init_randwerte(). The quadrature grade in the far-field
 *  is defined in the constant \link g_far \endlink.
 *
 *  For more information about the quadrature grade see Theorem 5.3.30. of
 *     Stefan A. Sauter, Christoph Schwab, "Boundary Element Methods",
 *     Springer Series in Computational Mathematics, Volume 39 2011.
 *
 */
static int LaplaceSingle_quadrature_order(double dist, double alpha, int M,
                                          int a_o) {
  int g;

  dist = (dist * (1 << M) < 1) ? -(M * log(2)) : log(dist);

  g = 0.5 * ((alpha + a_o - 1) * M * log(2) + (a_o - 2) * dist) /
      ((M + 2) * log(2) + dist);

  assert(g <= g_max);

  return g;
}

static void LaplaceSingle_interpol_kernel(double ***ker1, double ***ker2,
                                          const Spl::Patch &geom1,
                                          const Spl::Patch &geom2, vector2 p1,
                                          vector2 p2, int i, int j, int k,
                                          double kappa[2]) {
  vector3 f1;
  vector3 f2;
  vector3 n1;
  vector3 n2;

  f1 = geom1.f(p1);
  f2 = geom2.f(p2);
  n1 = geom1.n_f(p1);
  n2 = geom2.n_f(p2);

  ker1[0][0][k * i + j] =
      LaplaceSingleLayerKernel(f1, f2) * vector3_norm(n1) * vector3_norm(n2);

  return;
}

static void LaplaceSingle_insert_fmat(double **A1, double **A2, int i, int j,
                                      int nl, int a_bs, int a_bs2,
                                      double *dest) {
  int k;

  for (k = 0; k < a_bs; ++k)
    memcpy(&A1[0][(i * a_bs + k) * nl * a_bs + j * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
  return;
}

static void LaplaceSingle_insert_sym_fmat_diag(double **A, int i, int nl,
                                               int a_bs, int a_bs2,
                                               double *dest) {
  int k;

  for (k = 0; k < a_bs; ++k)
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + i * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
  return;
}

static void LaplaceSingle_insert_sym_fmat_offdiag(double **A, int i, int j,
                                                  int nl, int a_bs, int a_bs2,
                                                  double *dest) {
  int k;
  int l;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + j * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
    for (l = 0; l < a_bs; ++l)
      A[0][(j * a_bs + l) * nl * a_bs + i * a_bs + k] = dest[k * a_bs + l];
  }
  return;
}

static void LaplaceSingle_pot_eval(double *Pot, vector3 *R, int nr, double Qgw,
                                   vector3 y, vector3 n_y, double *rhol,
                                   double h, double *d, void *disc) {
  int a_bs = ((discretization *)disc)->a_bs;
  int j;
  double w;

  w = myddot(a_bs, rhol, d) / h * vector3_norm(n_y);
  for (j = 0; j < nr; ++j)
    Pot[j] += Qgw * w * LaplaceSingleLayerKernel(R[j], y);

  return;
}

static void LaplaceSingle_pot_normalize(double *Pot, int nr, double h) {
  int j;

  for (j = 0; j < nr; j++) Pot[j] *= (h * h);
  return;
}

/**
 *  \brief         Function collecting the pointers to the correct phi-functions
 *                 for the assembly of the moment matrices from a
 * discretization.
 *
 *  \param[in]     disc          Discretization, however, due to compiler
 *                               reasons, we have to put a void here. So be
 *                               Careful.
 *  \param[in]     n             Moment matrices to which cluster tree.
 *
 *  \return        Array of phi-functions. For each cluster tree, two
 *                 phi-function pointers are required. They are then tensorized
 *                 for the assembly of the moment matrices.
 *
 */
static void (**LaplaceSingle_Tmom_left_phi(void *voiddisc,
                                           int n))(double *, double, double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  assert(n == 0);

  phi = (void (**)(double *, double, double))calloc(
      2, sizeof(void (*)(double *, double, double)));
  phi[0] = disc->phi;
  phi[1] = disc->phi;

  return phi;
}

/**
 *  \brief         Function collecting the pointers to the correct phi-functions
 *                 for the assembly of the moment matrices from a
 * discretization.
 *
 *  \param[in]     disc          Discretization, however, due to compiler
 *                               reasons, we have to put a void here. So be
 *                               Careful.
 *  \param[in]     n             Moment matrices to which cluster tree.
 *
 *  \return        Array of phi-functions. For each cluster tree, two
 *                 phi-function pointers are required. They are then tensorized
 *                 for the assembly of the moment matrices.
 *
 */
static void (**LaplaceSingle_Tmom_right_phi(void *voiddisc,
                                            int n))(double *, double, double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  assert(n == 0);

  phi = (void (**)(double *, double, double))calloc(
      2, sizeof(void (*)(double *, double, double)));
  phi[0] = disc->phi;
  phi[1] = disc->phi;

  return phi;
}

static void LaplaceSingle_postproc(void *hmatfac, void *H) { return; }

Eigen::VectorXd LaplaceSingle::MatVecImplementation(
    void *Hv, const Eigen::VectorXd &x) const {
  ct_root *H = (ct_root *)Hv;
  Eigen::VectorXd y = Eigen::VectorXd::Zero(H->disc->na);
  H2l2_HtimesVsmallReal(H, (double *)x.data(), (double *)y.data());
  return y;
}

LaplaceSingle::LaplaceSingle(void) {
  _pde.pdetype = Laplace;
  _pde.name = "LaplaceSingle";
  _pde.nct = 1;
  _pde.np_max_fac = 1;
  _pde.Tmom_left_phi = &LaplaceSingle_Tmom_left_phi;
  _pde.Tmom_right_phi = &LaplaceSingle_Tmom_right_phi;
  _pde.symflags = "S";
  _pde.quadrature_accuracy = &LaplaceSingle_quadrature_accuracy;
  _pde.quadrature_bufsize = 1;
  _pde.g_far = &LaplaceSingle_g_far;
  _pde.g_pot = &LaplaceSingle_g_pot;
  _pde.init_randwerte = &LaplaceSingle_init_randwerte;
  _pde.free_randwerte = &LaplaceSingle_free_randwerte;
  _pde.IntPhi0 = &LaplaceSingle_IntPhi0;
  _pde.IntPhi1 = &LaplaceSingle_IntPhi1;
  _pde.IntPhi2 = &LaplaceSingle_IntPhi2;
  _pde.IntPhi3 = &LaplaceSingle_IntPhi3;
  _pde.IntPhi4 = &LaplaceSingle_IntPhi4;
  _pde.quadrature_order = &LaplaceSingle_quadrature_order;
  _pde.interpol_kernel = &LaplaceSingle_interpol_kernel;
  _pde.insert_fmat = &LaplaceSingle_insert_fmat;
  _pde.insert_sym_fmat_diag = &LaplaceSingle_insert_sym_fmat_diag;
  _pde.insert_sym_fmat_offdiag = &LaplaceSingle_insert_sym_fmat_offdiag;
  _pde.pot_eval = &LaplaceSingle_pot_eval;
  _pde.pot_normalize = &LaplaceSingle_pot_normalize;
  _pde.postproc = &LaplaceSingle_postproc;
};
}  // namespace Bembel
