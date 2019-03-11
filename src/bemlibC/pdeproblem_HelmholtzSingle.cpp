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

static inline void HelmholtzSingleKernel(double kappa[2], double d[2],
                                         vector3 x, vector3 y) {
  double r, a, b, e;
  vector3 c;
  c.x = x.x - y.x;
  c.y = x.y - y.y;
  c.z = x.z - y.z;
  r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
  e = exp(-kappa[1] * r);
  a = e * cos(kappa[0] * r);
  b = e * sin(kappa[0] * r);
  r *= 4 * M_PI;
  d[0] = a / r; /* real(Einfachschicht) */
  d[1] = b / r; /* imag(Einfachschicht) */
  return;
}

static inline int HelmholtzSingle_quadrature_accuracy(int a_o) {
  return 3 + 2 * (a_o - 1);
}

static inline int HelmholtzSingle_g_far(int a_o) { return a_o; }

static inline int HelmholtzSingle_g_pot(int a_o) { return a_o; }

static void HelmholtzSingle_init_randwerte(
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
static void HelmholtzSingle_free_randwerte(
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
static void HelmholtzSingle_IntPhi0(double *c, int no1, int no2, cubature *Q,
                                    randwerte **RW, void *disc) {
  int i, j;    /* Laufindizes */
  double w;    /* Quadraturgewicht */
  double d[2]; /* Ergebnis der Kernauswertung */
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;
  memset(c, 0, 2 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    for (j = 0; j < Q->nop; j++) {
      w = RW[no1][i].det_dChi * RW[no2][j].det_dChi;
      HelmholtzSingleKernel(kappa, d, RW[no1][i].Chi, RW[no2][j].Chi);
      Phi_times_Phi(&c[0], w * d[0], Q->xi[i], Q->xi[j]);
      Phi_times_Phi(&c[a_bs2], w * d[1], Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * No-Problem-Quadrature-Routine -> kanonisches Skalarprodukt
 *
 */
static void HelmholtzSingle_IntPhi1(double *c, vector2 s, vector2 t, double h,
                                    const Spl::Patch &Chi_s,
                                    const Spl::Patch &Chi_t, cubature *Q,
                                    void *disc) {
  int i, j;
  double d[2], w1, w2;
  vector2 xi, eta;
  vector3 x;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, 2 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = vector2_add(s, vector2_Smul(h, Q->xi[i]));
    x = Chi_s.f(xi);
    w1 = h * h * Q->w[i] * vector3_norm(Chi_s.n_f(xi));
    for (j = 0; j < Q->nop; j++) {
      eta = vector2_add(t, vector2_Smul(h, Q->xi[j]));
      w2 = Q->w[j] * vector3_norm(Chi_t.n_f(eta)) * w1;
      HelmholtzSingleKernel(kappa, d, x, Chi_t.f(eta));
      Phi_times_Phi(&c[0], w2 * d[0], Q->xi[i], Q->xi[j]);
      Phi_times_Phi(&c[a_bs2], w2 * d[1], Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * GLEICHE PATCHES [0,1]^2 -> kanonisches Skalarprodukt
 *
 */
static void HelmholtzSingle_IntPhi2(double *c, vector2 s, double h,
                                    const Spl::Patch &Chi_s, cubature *Q,
                                    void *disc) {
  int i, j;
  double d[2], r, w;
  double t1, t2, t3, t4;
  vector2 xi, eta, a, b;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, 2 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    r = h * h * Q->w[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
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
      w = r * Q->w[j] * vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(a), Chi_s.f(b));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t1, t2),
                    vector2_make(t3, t4));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t1, t2),
                    vector2_make(t3, t4));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t3, t4),
                    vector2_make(t1, t2));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t3, t4),
                    vector2_make(t1, t2));

      a.y = s.y + h * t4;
      b.y = s.y + h * t2;
      w = r * Q->w[j] * vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(a), Chi_s.f(b));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t1, t4),
                    vector2_make(t3, t2));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t1, t4),
                    vector2_make(t3, t2));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t3, t2),
                    vector2_make(t1, t4));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t3, t2),
                    vector2_make(t1, t4));

      a.x = s.x + h * t2;
      a.y = s.y + h * t1;
      b.x = s.x + h * t4;
      b.y = s.y + h * t3;
      w = r * Q->w[j] * vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(a), Chi_s.f(b));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t2, t1),
                    vector2_make(t4, t3));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t2, t1),
                    vector2_make(t4, t3));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t4, t3),
                    vector2_make(t2, t1));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t4, t3),
                    vector2_make(t2, t1));

      a.y = s.y + h * t3;
      b.y = s.y + h * t1;
      w = r * Q->w[j] * vector3_norm(Chi_s.n_f(a)) * vector3_norm(Chi_s.n_f(b));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(a), Chi_s.f(b));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t2, t3),
                    vector2_make(t4, t1));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t2, t3),
                    vector2_make(t4, t1));
      Phi_times_Phi(&c[0], w * d[0], vector2_make(t4, t1),
                    vector2_make(t2, t3));
      Phi_times_Phi(&c[a_bs2], w * d[1], vector2_make(t4, t1),
                    vector2_make(t2, t3));
    }
  }

  return;
}

/**
 * GEMEINSAME KANTE [0,1] -> kanonisches Skalarprodukt
 *
 */
static void HelmholtzSingle_IntPhi3(double *c, vector2 s, vector2 t, double h,
                                    int ind_s, int ind_t,
                                    const Spl::Patch &Chi_s,
                                    const Spl::Patch &Chi_t, cubature *Q,
                                    void *disc) {
  int i, j;
  double d[2], w, r, t1, t2, t3, t4;
  vector2 xi, eta, a, b, u, v;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;
  memset(c, 0, 2 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    r = h * h * xi.y * xi.y * Q->w[i];
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
      w = r * Q->w[j] * (1 - xi.y) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);

      a = Tau(1 - t1, eta.x, ind_s);
      b = Tau(1 - t2, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - xi.y) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);

      a = Tau(t3, xi.y, ind_s);
      b = Tau(t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);

      a = Tau(1 - t3, xi.y, ind_s);
      b = Tau(1 - t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);

      a = Tau(t4, eta.y, ind_s);
      b = Tau(t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);

      a = Tau(1 - t4, eta.y, ind_s);
      b = Tau(1 - t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x) * vector3_norm(Chi_s.n_f(u)) *
          vector3_norm(Chi_t.n_f(v));
      HelmholtzSingleKernel(kappa, d, Chi_s.f(u), Chi_t.f(v));
      Phi_times_Phi(&c[0], w * d[0], a, b);
      Phi_times_Phi(&c[a_bs2], w * d[1], a, b);
    }
  }

  return;
}

/**
 * GEMEINSAME ECKE IM NULLPUNKT -> kanonisches Skalarprodukt
 *
 */
static void HelmholtzSingle_IntPhi4(double *c, vector2 s, vector2 t, double h,
                                    int ind_s, int ind_t,
                                    const Spl::Patch &Chi_s,
                                    const Spl::Patch &Chi_t, cubature *Q,
                                    void *disc) {
  int i, j;
  double d[2], w1, w2, norm_x1, norm_x2, norm_y1, norm_y2;
  vector2 xi, eta, a, u, a1, a2, b1, b2;
  vector3 x1, x2, y1, y2, z;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  void (*Phi_times_Phi)(double *, double, vector2, vector2) =
      mydisc->Phi_times_Phi;

  memset(c, 0, 2 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    w1 = h * h * pow(xi.x, 3) * Q->w[i];
    xi.y *= xi.x;
    a1 = Tau(xi.x, xi.y, ind_s);
    a2 = Tau(xi.y, xi.x, ind_s);
    b1 = Tau(xi.x, xi.y, ind_t);
    b2 = Tau(xi.y, xi.x, ind_t);

    u = Kappa(s, a1, h);
    x1 = Chi_s.f(u);
    norm_x1 = vector3_norm(Chi_s.n_f(u));

    u = Kappa(s, a2, h);
    x2 = Chi_s.f(u);
    norm_x2 = vector3_norm(Chi_s.n_f(u));

    u = Kappa(t, b1, h);
    y1 = Chi_t.f(u);
    norm_y1 = vector3_norm(Chi_t.n_f(u));

    u = Kappa(t, b2, h);
    y2 = Chi_t.f(u);
    norm_y2 = vector3_norm(Chi_t.n_f(u));

    for (j = 0; j < Q->nop; j++) {
      eta = vector2_Smul(xi.x, Q->xi[j]);

      a = Tau(eta.x, eta.y, ind_t);
      u = Kappa(t, a, h);
      z = Chi_t.f(u);

      w2 = w1 * Q->w[j] * vector3_norm(Chi_t.n_f(u));
      HelmholtzSingleKernel(kappa, d, x1, z);
      Phi_times_Phi(&c[0], w2 * norm_x1 * d[0], a1, a);
      Phi_times_Phi(&c[a_bs2], w2 * norm_x1 * d[1], a1, a);
      HelmholtzSingleKernel(kappa, d, x2, z);
      Phi_times_Phi(&c[0], w2 * norm_x2 * d[0], a2, a);
      Phi_times_Phi(&c[a_bs2], w2 * norm_x2 * d[1], a2, a);

      a = Tau(eta.x, eta.y, ind_s);
      u = Kappa(s, a, h);
      z = Chi_s.f(u);

      w2 = w1 * Q->w[j] * vector3_norm(Chi_s.n_f(u));
      HelmholtzSingleKernel(kappa, d, z, y1);
      Phi_times_Phi(&c[0], w2 * norm_y1 * d[0], a, b1);
      Phi_times_Phi(&c[a_bs2], w2 * norm_y1 * d[1], a, b1);
      HelmholtzSingleKernel(kappa, d, z, y2);
      Phi_times_Phi(&c[0], w2 * norm_y2 * d[0], a, b2);
      Phi_times_Phi(&c[a_bs2], w2 * norm_y2 * d[1], a, b2);
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
static int HelmholtzSingle_quadrature_order(double dist, double alpha, int M,
                                            int a_o) {
  int g;

  dist = (dist * (1 << M) < 1) ? -(M * log(2)) : log(dist);

  g = 0.5 * ((alpha + a_o - 1) * M * log(2) + (a_o - 2) * dist) /
      ((M + 2) * log(2) + dist);

  assert(g <= g_max);

  return g;
}

static void HelmholtzSingle_interpol_kernel(double ***ker1, double ***ker2,
                                            const Spl::Patch &geom1,
                                            const Spl::Patch &geom2, vector2 p1,
                                            vector2 p2, int i, int j, int k,
                                            double kappa[2]) {
  double d[2];
  double norm;
  vector3 f1;
  vector3 f2;
  vector3 n1;
  vector3 n2;

  f1 = geom1.f(p1);
  f2 = geom2.f(p2);
  n1 = geom1.n_f(p1);
  n2 = geom2.n_f(p2);

  norm = vector3_norm(n1) * vector3_norm(n2);
  HelmholtzSingleKernel(kappa, d, f1, f2);
  ker1[0][0][k * i + j] = d[0] * norm;
  ker1[1][0][k * i + j] = d[1] * norm;

  return;
}

static void HelmholtzSingle_insert_fmat(double **A1, double **A2, int i, int j,
                                        int nl, int a_bs, int a_bs2,
                                        double *dest) {
  int k;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A1[0][(i * a_bs + k) * nl * a_bs + j * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
    memcpy(&A1[1][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + a_bs2, a_bs * sizeof(double));
  }
  return;
}

static void HelmholtzSingle_insert_sym_fmat_diag(double **A, int i, int nl,
                                                 int a_bs, int a_bs2,
                                                 double *dest) {
  int k;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + i * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
    memcpy(&A[1][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + a_bs2, a_bs * sizeof(double));
  }
  return;
}

static void HelmholtzSingle_insert_sym_fmat_offdiag(double **A, int i, int j,
                                                    int nl, int a_bs, int a_bs2,
                                                    double *dest) {
  int k;
  int l;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + j * a_bs], dest + a_bs * k,
           a_bs * sizeof(double));
    memcpy(&A[1][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + a_bs2, a_bs * sizeof(double));
    for (l = 0; l < a_bs; ++l) {
      A[0][(j * a_bs + l) * nl * a_bs + i * a_bs + k] = dest[k * a_bs + l];
      A[1][(j * a_bs + l) * nl * a_bs + i * a_bs + k] =
          dest[a_bs2 + k * a_bs + l];
    }
  }
  return;
}

static void HelmholtzSingle_pot_eval(double *Pot, vector3 *R, int nr,
                                     double Qgw, vector3 y, vector3 n_y,
                                     double *rhol, double h, double *d,
                                     void *disc) {
  int j;
  int a_bs = ((discretization *)disc)->a_bs;
  int N = a_bs * ((discretization *)disc)->mesh->nf;
  double c[2];
  double w[2];
  double *kappa = ((discretization *)disc)->pde->kappa;

  w[0] = Qgw * myddot(a_bs, rhol, d);
  w[1] = Qgw * myddot(a_bs, rhol + N, d);
  for (j = 0; j < nr; ++j) {
    HelmholtzSingleKernel(kappa, c, R[j], y);
    c[0] *= vector3_norm(n_y);
    c[1] *= vector3_norm(n_y);
    Pot[j] += w[0] * c[0] - w[1] * c[1];
    Pot[j + nr] += w[0] * c[1] + w[1] * c[0];
  }

  return;
}

static void HelmholtzSingle_pot_normalize(double *Pot, int nr, double h) {
  int j;

  for (j = 0; j < 2 * nr; j++) Pot[j] *= h;
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
static void (**HelmholtzSingle_Tmom_left_phi(void *voiddisc,
                                             int n))(double *, double, double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  phi = (void (**)(double *, double, double))calloc(
      4, sizeof(void (*)(double *, double, double)));
  switch (n) {
    case 0:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      break;
    case 1:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      break;
    default:
      assert(!"this should not happen");
  }

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
static void (**HelmholtzSingle_Tmom_right_phi(void *voiddisc, int n))(double *,
                                                                      double,
                                                                      double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  phi = (void (**)(double *, double, double))calloc(
      4, sizeof(void (*)(double *, double, double)));
  switch (n) {
    case 0:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      break;
    case 1:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      break;
    default:
      assert(!"this should not happen");
  }

  return phi;
}

static void HelmholtzSingle_postproc(void *hmatfac, void *H) { return; }

Eigen::VectorXcd HelmholtzSingle::MatVecImplementation(
    void *Hv, const Eigen::VectorXcd &x) const {
  ct_root *H = (ct_root *)Hv;
  int na = H->disc->na / 2;
  Eigen::VectorXd xr = Eigen::VectorXd(2 * na);
  Eigen::VectorXd yr = Eigen::VectorXd::Zero(2 * na);
  Eigen::VectorXcd y = Eigen::VectorXcd(na);
  xr << x.real(), x.imag();
  H2l2_HtimesVsmallComplex(H, (double *)xr.data(), (double *)yr.data());
  y.real() = yr.head(na);
  y.imag() = yr.tail(na);
  return y;
}

HelmholtzSingle::HelmholtzSingle(std::complex<double> in) {
  _pde.kappa[0] = in.real();
  _pde.kappa[1] = in.imag();
  _pde.pdetype = Helmholtz;
  _pde.name = "HelmholtzSingle";
  _pde.nct = 2;
  _pde.np_max_fac = 1;
  _pde.Tmom_left_phi = &HelmholtzSingle_Tmom_left_phi;
  _pde.Tmom_right_phi = &HelmholtzSingle_Tmom_right_phi;
  _pde.symflags = "SS";
  _pde.quadrature_accuracy = &HelmholtzSingle_quadrature_accuracy;
  _pde.quadrature_bufsize = 2;
  _pde.g_far = &HelmholtzSingle_g_far;
  _pde.g_pot = &HelmholtzSingle_g_pot;
  _pde.init_randwerte = &HelmholtzSingle_init_randwerte;
  _pde.free_randwerte = &HelmholtzSingle_free_randwerte;
  _pde.IntPhi0 = &HelmholtzSingle_IntPhi0;
  _pde.IntPhi1 = &HelmholtzSingle_IntPhi1;
  _pde.IntPhi2 = &HelmholtzSingle_IntPhi2;
  _pde.IntPhi3 = &HelmholtzSingle_IntPhi3;
  _pde.IntPhi4 = &HelmholtzSingle_IntPhi4;
  _pde.quadrature_order = &HelmholtzSingle_quadrature_order;
  _pde.interpol_kernel = &HelmholtzSingle_interpol_kernel;
  _pde.insert_fmat = &HelmholtzSingle_insert_fmat;
  _pde.insert_sym_fmat_diag = &HelmholtzSingle_insert_sym_fmat_diag;
  _pde.insert_sym_fmat_offdiag = &HelmholtzSingle_insert_sym_fmat_offdiag;
  _pde.pot_eval = &HelmholtzSingle_pot_eval;
  _pde.pot_normalize = &HelmholtzSingle_pot_normalize;
  _pde.postproc = &HelmholtzSingle_postproc;
};
}  // namespace Bembel
