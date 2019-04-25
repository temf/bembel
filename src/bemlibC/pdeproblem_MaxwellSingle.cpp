// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <array>
#include "H2_level2.h"
#include "cluster_tree.h"
#include "discretization.h"
#include "gram.h"
#include "pdeproblem.h"

namespace Bembel {

static inline void kernel(double kappa[2], double d[2], vector3 x, vector3 y) {
  double r, a, b, e;
  vector3 c;
  c.x = x.x - y.x;
  c.y = x.y - y.y;
  c.z = x.z - y.z;
  r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
  e = exp(kappa[1] * r);
  a = e * cos(-kappa[0] * r);
  b = e * sin(-kappa[0] * r);
  r *= 4 * M_PI;
  d[0] = a / r; /* real(Einfachschicht) */
  d[1] = b / r; /* imag(Einfachschicht) */
  return;
}

static inline void grad_kernel(double kappa[2], vector3 d[2], vector3 x,
                               vector3 y) {
  double r;
  double r2;
  double f;
  double e[2];
  vector3 c;
  c.x = x.x - y.x;
  c.y = x.y - y.y;
  c.z = x.z - y.z;
  r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
  r2 = r * r;
  f = exp(kappa[1] * r);
  e[0] = f * cos(-kappa[0] * r) / r2;
  e[1] = f * sin(-kappa[0] * r) / r2;
  d[0] = vector3_Smul(
      (kappa[0] * e[1] + kappa[1] * e[0] - e[0] / r) / (4 * M_PI), c);
  d[1] = vector3_Smul(
      (-kappa[0] * e[0] + kappa[1] * e[1] - e[1] / r) / (4 * M_PI), c);
  return;
}

static int MaxwellSingle_quadrature_accuracy(int a_o) {
  return 3 + 2 * (a_o - 1);
}

static inline int MaxwellSingle_g_far(int a_o) { return a_o; }

static inline int MaxwellSingle_g_pot(int a_o) { return a_o; }

static void MaxwellSingle_init_randwerte(
    randwerte ***RW,     /* zu berechnende * Randwerte *
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

  (*RW) = (randwerte **)calloc(p * n * n, sizeof(randwerte *));
  zi = 0;
  for (i1 = 0; i1 < p; i1++) {
    s.y = 0;
    for (i2 = 0; i2 < n; i2++) {
      s.x = 0;
      for (i3 = 0; i3 < n; i3++) {
        (*RW)[zi] = (randwerte *)calloc(Q->nop, sizeof(randwerte));
        for (k = 0; k < Q->nop; k++) {
          t = vector2_add(s, vector2_Smul(h, Q->xi[k]));
          (*RW)[zi][k].Chi = Chi[i1].f(t);
          const std::array<vector3, 2> df = Chi[i1].jacobian(t);
          (*RW)[zi][k].dx_Chi = vector3_Smul(Q->w[k], df[0]);
          (*RW)[zi][k].dy_Chi = vector3_Smul(Q->w[k], df[1]);
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
static void MaxwellSingle_free_randwerte(
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
static void MaxwellSingle_IntPhi0(double *c, int no1, int no2, cubature *Q,
                                  randwerte **RW, void *disc) {
  int i, j;    /* Laufindizes */
  double d[2]; /* Ergebnis der Kernauswertung */
  double phiw;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  int n = 1 << mydisc->mesh->M; /* n*n element per parameter
                                 * domain */
  double h = 1. / n;
  double *kappa = mydisc->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3) = mydisc->VPhi_scal_VPhi;
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2) =
      mydisc->Div_Phi_times_Div_Phi;
  memset(c, 0, 8 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    for (j = 0; j < Q->nop; j++) {
      kernel(kappa, d, RW[no1][i].Chi, RW[no2][j].Chi);
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], Q->xi[i], Q->xi[j], RW[no1][i].dx_Chi,
                     RW[no1][i].dy_Chi, RW[no2][j].dx_Chi, RW[no2][j].dy_Chi);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], Q->xi[i], Q->xi[j],
                     RW[no1][i].dx_Chi, RW[no1][i].dy_Chi, RW[no2][j].dx_Chi,
                     RW[no2][j].dy_Chi);
      phiw =
          -Q->w[i] * Q->w[j] / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            Q->xi[i], Q->xi[j]);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * No-Problem-Quadrature-Routine -> kanonisches Skalarprodukt
 *
 */
static void MaxwellSingle_IntPhi1(double *c, vector2 s, vector2 t, double h,
                                  const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double d[2], w1, w2;
  double phiw;
  vector2 xi, eta;
  vector3 x;
  vector3 y;
  std::array<vector3, 2> dx;
  std::array<vector3, 2> dy;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3) = mydisc->VPhi_scal_VPhi;
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2) =
      mydisc->Div_Phi_times_Div_Phi;

  memset(c, 0, 8 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = vector2_add(s, vector2_Smul(h, Q->xi[i]));
    x = Chi_s.f(xi);
    w1 = Q->w[i];
    dx = Chi_s.jacobian(xi);
    for (j = 0; j < Q->nop; j++) {
      eta = vector2_add(t, vector2_Smul(h, Q->xi[j]));
      y = Chi_t.f(eta);
      w2 = Q->w[j] * w1;
      dy = Chi_t.jacobian(eta);
      kernel(kappa, d, x, y);
      d[0] *= w2;
      d[1] *= w2;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], Q->xi[i], Q->xi[j], dx[0], dx[1],
                     dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], Q->xi[i], Q->xi[j], dx[0],
                     dx[1], dy[0], dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            Q->xi[i], Q->xi[j]);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            Q->xi[i], Q->xi[j]);
    }
  }
  return;
}

/**
 * GLEICHE PATCHES [0,1]^2 -> kanonisches Skalarprodukt
 *
 */
static void MaxwellSingle_IntPhi2(double *c, vector2 s, double h,
                                  const Spl::Patch &Chi_s, cubature *Q,
                                  void *disc) {
  int i, j;
  double phiw;
  double d[2], r, w;
  double t1, t2, t3, t4;
  vector2 xi, eta, a, b;
  vector3 x;
  vector3 y;
  std::array<vector3, 2> dx;
  std::array<vector3, 2> dy;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3) = mydisc->VPhi_scal_VPhi;
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2) =
      mydisc->Div_Phi_times_Div_Phi;

  memset(c, 0, 8 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    r = Q->w[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
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
      w = r * Q->w[j];
      x = Chi_s.f(a);
      dx = Chi_s.jacobian(a);
      y = Chi_s.f(b);
      dy = Chi_s.jacobian(b);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t1, t2),
                     vector2_make(t3, t4), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t1, t2),
                     vector2_make(t3, t4), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t3, t4),
                     vector2_make(t1, t2), dy[0], dy[1], dx[0], dx[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t3, t4),
                     vector2_make(t1, t2), dy[0], dy[1], dx[0], dx[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t1, t2), vector2_make(t3, t4));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t1, t2), vector2_make(t3, t4));
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t3, t4), vector2_make(t1, t2));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t3, t4), vector2_make(t1, t2));
      a.y = s.y + h * t4;
      b.y = s.y + h * t2;
      w = r * Q->w[j];
      x = Chi_s.f(a);
      dx = Chi_s.jacobian(a);
      y = Chi_s.f(b);
      dy = Chi_s.jacobian(b);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t1, t4),
                     vector2_make(t3, t2), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t1, t4),
                     vector2_make(t3, t2), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t3, t2),
                     vector2_make(t1, t4), dy[0], dy[1], dx[0], dx[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t3, t2),
                     vector2_make(t1, t4), dy[0], dy[1], dx[0], dx[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t1, t4), vector2_make(t3, t2));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t1, t4), vector2_make(t3, t2));
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t3, t2), vector2_make(t1, t4));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t3, t2), vector2_make(t1, t4));

      a.x = s.x + h * t2;
      a.y = s.y + h * t1;
      b.x = s.x + h * t4;
      b.y = s.y + h * t3;
      w = r * Q->w[j];
      x = Chi_s.f(a);
      dx = Chi_s.jacobian(a);
      y = Chi_s.f(b);
      dy = Chi_s.jacobian(b);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t2, t1),
                     vector2_make(t4, t3), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t2, t1),
                     vector2_make(t4, t3), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t4, t3),
                     vector2_make(t2, t1), dy[0], dy[1], dx[0], dx[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t4, t3),
                     vector2_make(t2, t1), dy[0], dy[1], dx[0], dx[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t2, t1), vector2_make(t4, t3));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t2, t1), vector2_make(t4, t3));
      Div_Phi_times_Div_Phi(&c[0 * a_bs2],
                            phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t4, t3), vector2_make(t2, t1));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t4, t3), vector2_make(t2, t1));

      a.y = s.y + h * t3;
      b.y = s.y + h * t1;
      w = r * Q->w[j];
      x = Chi_s.f(a);
      dx = Chi_s.jacobian(a);
      y = Chi_s.f(b);
      dy = Chi_s.jacobian(b);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t2, t3),
                     vector2_make(t4, t1), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t2, t3),
                     vector2_make(t4, t1), dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[0], phiw * d[0], vector2_make(t4, t1),
                     vector2_make(t2, t3), dy[0], dy[1], dx[0], dx[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], vector2_make(t4, t1),
                     vector2_make(t2, t3), dy[0], dy[1], dx[0], dx[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t2, t3), vector2_make(t4, t1));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t2, t3), vector2_make(t4, t1));
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            vector2_make(t4, t1), vector2_make(t2, t3));
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]),
                            vector2_make(t4, t1), vector2_make(t2, t3));
    }
  }

  return;
}

/**
 * GEMEINSAME KANTE [0,1] -> kanonisches Skalarprodukt
 *
 */
static void MaxwellSingle_IntPhi3(double *c, vector2 s, vector2 t, double h,
                                  int ind_s, int ind_t, const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double phiw;
  double d[2], w, r, t1, t2, t3, t4;
  vector2 xi, eta, a, b, u, v;
  vector3 x;
  std::array<vector3, 2> dx;
  vector3 y;
  std::array<vector3, 2> dy;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3) = mydisc->VPhi_scal_VPhi;
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2) =
      mydisc->Div_Phi_times_Div_Phi;
  memset(c, 0, 8 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    r = xi.y * xi.y * Q->w[i];
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
      w = r * Q->w[j] * (1 - xi.y);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);

      a = Tau(1 - t1, eta.x, ind_s);
      b = Tau(1 - t2, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - xi.y);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0 * a_bs2],
                            phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]), a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);

      a = Tau(t3, xi.y, ind_s);
      b = Tau(t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);

      a = Tau(1 - t3, xi.y, ind_s);
      b = Tau(1 - t4, eta.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);

      a = Tau(t4, eta.y, ind_s);
      b = Tau(t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);

      a = Tau(1 - t4, eta.y, ind_s);
      b = Tau(1 - t3, xi.y, ind_t);
      u = Kappa(s, a, h);
      v = Kappa(t, b, h);
      w = r * Q->w[j] * (1 - eta.x);
      x = Chi_s.f(u);
      dx = Chi_s.jacobian(u);
      y = Chi_t.f(v);
      dy = Chi_t.jacobian(v);
      kernel(kappa, d, x, y);
      d[0] *= w;
      d[1] *= w;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b, dx[0], dx[1], dy[0], dy[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b, dx[0], dx[1], dy[0],
                     dy[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b);
      Div_Phi_times_Div_Phi(&c[4 * a_bs2],
                            phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b);
    }
  }

  return;
}

/**
 * GEMEINSAME ECKE IM NULLPUNKT -> kanonisches Skalarprodukt
 *
 */
static void MaxwellSingle_IntPhi4(double *c, vector2 s, vector2 t, double h,
                                  int ind_s, int ind_t, const Spl::Patch &Chi_s,
                                  const Spl::Patch &Chi_t, cubature *Q,
                                  void *disc) {
  int i, j;
  double phiw;
  double d[2], w1, w2;
  vector2 xi, eta, a, u, a1, a2, b1, b2;
  vector3 x1, x2, y1, y2, z;
  std::array<vector3, 2> dx1, dx2, dy1, dy2, dz;
  discretization *mydisc = (discretization *)disc;
  int a_bs2 = mydisc->a_bs2;
  double *kappa = mydisc->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3) = mydisc->VPhi_scal_VPhi;
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2) =
      mydisc->Div_Phi_times_Div_Phi;
  memset(c, 0, 8 * a_bs2 * sizeof(double));
  for (i = 0; i < Q->nop; i++) {
    xi = Q->xi[i];
    w1 = pow(xi.x, 3) * Q->w[i];
    xi.y *= xi.x;
    a1 = Tau(xi.x, xi.y, ind_s);
    a2 = Tau(xi.y, xi.x, ind_s);
    b1 = Tau(xi.x, xi.y, ind_t);
    b2 = Tau(xi.y, xi.x, ind_t);

    u = Kappa(s, a1, h);
    x1 = Chi_s.f(u);
    dx1 = Chi_s.jacobian(u);
    u = Kappa(s, a2, h);
    x2 = Chi_s.f(u);
    dx2 = Chi_s.jacobian(u);
    u = Kappa(t, b1, h);
    y1 = Chi_t.f(u);
    dy1 = Chi_t.jacobian(u);
    u = Kappa(t, b2, h);
    y2 = Chi_t.f(u);
    dy2 = Chi_t.jacobian(u);

    for (j = 0; j < Q->nop; j++) {
      eta = vector2_Smul(xi.x, Q->xi[j]);

      a = Tau(eta.x, eta.y, ind_t);
      u = Kappa(t, a, h);
      z = Chi_t.f(u);
      w2 = w1 * Q->w[j];
      dz = Chi_t.jacobian(u);
      kernel(kappa, d, x1, z);
      d[0] *= w2;
      d[1] *= w2;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a1, a, dx1[0], dx1[1], dz[0], dz[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a1, a, dx1[0], dx1[1], dz[0],
                     dz[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a1, a);
      Div_Phi_times_Div_Phi(
          &c[4 * a_bs2], phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a1, a);
      kernel(kappa, d, x2, z);
      d[0] *= w2;
      d[1] *= w2;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a2, a, dx2[0], dx2[1], dz[0], dz[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a2, a, dx2[0], dx2[1], dz[0],
                     dz[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a2, a);
      Div_Phi_times_Div_Phi(
          &c[4 * a_bs2], phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a2, a);
      a = Tau(eta.x, eta.y, ind_s);
      u = Kappa(s, a, h);
      z = Chi_s.f(u);
      w2 = w1 * Q->w[j];
      dz = Chi_s.jacobian(u);
      kernel(kappa, d, z, y1);
      d[0] *= w2;
      d[1] *= w2;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b1, dz[0], dz[1], dy1[0], dy1[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b1, dz[0], dz[1], dy1[0],
                     dy1[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b1);
      Div_Phi_times_Div_Phi(
          &c[4 * a_bs2], phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b1);
      kernel(kappa, d, z, y2);
      d[0] *= w2;
      d[1] *= w2;
      phiw = h * h;
      VPhi_scal_VPhi(&c[0], phiw * d[0], a, b2, dz[0], dz[1], dy2[0], dy2[1]);
      VPhi_scal_VPhi(&c[4 * a_bs2], phiw * d[1], a, b2, dz[0], dz[1], dy2[0],
                     dy2[1]);
      phiw = -1. / (kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1]);
      Div_Phi_times_Div_Phi(&c[0], phiw * (kappa2[0] * d[0] + kappa2[1] * d[1]),
                            a, b2);
      Div_Phi_times_Div_Phi(
          &c[4 * a_bs2], phiw * (kappa2[0] * d[1] - kappa2[1] * d[0]), a, b2);
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
static int MaxwellSingle_quadrature_order(double dist, double alpha, int M,
                                          int a_o) {
  int g;

  dist = (dist * (1 << M) < 1) ? -(M * log(2)) : log(dist);

  g = 0.5 * ((alpha + a_o - 1) * M * log(2) + (a_o - 2) * dist) /
      ((M + 2) * log(2) + dist);

  assert(g <= g_max);

  return g;
}

static void MaxwellSingle_interpol_kernel(double ***ker1, double ***ker2,
                                          const Spl::Patch &geom1,
                                          const Spl::Patch &geom2, vector2 p1,
                                          vector2 p2, int i, int j, int k,
                                          double kappa[2]) {
  int np2 = k / 2;
  double w;
  double d[2];
  double dh[2];
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  double kappa2abs2 = kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1];
  vector3 f1;
  vector3 f2;
  vector3 d1_dx;
  vector3 d1_dy;
  vector3 d2_dx;
  vector3 d2_dy;

  f1 = geom1.f(p1);
  f2 = geom2.f(p2);
  d1_dx = geom1.df_dx(p1);
  d1_dy = geom1.df_dy(p1);
  d2_dx = geom2.df_dx(p2);
  d2_dy = geom2.df_dy(p2);
  kernel(kappa, d, f1, f2);

  /* first part: scalar product of ansatz functions */
  w = vector3_skalp(d1_dx, d2_dx);
  ker1[0][0][k * i + j] = d[0] * w; /* lower matrix */
  ker1[3][0][k * i + j] = d[1] * w; /* lower matrix */
  w = vector3_skalp(d1_dx, d2_dy);
  ker1[1][0][k * i + j] = d[0] * w /* lower matrix */;
  ker1[4][0][k * i + j] = d[1] * w; /* lower matrix */
  w = vector3_skalp(d1_dy, d2_dx);
  ker2[1][0][k * j + i] = d[0] * w; /* upper matrix, different indices */
  ker2[4][0][k * j + i] = d[1] * w; /* upper matrix, different indices */
  w = vector3_skalp(d1_dy, d2_dy);
  ker1[2][0][k * i + j] = d[0] * w; /* lower matrix */
  ker1[5][0][k * i + j] = d[1] * w; /* lower matrix */

  /* second part: surface divergence part of ansatz functions */
  dh[0] = d[0];
  dh[1] = d[1];
  d[0] = -(dh[0] * kappa2[0] + dh[1] * kappa2[1]) / kappa2abs2;
  d[1] = -(dh[1] * kappa2[0] - dh[0] * kappa2[1]) / kappa2abs2;
  ker1[0][0][k * (i + np2) + j + np2] = d[0]; /* lower matrix */
  ker1[1][0][k * (i + np2) + j + np2] = d[0]; /* lower matrix */
  ker1[2][0][k * (i + np2) + j + np2] = d[0]; /* lower matrix */
  ker1[3][0][k * (i + np2) + j + np2] = d[1]; /* lower matrix */
  ker1[4][0][k * (i + np2) + j + np2] = d[1]; /* lower matrix */
  ker1[5][0][k * (i + np2) + j + np2] = d[1]; /* lower matrix */
  ker2[1][0][k * (j + np2) + i + np2] =
      d[0]; /* upper matrix, different indices */
  ker2[4][0][k * (j + np2) + i + np2] =
      d[1]; /* upper matrix, different indices */

  ker1[0][0][k * (i + np2) + j] = 0;
  ker1[1][0][k * (i + np2) + j] = 0;
  ker2[1][0][k * (i + np2) + j] = 0;
  ker1[2][0][k * (i + np2) + j] = 0;
  ker1[3][0][k * (i + np2) + j] = 0;
  ker1[4][0][k * (i + np2) + j] = 0;
  ker1[5][0][k * (i + np2) + j] = 0;
  ker2[4][0][k * (i + np2) + j] = 0;

  ker1[0][0][k * i + j + np2] = 0;
  ker1[1][0][k * i + j + np2] = 0;
  ker2[1][0][k * i + j + np2] = 0;
  ker1[2][0][k * i + j + np2] = 0;
  ker1[3][0][k * i + j + np2] = 0;
  ker1[4][0][k * i + j + np2] = 0;
  ker1[5][0][k * i + j + np2] = 0;
  ker2[4][0][k * i + j + np2] = 0;

  return;
}

static void MaxwellSingle_insert_fmat(double **A1, double **A2, int i, int j,
                                      int nl, int a_bs, int a_bs2,
                                      double *dest) {
  int k;
  int m;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A1[0][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 0 * a_bs2, a_bs * sizeof(double));
    memcpy(&A1[1][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 1 * a_bs2, a_bs * sizeof(double));
    memcpy(&A1[2][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 3 * a_bs2, a_bs * sizeof(double));
    memcpy(&A1[3][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 4 * a_bs2, a_bs * sizeof(double));
    memcpy(&A1[4][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 5 * a_bs2, a_bs * sizeof(double));
    memcpy(&A1[5][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 7 * a_bs2, a_bs * sizeof(double));
    for (m = 0; m < a_bs; ++m) {
      A2[1][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[2 * a_bs2 + k * a_bs + m];
      A2[4][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[6 * a_bs2 + k * a_bs + m];
    }
  }
  return;
}

static void MaxwellSingle_insert_sym_fmat_diag(double **A, int i, int nl,
                                               int a_bs, int a_bs2,
                                               double *dest) {
  int k;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 0 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[1][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 1 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[2][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 3 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[3][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 4 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[4][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 5 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[5][(i * a_bs + k) * nl * a_bs + i * a_bs],
           dest + a_bs * k + 7 * a_bs2, a_bs * sizeof(double));
  }
  return;
}

static void MaxwellSingle_insert_sym_fmat_offdiag(double **A, int i, int j,
                                                  int nl, int a_bs, int a_bs2,
                                                  double *dest) {
  int k;
  int m;

  for (k = 0; k < a_bs; ++k) {
    memcpy(&A[0][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 0 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[1][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 1 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[2][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 3 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[3][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 4 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[4][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 5 * a_bs2, a_bs * sizeof(double));
    memcpy(&A[5][(i * a_bs + k) * nl * a_bs + j * a_bs],
           dest + a_bs * k + 7 * a_bs2, a_bs * sizeof(double));
    for (m = 0; m < a_bs; ++m) {
      A[0][(j * a_bs + m) * nl * a_bs + i * a_bs + k] = dest[k * a_bs + m];
      A[2][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[3 * a_bs2 + k * a_bs + m];
      A[3][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[4 * a_bs2 + k * a_bs + m];
      A[5][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[7 * a_bs2 + k * a_bs + m];

      A[1][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[2 * a_bs2 + k * a_bs + m];
      A[4][(j * a_bs + m) * nl * a_bs + i * a_bs + k] =
          dest[6 * a_bs2 + k * a_bs + m];
    }
  }
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
 */
void potMaxwell(double *rho, double **Pot, vector3 *R, int nr,
                discretization *disc, pdeproblem *pde) {
  meshdata *mesh = disc->mesh;
  int i;                      /* increment variable */
  int j;                      /* increment variable */
  int k;                      /* rotation and increment variable */
  int g = disc->g_pot;        /* quadrature degree */
  int p = mesh->geom->size(); /* number of parameter domains */
  int n = 1 << mesh->M;       /* n*n element per parameter domain */
  int nf = p * n * n;         /* total number of element */
  int a_bs = ((discretization *)disc)->a_bs;
  int N = a_bs * ((discretization *)disc)->mesh->nf;
  int zi; /* row index of the element */
  double wr;
  double wc;
  double c[2];
  vector3 e[2];
  double w[4];
  double w_div[2];
  double w_divh[2];
  double h = 1. / n; /* step size */
  double *kappa = ((discretization *)disc)->pde->kappa;
  double kappa2[2] = {kappa[0] * kappa[0] - kappa[1] * kappa[1],
                      2 * kappa[0] * kappa[1]};
  double kappa2abs2 = kappa2[0] * kappa2[0] + kappa2[1] * kappa2[1];
  double *d;
  double *d_dx;
  double *d_dy;
  double *rhol;
  double *mypot;
  vector2 s; /* left lower corner on square for zi */
  vector2 t; /* quadrature point on square */
  vector3 y; /* point on surface */
  std::array<vector3, 2> Chi_df;
  cubature *Q;                       /* quadrature formulas */
  const geometry &Chi = *mesh->geom; /* parametrizations */
  et_node *pE = mesh->E.patch[0];    /* pointer to the first element
                                      * on the lowest level of the
                                      * element tree */

  /*
   * initialize geometry and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1); /* Kubatur-Formeln */

  /*
   * allocate memory
   */
  (*Pot) = (double *)calloc(6 * nr, sizeof(double));
  rhol = (double *)calloc(6 * a_bs * nf, sizeof(double));

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
#pragma omp parallel default(none) private( \
    j, d, d_dx, d_dy, zi, s, k, t, y, mypot, Chi_df, w, w_div, w_divh, wr, wc, c, e)\
    shared(nr, a_bs, nf, pE, h, Q, g, Chi, disc, rhol, N, kappa2, kappa2abs2, R, kappa, Pot)
  {
    mypot = (double *)calloc(6 * nr, sizeof(double));
    d = (double *)calloc(a_bs, sizeof(double));
    d_dx = (double *)calloc(a_bs, sizeof(double));
    d_dy = (double *)calloc(a_bs, sizeof(double));
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
        /*
         * compute geometry informations
         */
        t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
        y = Chi[pE[zi].patch].f(t);

        Chi_df = Chi[pE[zi].patch].jacobian(t);

        /*
         * evaluate Galerkin solution on unit square, scale with quadrature
         * weight
         */
        memset(d, 0, a_bs * sizeof(double));
        disc->phiphi(d, Q[g].xi[k]);
        w[0] = Q[g].w[k] * myddot(a_bs, rhol + a_bs * zi, d);
        w[1] = Q[g].w[k] * myddot(a_bs, rhol + a_bs * zi + N, d);
        w[2] = Q[g].w[k] * myddot(a_bs, rhol + a_bs * zi + 2 * N, d);
        w[3] = Q[g].w[k] * myddot(a_bs, rhol + a_bs * zi + 3 * N, d);

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
        w_divh[0] = w_div[0];
        w_divh[1] = w_div[1];
        w_div[0] = (w_divh[0] * kappa2[0] + w_divh[1] * kappa2[1]) / kappa2abs2;
        w_div[1] = (w_divh[1] * kappa2[0] - w_divh[0] * kappa2[1]) / kappa2abs2;
        w_div[0] *= Q[g].w[k] / h;
        w_div[1] *= Q[g].w[k] / h;

        for (j = 0; j < nr; ++j) {
          /*
           * first part: with scalar fundamental solution
           */
          kernel(kappa, c, R[j], y);
          wr = w[0] * c[0] - w[1] * c[1];
          wc = w[0] * c[1] + w[1] * c[0];
          mypot[j] += Chi_df[0].x * wr;
          mypot[j + nr] += Chi_df[0].x * wc;
          mypot[j + 2 * nr] += Chi_df[0].y * wr;
          mypot[j + 3 * nr] += Chi_df[0].y * wc;
          mypot[j + 4 * nr] += Chi_df[0].z * wr;
          mypot[j + 5 * nr] += Chi_df[0].z * wc;
          wr = w[2] * c[0] - w[3] * c[1];
          wc = w[2] * c[1] + w[3] * c[0];
          mypot[j] += Chi_df[1].x * wr;
          mypot[j + nr] += Chi_df[1].x * wc;
          mypot[j + 2 * nr] += Chi_df[1].y * wr;
          mypot[j + 3 * nr] += Chi_df[1].y * wc;
          mypot[j + 4 * nr] += Chi_df[1].z * wr;
          mypot[j + 5 * nr] += Chi_df[1].z * wc;

          /*
           * second part: with vector valued gradient of fundamental solution
           */
          grad_kernel(kappa, e, R[j], y);
          mypot[j] += w_div[0] * e[0].x - w_div[1] * e[1].x;
          mypot[j + nr] += w_div[1] * e[0].x + w_div[0] * e[1].x;
          mypot[j + 2 * nr] += w_div[0] * e[0].y - w_div[1] * e[1].y;
          mypot[j + 3 * nr] += w_div[1] * e[0].y + w_div[0] * e[1].y;
          mypot[j + 4 * nr] += w_div[0] * e[0].z - w_div[1] * e[1].z;
          mypot[j + 5 * nr] += w_div[1] * e[0].z + w_div[0] * e[1].z;
        }
      }
    }
    free(d);
    free(d_dx);
    free(d_dy);
#pragma omp critical
    {
      for (zi = 0; zi < 6 * nr; ++zi) (*Pot)[zi] += mypot[zi];
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
 *  \author        Juergen Doelz
 */
static void MaxwellSingle_pot_eval(double *Pot, vector3 *R, int nr, double Qgw,
                                   vector3 y, vector3 n_y, double *rhol,
                                   double h, double *d, void *disc) {
  (void)Pot;
  (void)R;
  (void)nr;
  (void)Qgw;
  (void)y;
  (void)n_y;
  (void)rhol;
  (void)h;
  (void)d;
  (void)disc;
  assert(!"Cannot use this function. Use potMaxwell instead.");
  return;
}

static void MaxwellSingle_pot_normalize(double *Pot, int nr, double h) {
  int j;

  for (j = 0; j < 6 * nr; j++) Pot[j] *= h;
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
static void (**MaxwellSingle_Tmom_left_phi(void *voiddisc, int n))(double *, double,
                                                            double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  phi = (void (**)(double *, double, double))calloc(
      4, sizeof(void (*)(double *, double, double)));
  switch (n) {
    case 0:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 1:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 2:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
      break;
    case 3:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 4:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 5:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
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
static void (**MaxwellSingle_Tmom_right_phi(void *voiddisc, int n))(double *, double,
                                                             double) {
  discretization *disc = (discretization *)voiddisc;
  void (**phi)(double *, double, double);

  phi = (void (**)(double *, double, double))calloc(
      4, sizeof(void (*)(double *, double, double)));
  switch (n) {
    case 0:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 1:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
      break;
    case 2:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
      break;
    case 3:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi_dx;
      phi[3] = disc->phi;
      break;
    case 4:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
      break;
    case 5:
      phi[0] = disc->phi;
      phi[1] = disc->phi;
      phi[2] = disc->phi;
      phi[3] = disc->phi_dx;
      break;
    default:
      assert(!"this should not happen");
  }

  return phi;
}

static void MaxwellSingle_postproc(void *hmatfac, void *H) {
  ct_root **Hmat = (ct_root **)H;
  int i;
  int k;
  int aktn = 0;
  int a_bs = Hmat[0]->disc->a_bs;
  int minlvl = Hmat[0]->hmatset->min_bsize;
  int M = Hmat[0]->disc->mesh->M;
  int n = 1 << M;
  int np_max = ((hmatrixfactory *)hmatfac)->hmatset->np_max;
  int np2 = np_max * np_max;

  if (!Hmat[0][0].Tmom_left) return;

  for (k = 0; k < 6; ++k) {
    /*
     * apply Moment Tmom to vector x and store everything in transx[0]
     */
    aktn = a_bs * (1 << (2 * (minlvl)));
    for (i = np2; i < 2 * np2; ++i) {
      mydscal(aktn, n, Hmat[0][k].Tmom_left[i]);
      mydscal(aktn, n, Hmat[0][k].Tmom_right[i]);
    }
  }

  return;
}

Eigen::VectorXcd MaxwellSingle::MatVecImplementation(
    void *Hv, const Eigen::VectorXcd &x) const {
  // std::cout << "DOing mat Vec" << std::endl;
  ct_root *H = (ct_root *)Hv;
  int na = H->disc->na / 2;

  Eigen::VectorXd xr = Eigen::VectorXd(2 * na);
  Eigen::VectorXd yr = Eigen::VectorXd::Zero(2 * na);
  Eigen::VectorXcd y = Eigen::VectorXcd(na);
  xr << x.real(), x.imag();
  H2l2_HtimesVsmallMaxwell(H, (double *)xr.data(), (double *)yr.data());
  y.real() = yr.head(na);
  y.imag() = yr.tail(na);
  return y;
}

MaxwellSingle::MaxwellSingle(std::complex<double> in) {
  _pde.kappa[0] = in.real();
  _pde.kappa[1] = in.imag();
  _pde.pdetype = Maxwell;
  _pde.name = "MaxwellSingle";
  _pde.nct = 6;
  _pde.np_max_fac = 2;
  _pde.Tmom_left_phi = &MaxwellSingle_Tmom_left_phi;
  _pde.Tmom_right_phi = &MaxwellSingle_Tmom_right_phi;
  _pde.symflags = "SNSSNS";
  _pde.quadrature_accuracy = &MaxwellSingle_quadrature_accuracy;
  _pde.quadrature_bufsize = 8;
  _pde.g_far = &MaxwellSingle_g_far;
  _pde.g_pot = &MaxwellSingle_g_pot;
  _pde.init_randwerte = &MaxwellSingle_init_randwerte;
  _pde.free_randwerte = &MaxwellSingle_free_randwerte;
  _pde.IntPhi0 = &MaxwellSingle_IntPhi0;
  _pde.IntPhi1 = &MaxwellSingle_IntPhi1;
  _pde.IntPhi2 = &MaxwellSingle_IntPhi2;
  _pde.IntPhi3 = &MaxwellSingle_IntPhi3;
  _pde.IntPhi4 = &MaxwellSingle_IntPhi4;
  _pde.quadrature_order = &MaxwellSingle_quadrature_order;
  _pde.interpol_kernel = &MaxwellSingle_interpol_kernel;
  _pde.insert_fmat = &MaxwellSingle_insert_fmat;
  _pde.insert_sym_fmat_diag = &MaxwellSingle_insert_sym_fmat_diag;
  _pde.insert_sym_fmat_offdiag = &MaxwellSingle_insert_sym_fmat_offdiag;
  _pde.pot_eval = &MaxwellSingle_pot_eval;
  _pde.pot_normalize = &MaxwellSingle_pot_normalize;
  _pde.postproc = &MaxwellSingle_postproc;
}
}  // namespace Bembel
