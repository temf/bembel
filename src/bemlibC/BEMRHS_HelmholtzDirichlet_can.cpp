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
 *
 */
Eigen::VectorXcd computeRhs(
    Discretization<HelmholtzSingle> &ddisc,
    std::function<std::complex<double>(Eigen::Vector3d, std::complex<double>)>
        fun) {
  double *rhs;

  discretization *disc = &(ddisc.get_disc());

  meshdata *mesh = disc->mesh;
  int i;                            /* increment variable */
  int k;                            /* increment variable */
  int g = disc->g_far;              /* quadrature degree */
  const int p = mesh->geom->size(); /* dumber of parameter domains */
  int n = 1 << disc->mesh->M;       /* n*n patches per parameter
                                     * domain */
  int nf = p * n * n;               /* number of elements */
  int zi;                           /* row index of the element */
  double c[2];                      /* output of quadrature routines */
  double h = 1. / n;                /* step size */
  double w;
  double *kappa = disc->pde->kappa;
  double *d;
  double *rhsl;
  double *rhslc;
  vector2 s;   /* left lower corner on square for zi */
  vector2 t;   /* quadrature point on square */
  cubature *Q; /* quadrature formulas */
  const geometry &Chi = *mesh->geom;
  et_node *pE = mesh->E.patch[0]; /* pointer to the first element
                                   * on the */

  auto local_data = [fun](double out[2], vector3 &&in, double kapl[2]) {
    std::complex<double> tmp = fun(Eigen::Vector3d(in.x, in.y, in.z),
                                   std::complex<double>(kapl[0], kapl[1]));
    out[0] = tmp.real();
    out[1] = tmp.imag();
    return;
  };

  /*
   * initialize and quadrature formulas
   */
  init_Gauss_Square(&Q, g + 1);

  /*
   * allocate memory
   */
  rhsl = (double *)calloc(nf * 2 * disc->a_bs, sizeof(double));
  rhslc = rhsl + nf * disc->a_bs;

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
      local_data(c, Chi[pE[zi].patch].f(t), kappa);
      w = h * Q[g].w[k] * vector3_norm(Chi[pE[zi].patch].n_f(t));
      c[0] *= w;
      c[1] *= w;
      memset(d, 0, disc->a_bs * sizeof(double));
      disc->phiphi(d, Q[g].xi[k]);
      mydaxpy(disc->a_bs, c[0], d, rhsl + disc->a_bs * zi);
      mydaxpy(disc->a_bs, c[1], d, rhslc + disc->a_bs * zi);
    }
  }

  rhs = (double *)calloc(disc->na, sizeof(double));
  proj_restr_et(disc, rhsl, rhs);

  const int sz = disc->na / 2;
  Eigen::VectorXcd out(sz);

  for (int i = 0; i < sz; i++) {
    out(i) = std::complex<double>(rhs[i], rhs[i + sz]);
  }

  /*
   * free memory
   */
  free_Gauss_Square(&Q, g + 1);
  free(rhsl);
  free(d);
  free(rhs);
  return out;
}
}  // namespace Rhs
}  // namespace Bembel
