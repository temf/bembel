// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_DISCRETIZATION__
#define __BEMBEL_C_DISCRETIZATION__

#include "PDEproblem.hpp"
#include "meshdata.h"
#include "phi.h"
#include "spline/patch.h"

namespace Bembel {
struct _proj_info {
  std::vector<double> vals;
  std::vector<int> rows;
  std::vector<int> cols;
  int num_c_dof;
  int num_dc_dof;
};  // namespace Bembelstruct_proj_info

struct discretization {
  pdeproblem *pde;
  meshdata *mesh;
  sparse *T;
  int na;
  int real_na;
  int a_o;
  int a_bs;
  int a_bs2;

  int quadrature_accuracy;
  int g_far;
  int g_pot;

  void (*phi)(double *, double, double);
  void (*phi_dx)(double *, double, double);

  void (*phiphi)(double *, vector2);
  void (*phiphi_dx)(double *, vector2);
  void (*phiphi_dy)(double *, vector2);

  void (*Phi_times_Phi)(double *, double, vector2, vector2);
  void (*VPhi_scal_VPhi)(double *, double, vector2, vector2, vector3, vector3,
                         vector3, vector3);
  void (*Div_Phi_times_Div_Phi)(double *, double, vector2, vector2);
  void (*Curl_Phi_times_Curl_Phi)(double *, double, vector2, vector2, vector3,
                                  vector3, vector3, vector3);
};

void proj_distr_et(discretization *, double *, double *);
void proj_restr_et(discretization *, double *, double *);

void free_discretization(discretization *);

discretization get_discretization_0P(pdeproblem *pde, meshdata *mesh);
discretization get_discretization_1P(pdeproblem *pde, meshdata *mesh);
discretization get_discretization_2B(pdeproblem *pde, meshdata *mesh);

discretization get_discretization_NB(int deg, int kntrep, pdeproblem *pde,
                                     meshdata *mesh);
discretization get_discretization_NB_Laplace_cont(const int deg,
                                                  const int kntrep,
                                                  pdeproblem *pde,
                                                  meshdata *mesh);

int init_projector_Laplace(discretization *disc, et_node *E, const int M,
                           const int pp1, const int kntrep);
int init_projector_Laplace_cont(discretization *disc, et_node *E, const int M,
                                const int pp1, const int kntrep);
int init_projector_Helmholtz(discretization *disc, et_node *E, const int M,
                             const int pp1, const int kntrep);
int init_projector_Maxwell(discretization *disc, et_node *E, const int M,
                           const int pp1, const int kntrep);
_proj_info init_projector_base(discretization *disc, const int pp1x,
                               const int pp1y, const int kntrep);

inline int init_projector(discretization *disc, et_node *E, const int M,
                          const int pp1, const int kntrep) {
  switch (disc->pde->pdetype) {
    case (Laplace):
      return init_projector_Laplace(disc, E, M, pp1, kntrep);
    case (Helmholtz):
      return init_projector_Helmholtz(disc, E, M, pp1, kntrep);
    case (Maxwell):
      return init_projector_Maxwell(disc, E, M, pp1, kntrep);
    default:
      return 0;
  }
}

inline void free_projector(sparse *T) {
  free_sparse(T);
  free(T);

  return;
}
}  // namespace Bembel
#endif
