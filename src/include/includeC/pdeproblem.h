// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_PDEPROBLEM__
#define __BEMBEL_C_PDEPROBLEM__

#include <stdlib.h>
#include <string>
#include "constants.h"
#include "cubature.h"
#include "data.h"
#include "element_tree.h"
#include "gamma.h"
#include "mycblas.h"
#include "myvector.h"
#include "sparse.h"
#include "trafos.h"

namespace Bembel {
typedef struct {
  vector3 Chi;
  vector3 n_Chi;
  vector3 dx_Chi;
  vector3 dy_Chi;
  double det_dChi;
} randwerte;

enum Pdetype { Laplace, Helmholtz, Maxwell };

struct pdeproblem {
  Pdetype pdetype;
  std::string name;
  std::string symflags;
  int nct;

  int np_max_fac;
  void (**(*Tmom_left_phi)(void *, int))(double *, double, double);
  void (**(*Tmom_right_phi)(void *, int))(double *, double, double);

  double kappa[2];

  int (*quadrature_accuracy)(int);
  int quadrature_bufsize;

  int (*g_far)(int);
  int (*g_pot)(int);

  /*
   * quadrature routines
   */
  void (*init_randwerte)(randwerte ***, cubature *, const geometry &, int);
  void (*free_randwerte)(randwerte ***, int, const int);
  void (*IntPhi0)(double *, int, int, cubature *, randwerte **, void *);
  void (*IntPhi1)(double *, vector2, vector2, double, const Spl::Patch &,
                  const Spl::Patch &, cubature *, void *);
  void (*IntPhi2)(double *, vector2, double, const Spl::Patch &, cubature *,
                  void *);
  void (*IntPhi3)(double *, vector2, vector2, double, int, int,
                  const Spl::Patch &, const Spl::Patch &, cubature *, void *);
  void (*IntPhi4)(double *, vector2, vector2, double, int, int,
                  const Spl::Patch &, const Spl::Patch &, cubature *, void *);
  int (*quadrature_order)(double, double, int, int);

  /*
   * tree leaf routines
   */
  void (*interpol_kernel)(double ***, double ***, const Spl::Patch &,
                          const Spl::Patch &, vector2, vector2, int, int, int,
                          double[2]);
  void (*insert_fmat)(double **, double **, int, int, int, int, int, double *);
  void (*insert_sym_fmat_diag)(double **, int, int, int, int, double *);
  void (*insert_sym_fmat_offdiag)(double **, int, int, int, int, int, double *);

  /*
   * potential evaluation
   */
  void (*pot_eval)(double *, vector3 *, int, double, vector3, vector3, double *,
                   double, double *, void *);
  void (*pot_normalize)(double *, int, double);

  /*
   * postprocessing
   */
  void (*postproc)(void *, void *);
  //
};

void pdeproblem_print(pdeproblem *pde);
}  // namespace Bembel
#endif
