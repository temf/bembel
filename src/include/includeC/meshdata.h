// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_MESHDATA__
#define __BEMBEL_C_MESHDATA__

#include "Geometry.hpp"
#include "element_tree.h"
namespace Bembel {
typedef struct {
  int M;
  int p;
  int n;
  int np;
  int nf;
  int na;
  int **F;
  double *R;
  vector3 *D;
  vector3 *P;

  const geometry *geom;

  et_root E;

} meshdata;

meshdata get_meshdata(const geometry &geom, int M);

void free_meshdata(meshdata *mesh);
}  // namespace Bembel
#endif
