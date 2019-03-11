// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include "meshdata.h"

namespace Bembel {
meshdata get_meshdata(const geometry &geom, int M) {
  meshdata mesh;

  mesh.geom = &geom;
  mesh.M = M;
  mesh.p = geom.size();
  mesh.n = 1 << mesh.M;
  mesh.nf = mesh.p * mesh.n * mesh.n;
  mesh.np = gennet(&mesh.P, &mesh.F, mesh.M, geom);
  mesh.na = init_element_tree(&mesh.E, mesh.P, mesh.F, mesh.p, mesh.M, mesh.np);
  dist(mesh.P, mesh.F, &mesh.D, &mesh.R, mesh.M, mesh.p);

#ifdef _BEMBEL_PRINT_INFO
  printf("                           %dth level\n", mesh.M);
  printf("                           %d elements\n", mesh.nf);
  printf("                           %d points\n", mesh.np);
#endif
  return mesh;
}

void free_meshdata(meshdata *mesh) {
  free(mesh->P);
  free(mesh->D);
  free(mesh->R);
  free_element_tree(&mesh->E);
  free_patchlist(&mesh->F, mesh->nf);

  return;
}
}  // namespace Bembel