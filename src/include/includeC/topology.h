// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_TOPOLOGY__
#define __BEMBEL_C_TOPOLOGY__

#include "Geometry.hpp"
#include "myvector.h"
namespace Bembel {
void unify(vector3 *d, double *r, vector3 d1, double r1, vector3 d2, double r2);
/*
 * bildet die Vereinigung K(d,r) = K(d1,r1) \cup K(d2,r2)
 */

int gennet(vector3 **P, int ***F, int m, const geometry &geom);
/*
 * berechnet Punkt- und Patchliste
 */

void dist(vector3 *P, int **F, vector3 **D, double **R, int m, const int p);
/*
 * berechnet Einschliessungen der Patches
 */

void free_patchlist(int ***F, int nf);
/*
 * gibt den Speicherplatz der Patchliste frei
 */

int search_point(vector3 x, vector3 *P, int *pz);

void init_grid(vector3 **P, int ***F, int m, int *np, int *nf,
               const geometry &Chi);

void refine_grid(vector3 **P, int ***F, int m, int *np, int *nf,
                 const geometry &Chi);
}  // namespace Bembel
#endif
