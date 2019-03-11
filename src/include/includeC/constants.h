// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_CONSTANTS_
#define __BEMBEL_C_CONSTANTS_

namespace Bembel {

#if !defined(M_PI)
#define M_PI 3.1415926535897932385
#endif

extern const double tol_points;
extern const double tol_solver;

/*
 * Hmatrix arithmetic properties
 */
/*
 * extern const double arithmetics_prec; extern int arithmetics_kmax; extern
 * int np_max;
 */
}  // namespace Bembel
#endif
