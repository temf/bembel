// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "constants.h"

namespace Bembel {
/*******************************************************************************
 *  Quadrature properties *
 *******************************************************************************/

/**
 *  \brief         Tolerance for point comparison.
 *
 *  The tolerance is used to determine if two patches lie next to each other.
 */
const double tol_points = 1e-6;

/**
 *  \brief         Tolerance for solvers.
 */
const double tol_solver = 1e-8;

/*******************************************************************************
 *  Block cluster tree properties *
 *******************************************************************************/

/**
 *  \brief         If M-l<= min_bsize, then the whole block shall be a full
 *                 matrix.
 *
 *  Can be modified by the program at runtime.
 */
int min_bsize = 1;
}  // namespace Bembel