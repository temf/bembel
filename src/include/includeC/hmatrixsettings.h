// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_hmatrixsettings__
#define __BEMBEL_C_hmatrixsettings__

#include "discretization.h"
#include "pdeproblem.h"
namespace Bembel {
typedef struct {
  /**
   *  \brief         Cut-off tolerance for the singular values in recompression
   *                 of rk-matrices.
   */
  double arithmetics_prec;

  /**
   *  \brief         Maximal rank of the Hmatrix after Hmatrix arithmetics.
   *                 Usually, this variable will be overwritten by the choice of
   *                 the user at runtime.
   */
  int arithmetics_kmax;

  /**
   *  \brief         Fix polynomial degree for multipole interpolation of the
   *                 kernel.
   */
  int np_max;

  int min_bsize;

  double eta;

} hmatrixsettings;

hmatrixsettings get_hmatrixsettings(int np_max, discretization *disc);
void hmatrixsettings_print(hmatrixsettings *hmatset);
}  // namespace Bembel
#endif
