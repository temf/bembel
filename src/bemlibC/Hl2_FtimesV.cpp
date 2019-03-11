// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "H_level2.h"

namespace Bembel {
/**
 * \brief y+=F*x or y+=F^T*x
 *
 *  \param[in]     trans    If trans=='T', then pf will be transposed. For all
 *                          other values this variable has no influence.
 *  \param[in]     pf       full matrix which will be applied to x
 *  \param[in]     px       vector to which pf is applied to
 *  \param[in,out] py       vector to which the product of pf and x is added
 */
int Hl2_FtimesV(char trans, ct_fmat *pf, double *px, double *py) {
  myqdgemv(trans, pf->bs, pf->A, px, py);

  return 0;
}
}  // namespace Bembel