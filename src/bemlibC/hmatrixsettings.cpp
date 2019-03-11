// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "hmatrixsettings.h"

namespace Bembel {
hmatrixsettings get_hmatrixsettings(int np_max, discretization *disc) {
  hmatrixsettings hmatset;

  hmatset.arithmetics_prec = 1e-12;
  hmatset.np_max = np_max;
  hmatset.eta = 1.6;
  hmatset.min_bsize =
      (int)(log((disc->pde->np_max_fac * hmatset.np_max * hmatset.np_max - 1) /
                disc->a_bs) /
            log(4));
  hmatset.min_bsize = hmatset.min_bsize < 1 ? 1 : hmatset.min_bsize;
  hmatset.arithmetics_kmax =
      disc->pde->np_max_fac * hmatset.np_max * hmatset.np_max;

  return hmatset;
}

void hmatrixsettings_print(hmatrixsettings *hmatset) {
  printf("%26s arithmetics_prec is %g\n", "", hmatset->arithmetics_prec);
  printf("%26s np_max is %d\n", "", hmatset->np_max);
  printf("%26s eta is %g\n", "", hmatset->eta);
  printf("%26s min_bsize is %d\n", "", hmatset->min_bsize);
  printf("%26s arithmetics_kmax is %d\n", "", hmatset->arithmetics_kmax);
  return;
}
}  // namespace Bembel