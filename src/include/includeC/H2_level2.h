// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_H2LEVEL2_
#define __BEMBEL_C_H2LEVEL2_

#include <stdio.h>
#include <string.h>

#include "H2transforms.h"
#include "cluster_tree.h"
#include "constants.h"
#include "mycblas.h"
#include "sparse.h"

namespace Bembel {
int Hl2_FtimesV(char trans, ct_fmat *pf, double *px, double *py);
int H2l2_RtimesV(char trans, ct_rkmat *prk, double *px, double *py);
int H2l2_HtimesV(char trans, ct_node *H, et_node *nc, et_node *nr,
                 H2tree *transx, H2tree *transy, double *x, double *y, int p,
                 int nf, int M, int min_bsize);
int H2l2_HtimesVbig(char trans, ct_root *H, double *x, double *y);
int H2l2_HtimesVsmall(char trans, ct_root *H, double *x, double *y);
int H2l2_HtimesVsmallReal(ct_root *H, double *x, double *y);
int H2l2_HtimesVsmallComplex(ct_root *H, double *x, double *y);
int H2l2_HtimesVsmallMaxwell(ct_root *H, double *x, double *y);
int H2l2_HtimesVbigMaxwell(ct_root *H, double *x, double *y);
int Hl2_getHdiag(ct_node *H, double *d, int p, int nf);
int Hl2_getHinvdiag(ct_node *H, double *d, int p, int nf);
int Hl2_updateHdiag(ct_node *H, double *d, int p, int nf);
}  // namespace Bembel
#endif
