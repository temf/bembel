// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_H_level2__
#define __BEMBEL_C_H_level2__

#include <stdio.h>
#include <string.h>

#include "cluster_tree.h"
#include "constants.h"
#include "mycblas.h"
#include "sparse.h"
namespace Bembel {
int Hl2_FtimesV(char trans, ct_fmat *pf, double *px, double *py);
int Hl2_RtimesV(char trans, ct_rkmat *prk, double *px, double *py);
int Hl2_HtimesV(char trans, ct_node *H, double *x, double *y, int p, int nf);
int Hl2_HtimesVsmall(char trans, ct_node *H, double *x, double *y, int p, int M,
                     et_node *E);
int Hl2_extHtimesV(char trans, ct_node **pH1, ct_node **pH2, int nblocks,
                   double *x, double *y, double *swap, int p, int nf);
int Hl2_getHdiag(ct_node *H, double *d, int p, int nf);
int Hl2_getElementDiag(ct_node *H, double *d, int p, int N, int a_bs);
int EDtimesV(double *G, double *x, double *y, int nf, int a_bs);
int Hl2_getHinvdiag(ct_node *H, double *d, int p, int nf);
int Hl2_updateHdiag(ct_node *H, double *d, int p, int nf);
}  // namespace Bembel
#endif