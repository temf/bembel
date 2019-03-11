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
 * \brief y+=H*x or y+=H^T*x
 *
 *  \param[in]     trans    Character of matrix, possible are 'N', 'T',
 *                          'L', 'U' (which is L^T) and 'S' (symmetrix).
 *  \param[in]     H        Hmatrix which will be applied to x
 *  \param[in]     x        vector to which H is applied to
 *  \param[in,out] y        vector to which the product of H and x is added
 *  \param[in]     p        number of children on the highest level of current
 *                          block cluster tree
 *  \param[in]     nf       block size of matrix H
 *
 *  If compiled with OpenMP, then the multiplication will be performed parallel.
 *  The number of threads is bounded by p or by the maximum number of allowed
 *  threads.
 */
int Hl2_HtimesV(char trans, ct_node *H, double *x, double *y, int p, int nf) {
  int i = 0;
  int j = 0;
  int n = nf / p;
  double *px = NULL;
  double *py = NULL;

  if (H->pch) {
    switch (trans) {
      case 'N':
#pragma omp parallel for default(none) private(i, j, px, py) \
    shared(p, n, H, x, y, trans)
        for (i = 0; i < p; ++i) {
          py = y + i * n;
          for (j = 0; j < p; ++j) {
            px = x + j * n;
            Hl2_HtimesV('N', H->pch[i * p + j], px, py, 4, n);
          }
        }
        break;
      case 'T':
#pragma omp parallel for default(none) private(i, j, px, py) \
    shared(p, n, H, x, y, trans)
        for (i = 0; i < p; ++i) {
          py = y + i * n;
          for (j = 0; j < p; ++j) {
            px = x + j * n;
            Hl2_HtimesV('T', H->pch[j * p + i], px, py, 4, n);
          }
        }
        break;
      case 'L':
#pragma omp parallel for default(none) private(i, j, px, py) \
    shared(p, n, H, x, y, trans)
        for (i = 0; i < p; ++i) {
          py = y + i * n;
          for (j = 0; j < i; ++j) {
            px = x + j * n;
            Hl2_HtimesV('N', H->pch[i * p + j], px, py, 4, n);
          }
          px = x + i * n;
          Hl2_HtimesV('L', H->pch[i * p + i], px, py, 4, n);
        }
        break;
      case 'U':
#pragma omp parallel for default(none) private(i, j, px, py) \
    shared(p, n, H, x, y, trans)
        for (i = 0; i < p; ++i) {
          py = y + i * n;
          px = x + i * n;
          Hl2_HtimesV('U', H->pch[i * p + i], px, py, 4, n);
          for (j = i + 1; j < p; ++j) {
            px = x + j * n;
            Hl2_HtimesV('T', H->pch[j * p + i], px, py, 4, n);
          }
        }
        break;
      case 'S':
#pragma omp parallel for default(none) private(i, j, px, py) \
    shared(p, n, H, x, y)
        for (i = 0; i < p; ++i) {
          py = y + i * n;
          for (j = 0; j < i; ++j) {
            px = x + j * n;
            Hl2_HtimesV('N', H->pch[i * p + j], px, py, 4, n);
          }
          px = x + i * n;
          Hl2_HtimesV('S', H->pch[i * p + i], px, py, 4, n);
          for (j = i + 1; j < p; ++j) {
            px = x + j * n;
            Hl2_HtimesV('T', H->pch[j * p + i], px, py, 4, n);
          }
        }
        break;
      default:
        printf("Hl2_HtimesV says: This shouldn't happen!!!\n");
        exit(1);
        break;
    }
  } else if (H->rkmat) {
    Hl2_RtimesV((trans == 'T' || trans == 'U') ? 'T' : 'N', H->rkmat, x, y);
  } else {
    Hl2_FtimesV((trans == 'T' || trans == 'U') ? 'T' : 'N', H->fmat, x, y);
  }

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *  \param[in]     p              Number of patches of H.
 *  \param[in]     M              Number of levels.
 *  \param[in]     E              Element list.
 */
int Hl2_HtimesVsmall(char trans, ct_node *H, double *x, double *y, int p, int M,
                     et_node *E) {
  int n = 1 << M;     /* n*n elements per patch */
  int nf = p * n * n; /* block size of H */
  int N = a_bs * nf;  /* size of big vectors */
  double *lx;         /* source big vector */
  double *ly;         /* target big vector */

  /*
   * get temporary memory
   */
  lx = calloc(N, sizeof(double));
  ly = calloc(N, sizeof(double));

  /*
   * distribute vectors
   */
  proj_distr_et(x, lx, E, M);

  /*
   * apply big matrix to big vectors
   */
  Hl2_HtimesV(trans, H, lx, ly, p, N);

  /*
   * restrict destination vector
   */
  proj_restr_et(ly, y, E, M);

  free(lx);
  free(ly);

  return 0;
}
}  // namespace Bembel