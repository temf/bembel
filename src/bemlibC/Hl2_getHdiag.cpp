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
 * \brief Evaluates the diagonal of H-matrix and stores the result in the
 *        vector d.
 *
 *  \param[in]     H        Hmatrix from which the diagonal shall be taken
 *  \param[out]    d        vector in which the diagonal will be stored
 *  \param[in]     p        number of children on the highest level of current
 *                          block cluster tree
 *  \param[in]     nf       block size of matrix H
 *
 */
int Hl2_getHdiag(ct_node *H, double *d, int p, int N) {
  int i = 0;
  int n = N / p;
  double *pd = NULL;

  if (H->pch) {
    for (i = 0; i < p; ++i) {
      pd = d + i * n;
      Hl2_getHdiag(H->pch[i * p + i], pd, 4, n);
    }
  } else if (H->fmat) {
    assert(N == H->fmat->bs);
    for (i = 0; i < N; ++i) d[i] = H->fmat->A[i * H->fmat->bs + i];
  } else {
    assert(!"Found an rk-matrix on the diagonal. Is this possible?");
  }

  return 0;
}

/**
 * \brief Evaluates the ElementDiagonal of H-matrix and stores the result in the
 *        vector d.
 *
 *  \param[in]     H        Hmatrix from which the diagonal shall be taken
 *  \param[out]    d        vector in which the diagonal will be stored
 *  \param[in]     p        number of children on the highest level of current
 *                          block cluster tree
 *  \param[in]     nf       block size of matrix H
 *  \param[in]     a_bs     a_bs...
 */
int Hl2_getElementDiag(ct_node *H, double *d, int p, int N, int a_bs) {
  int i = 0;
  int j = 0;
  int a_bs2 = a_bs * a_bs;
  int n = N / p;
  double *pd = NULL;

  if (H->pch) {
    for (i = 0; i < p; ++i) {
      pd = d + i * n * a_bs;
      Hl2_getElementDiag(H->pch[i * p + i], pd, 4, n, a_bs);
    }
  } else if (H->fmat) {
    assert(N == H->fmat->bs);
    for (i = 0; i < N / a_bs; ++i)
      for (j = 0; j < a_bs; ++j)
        memcpy(d + i * a_bs2 + j * a_bs,
               H->fmat->A + (i * N + i) * a_bs + j * N, a_bs * sizeof(double));
  } else {
    assert(!"Found an rk-matrix on the diagonal. Is this possible?");
  }

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x for element diagonal matrices.
 *
 *  \param[in]     G              Vector with entries of matrix.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *  \param[in]     nf             Number of elements.
 */
int EDtimesV(double *G, double *x, double *y, int nf, int a_bs) {
  int i;
  int a_bs2 = a_bs * a_bs;

  /* apply big matrix to big vectors */
  for (i = 0; i < nf; ++i)
    myqdgemv('N', a_bs, G + i * a_bs2, x + i * a_bs, y + i * a_bs);

  return 0;
}
}  // namespace Bembel