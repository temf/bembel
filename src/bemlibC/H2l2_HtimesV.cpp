// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "H2_level2.h"

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
int H2l2_HtimesV(char trans, ct_node *H, et_node *nc, et_node *nr,
                 H2tree *transx, H2tree *transy, double *x, double *y, int p,
                 int nf, int M, int rank) {
  int i = 0;
  int j = 0;
  int n = nf / p;

  if (H->pch) switch (trans) {
      case 'N':
        for (i = 0; i < p; ++i)
          for (j = 0; j < p; ++j) {
            {
              H2l2_HtimesV('N', H->pch[i * p + j], nc->son[i], nr->son[j],
                           &transx->sons[j], &transy->sons[i], x, y, 4, n, M,
                           rank);
            }
          }
        break;
      case 'T':
        for (i = 0; i < p; ++i)
          for (j = 0; j < p; ++j) {
            {
              H2l2_HtimesV('T', H->pch[j * p + i], nc->son[i], nr->son[j],
                           &transx->sons[j], &transy->sons[i], x, y, 4, n, M,
                           rank);
            }
          }
        break;
      case 'S':
        for (i = 0; i < p; ++i) {
          for (j = 0; j < i; ++j) {
            H2l2_HtimesV('N', H->pch[i * p + j], nc->son[i], nr->son[j],
                         transx->sons + j, transy->sons + i, x, y, 4, n, M,
                         rank);
          }
          H2l2_HtimesV('S', H->pch[i * p + i], nc->son[i], nr->son[i],
                       transx->sons + i, transy->sons + i, x, y, 4, n, M, rank);
          for (j = i + 1; j < p; ++j) {
            H2l2_HtimesV('T', H->pch[j * p + i], nc->son[i], nr->son[j],
                         transx->sons + j, transy->sons + i, x, y, 4, n, M,
                         rank);
          }
        }
        break;
      default:
        assert(!("This shouldn't happen!!!"));
        break;
    }

  else if (H->rkmat) {
    H2l2_RtimesV((trans == 'T' || trans == 'U') ? 'T' : 'N', H->rkmat,
                 transx->data, transy->data);
  } else {
    Hl2_FtimesV((trans == 'T' || trans == 'U') ? 'T' : 'N', H->fmat,
                &(x[nf * nr->number]), &(y[nf * nc->number]));
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
 */
int H2l2_HtimesVbig(char trans, ct_root *H, double *x, double *y) {
  meshdata *mesh = H->disc->mesh;
  int i;
  int j;
  const int n = mesh->n; /* n*n elements per patch */
  const int p = mesh->geom->size();
  const int nf = p * n * n;         /* block size of H */
  const int N = H->disc->a_bs * nf; /* size of big vectors */
  const int M = mesh->M;
  const int rank = H->rank;
  int depth = M - H->hmatset->min_bsize;
  H2tree *transx; /* H2 splitted big source vector */
  H2tree *transy; /* H2 splitted big target vector */

  if (depth < 0) depth = 0;

  /*
   * Forward transform
   */
  /*memset(transx, 0, sizeof(H2tree));
  memset(transy, 0, sizeof(H2tree));*/
  transx = (H2tree *)calloc(1, sizeof(H2tree));
  transy = (H2tree *)calloc(1, sizeof(H2tree));
  new_H2tree(transx, depth + 1, p, rank);
  new_H2tree(transy, depth + 1, p, rank);
  ForwardTransform(H, x, transx);

  /*
   * apply big matrix to big vectors
   */
#if 0
  H2l2_HtimesV(trans, H->root, mesh->E.patch, mesh->E.patch, transx, transy, x,
               y, H->p, N, M, rank);
#else
  switch (trans) {
    case 'N':
#pragma omp parallel for default(shared) private(i, j)
      for (i = 0; i < p; ++i)
        for (j = 0; j < p; ++j) {
          {
            H2l2_HtimesV('N', H->root->pch[i * p + j], mesh->E.patch[i],
                         mesh->E.patch[j], &transx->sons[j], &transy->sons[i],
                         x, y, 4, N / p, M, rank);
          }
        }
      break;
    case 'T':
#pragma omp parallel for default(shared) private(i, j)
      for (i = 0; i < p; ++i)
        for (j = 0; j < p; ++j) {
          {
            H2l2_HtimesV('T', H->root->pch[j * p + i], mesh->E.patch[i],
                         mesh->E.patch[j], &transx->sons[j], &transy->sons[i],
                         x, y, 4, N / p, M, rank);
          }
        }
      break;
    case 'S':
#pragma omp parallel for default(shared) private(i, j)
      for (i = 0; i < p; ++i) {
        for (j = 0; j < i; ++j) {
          H2l2_HtimesV('N', H->root->pch[i * p + j], mesh->E.patch[i],
                       mesh->E.patch[j], transx->sons + j, transy->sons + i, x,
                       y, 4, N / p, M, rank);
        }
        H2l2_HtimesV('S', H->root->pch[i * p + i], mesh->E.patch[i],
                     mesh->E.patch[i], transx->sons + i, transy->sons + i, x, y,
                     4, N / p, M, rank);
        for (j = i + 1; j < p; ++j) {
          H2l2_HtimesV('T', H->root->pch[j * p + i], mesh->E.patch[i],
                       mesh->E.patch[j], transx->sons + j, transy->sons + i, x,
                       y, 4, N / p, M, rank);
        }
      }
      break;
    default:
      assert(!("This shouldn't happen!!!"));
      break;
  }
#endif

  /*
   * Backward transform
   */
  BackwardTransform(H, y, transy);
  free_H2tree(transx);
  free_H2tree(transy);

  free(transx);
  free(transy);

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *
 */
int H2l2_HtimesVbigComplex(ct_root *H, double *x, double *y) {
  int N = H->disc->mesh->nf * H->disc->a_bs;
  double *xc = x + N;
  double *yc = y + N;
  double *temp;
  ct_root *Hc = H + 1;

  /*
   * R -= C*C
   */
  temp = (double *)calloc(N, sizeof(double));
  H2l2_HtimesVbig(Hc->sym, Hc, xc, temp);
  mydaxpy(N, -1., temp, y);
  free(temp);
  /*
   * R+=R*R
   */
  H2l2_HtimesVbig(H->sym, H, x, y);
  /*
   * C+=C*R
   */
  H2l2_HtimesVbig(Hc->sym, Hc, x, yc);
  /*
   * C+=R*C
   */
  H2l2_HtimesVbig(H->sym, H, xc, yc);

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 */
int H2l2_HtimesVsmall(char trans, ct_root *H, double *x, double *y) {
  meshdata *mesh = H->disc->mesh;
  const int n = mesh->n; /* n*n elements per patch */
  const int p = mesh->geom->size();
  const int nf = p * n * n;         /* block size of H */
  const int N = H->disc->a_bs * nf; /* size of big vectors */
  double *lx;                       /* source big vector */
  double *ly;                       /* target big vector */

  /*
   * transform short vectors to long vectors
   */
  lx = (double *)calloc(N, sizeof(double));
  ly = (double *)calloc(N, sizeof(double));
  proj_distr_et(H->disc, x, lx);

  H2l2_HtimesVbig(trans, H, lx, ly);

  /*
   * transform long vector to short vector
   */
  proj_restr_et(H->disc, ly, y);

  free(lx);
  free(ly);

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 */
int H2l2_HtimesVsmallReal(ct_root *H, double *x, double *y) {
  H2l2_HtimesVsmall(H->sym, H, x, y);
  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *
 */
int H2l2_HtimesVsmallComplex(ct_root *H, double *x, double *y) {
  meshdata *mesh = H->disc->mesh;
  const int n = mesh->n; /* n*n elements per patch */
  const int p = mesh->geom->size();
  const int nf = p * n * n;         /* block size of H */
  const int N = H->disc->a_bs * nf; /* size of big vectors */
  double *lx;
  double *ly;
  double *lxc;
  double *lyc;
  double *temp;
  ct_root *Hc = H + 1;

  /*
   * transform short vectors to long vectors
   */
  lx = (double *)calloc(2 * N, sizeof(double));
  ly = (double *)calloc(2 * N, sizeof(double));
  proj_distr_et(H->disc, x, lx);
  lxc = lx + N;
  lyc = ly + N;

  /*
   * R -= C*C
   */
  temp = (double *)calloc(N, sizeof(double));
  H2l2_HtimesVbig(Hc->sym, Hc, lxc, temp);
  mydaxpy(N, -1., temp, ly);
  free(temp);
  /*
   * R+=R*R
   */
  H2l2_HtimesVbig(H->sym, H, lx, ly);
  /*
   * C+=C*R
   */
  H2l2_HtimesVbig(Hc->sym, Hc, lx, lyc);
  /*
   * C+=R*C
   */
  H2l2_HtimesVbig(H->sym, H, lxc, lyc);

  /*
   * transform long vector to short vector
   */
  proj_restr_et(H->disc, ly, y);

  free(lx);
  free(ly);

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *
 */
int H2l2_HtimesVbigMaxwell(ct_root *H, double *x1, double *y1) {
  int N = H->disc->mesh->nf * H->disc->a_bs;
  double *x2; /* source big vector */
  double *y2; /* target big vector */
  ct_root Hul[2] = {H[0], H[3]};
  ct_root Hur[2] = {H[1], H[4]};
  ct_root Hll[2] = {H[1], H[4]};
  ct_root Hlr[2] = {H[2], H[5]};

  /* this is a bit of a hack, but it is only a struct which is copied, not the
   * tree behind it */
  Hll[0].sym = 'T';
  // Hll[0].Tmom_right = Hur[0].Tmom_left;
  // Hll[0].Tmom_left = Hur[0].Tmom_right;
  Hll[1].sym = 'T';
  // Hll[1].Tmom_right = Hur[1].Tmom_left;
  // Hll[1].Tmom_left = Hur[1].Tmom_right;

  x2 = x1 + 2 * N;
  y2 = y1 + 2 * N;

  H2l2_HtimesVbigComplex(Hul, x1, y1);
  H2l2_HtimesVbigComplex(Hur, x2, y1);
  H2l2_HtimesVbigComplex(Hll, x1, y2);
  H2l2_HtimesVbigComplex(Hlr, x2, y2);

  return 0;
}

/**
 *  \brief      Matrix-vector multiplication V*x.
 *
 *  \param[in]     trans          'N' or 'T'.
 *  \param[in]     H              Hmatrix for V.
 *  \param[in]     x              Vector to which the operator is applied.
 *  \param[in,out] y              Vector to which the result will be added.
 *
 */
int H2l2_HtimesVsmallMaxwell(ct_root *H, double *x, double *y) {
  int N = H->disc->mesh->nf * H->disc->a_bs;
  double *lx1; /* source big vector */
  double *ly1; /* target big vector */

  /*
   * transform short vectors to long vectors
   */
  lx1 = (double *)calloc(4 * N, sizeof(double));
  ly1 = (double *)calloc(4 * N, sizeof(double));
  proj_distr_et(H->disc, x, lx1);

  H2l2_HtimesVbigMaxwell(H, lx1, ly1);

  /*
   * transform long vector to short vector
   */
  proj_restr_et(H->disc, ly1, y);

  free(lx1);
  free(ly1);

  return 0;
}
}  // namespace Bembel