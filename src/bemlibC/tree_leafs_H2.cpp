// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "tree_leafs.h"

namespace Bembel {
/**
 *  \brief         Initializes the rkmatrices. Uses the symmetric structure
 *                 of the cluster tree to increase computation speed.
 *
 *  \param[in]     P              Point list.
 *  \param[in,out] prk1           Rkmatrix below the diagonal.
 *  \param[in,out] prk2           Rkmatrix above the diagonal.
 *  \param[in]     pC1            Cluster 1 for matrix below diagonal.
 *  \param[in]     pC2            Cluster 2 for matrix below diagonal.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 */
int init_rkmat(hmatrixfactory *hmatfac, vector3 *P, ct_rkmat **prk1,
               ct_rkmat **prk2, et_node *pC1, et_node *pC2, int M) {
  int ell = 0;
  int lind1 = 0;
  int lind2 = 0;
  int cind1 = 0;
  int cind2 = 0;
  int np_max = hmatfac->hmatset->np_max;
  pdeproblem *pde = hmatfac->disc->pde;
  int np_max_fac = pde->np_max_fac;
  double **(A1[pde->nct]);
  double **(A2[pde->nct]);
  et_node *pp = NULL;

#pragma omp atomic
  hmatfac->assemrkmats += 1;

  pp = pC1;
  while (pp->son[0]) pp = pp->son[0];
  lind1 = pp->number;
  pp = pC1;
  while (pp->son[3]) pp = pp->son[3];
  lind2 = pp->number;
  pp = pC2;
  while (pp->son[0]) pp = pp->son[0];
  cind1 = pp->number;
  pp = pC2;
  while (pp->son[3]) pp = pp->son[3];
  cind2 = pp->number;
  pp -= cind2;

  for (ell = 0; ell < pde->nct; ++ell) {
    prk1[ell]->lind1 = lind1;
    prk1[ell]->cind1 = cind1;
    prk1[ell]->bs = lind2 - lind1 + 1;
    prk1[ell]->k = np_max * np_max * np_max_fac;
    prk1[ell]->ker = (ct_ker *)calloc(1, sizeof(ct_ker));
    A1[ell] = &(prk1[ell]->ker->A);
    if (pde->symflags[ell] == 'N') {
      prk2[ell]->lind1 = lind1;
      prk2[ell]->cind1 = cind1;
      prk2[ell]->bs = lind2 - lind1 + 1;
      prk2[ell]->k = np_max * np_max * np_max_fac;
      prk2[ell]->ker = (ct_ker *)calloc(1, sizeof(ct_ker));
      A2[ell] = &(prk2[ell]->ker->A);
    }
  }

  interpol_kernel(hmatfac, A1, A2, pC1, pC2, np_max);

  return 0;
}

/**
 *  \brief         Initializes the full matrices away from the diagonal
 *                 in the Hmatrix.
 *
 *  \param[in]     P              Point list.
 *  \param[in,out] pf1            Full matrix below the diagonal.
 *  \param[in,out] pf2            Full matrix above the diagonal.
 *  \param[in]     pC1            Cluster 1 for matrix below diagonal.
 *  \param[in]     pC2            Cluster 2 for matrix below diagonal.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 */
int init_fmat(hmatrixfactory *hmatfac, vector3 *P, ct_fmat **pf1, ct_fmat **pf2,
              et_node *pC1, et_node *pC2, int M) {
  int ell = 0;
  int lind1 = 0;
  int lind2 = 0;
  int cind1 = 0;
  int cind2 = 0;
  int nl = 0;
  int nl2 = 0;
  int i = 0;
  int j = 0;
  int a_bs = hmatfac->disc->a_bs;
  int a_bs2 = hmatfac->disc->a_bs2;
  pdeproblem *pde = hmatfac->disc->pde;
  double *dest;
  double *A1[pde->nct];
  double *A2[pde->nct];
  et_node *pp = NULL;

#pragma omp atomic
  hmatfac->assemfmats += 1;

  pp = pC1;
  while (pp->son[0]) pp = pp->son[0];
  lind1 = pp->number;
  pp = pC1;
  while (pp->son[3]) pp = pp->son[3];
  lind2 = pp->number;
  pp = pC2;
  while (pp->son[0]) pp = pp->son[0];
  cind1 = pp->number;
  pp = pC2;
  while (pp->son[3]) pp = pp->son[3];
  cind2 = pp->number;
  pp -= cind2;
  nl = lind2 - lind1 + 1;
  nl2 = nl * nl;

  for (ell = 0; ell < pde->nct; ++ell) {
    pf1[ell]->bs = a_bs * nl;
    pf1[ell]->A = A1[ell] = (double *)malloc(a_bs2 * nl2 * sizeof(double));
    if (pde->symflags[ell] == 'N') {
      pf2[ell]->bs = a_bs * nl;
      pf2[ell]->A = A2[ell] = (double *)malloc(a_bs2 * nl2 * sizeof(double));
    }
  }

  dest = (double *)malloc(pde->quadrature_bufsize * a_bs2 * sizeof(double));

  for (i = 0; i < nl; ++i)
    for (j = 0; j < nl; ++j) {
      memset(dest, 0, pde->quadrature_bufsize * a_bs2 * sizeof(double));
      get_stiff(hmatfac, P, dest, pp, M, i + lind1, j + cind1);
      pde->insert_fmat(A1, A2, i, j, nl, a_bs, a_bs2, dest);
    }

  free(dest);

  return 0;
}

/**
 *  \brief         Initializes the full matrices on the diagonal
 *                 of the Hmatrix.
 *
 *  \param[in]     P              Point list.
 *  \param[in,out] pf1            Full matrix.
 *  \param[in]     pC1            Cluster 1.
 *  \param[in]     pC2            Cluster 2.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 */
int init_sym_fmat(hmatrixfactory *hmatfac, vector3 *P, ct_fmat **pf1,
                  et_node *pC1, et_node *pC2, int M) {
  int lind1 = 0;
  int lind2 = 0;
  int cind1 = 0;
  int cind2 = 0;
  int nl = 0;
  int nl2 = 0;
  int ell = 0;
  int i = 0;
  int j = 0;
  int a_bs = hmatfac->disc->a_bs;
  int a_bs2 = hmatfac->disc->a_bs2;
  pdeproblem *pde = hmatfac->disc->pde;
  double *dest;
  double *A[pde->nct];
  et_node *pp = NULL;

#pragma omp atomic
  hmatfac->assemsfmats += 1;

  pp = pC1;
  while (pp->son[0]) pp = pp->son[0];
  lind1 = pp->number;
  pp = pC1;
  while (pp->son[3]) pp = pp->son[3];
  lind2 = pp->number;
  pp = pC2;
  while (pp->son[0]) pp = pp->son[0];
  cind1 = pp->number;
  pp = pC2;
  while (pp->son[3]) pp = pp->son[3];
  cind2 = pp->number;
  pp -= cind2;
  nl = lind2 - lind1 + 1;
  nl2 = nl * nl;

  for (ell = 0; ell < pde->nct; ++ell) {
    pf1[ell]->bs = a_bs * nl;
    pf1[ell]->A = A[ell] = (double *)malloc(a_bs2 * nl2 * sizeof(double));
  }

  dest = (double *)malloc(pde->quadrature_bufsize * a_bs2 * sizeof(double));

  for (i = 0; i < nl; ++i) {
    memset(dest, 0, pde->quadrature_bufsize * a_bs2 * sizeof(double));
    get_stiff(hmatfac, P, dest, pp, M, i + lind1, i + cind1);
    pde->insert_sym_fmat_diag(A, i, nl, a_bs, a_bs2, dest);
    for (j = 0; j < i; ++j) {
      memset(dest, 0, pde->quadrature_bufsize * a_bs2 * sizeof(double));
      get_stiff(hmatfac, P, dest, pp, M, i + lind1, j + cind1);
      pde->insert_sym_fmat_offdiag(A, i, j, nl, a_bs, a_bs2, dest);
    }
  }

  free(dest);

  return 0;
}

/**
 *  \brief         Frees an rkmat created by init_rkmat().
 *
 *  \param[in]     prk            rkmatrix to free.
 *
 */
int free_rkmat(ct_rkmat *prk) {
  prk->l = NULL;
  prk->r = NULL;
  free(prk->ker->A);
  free(prk->ker);
  prk->ker = NULL;

  return 0;
}

/**
 *  \brief         Frees a full matrix created by init_fmat().
 *
 *  \param[in]     pf             Fmat to free.
 *
 */
int free_fmat(ct_fmat *pf) {
  free(pf->A);
  return 0;
}
}  // namespace Bembel