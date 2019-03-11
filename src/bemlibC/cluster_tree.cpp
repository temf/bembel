// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "cluster_tree.h"

#include "tree_leafs.h"

namespace Bembel {
/**
 *  \brief         Checks the admissibility condition for two clusters.
 *
 *  \param[in]     C1             Cluster 1.
 *  \param[in]     C2             Cluster 2.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 *  \return        2, if min_bsize is reached, 0, if clusters have similarities
 *                 and 1 if product cluster can be low rank approximated.
 *
 */
int compare_cluster(et_node *C1, et_node *C2, int M, hmatrixsettings *hmatset) {
  double dist = sqrt((C1->midpoint.x - C2->midpoint.x) *
                         (C1->midpoint.x - C2->midpoint.x) +
                     (C1->midpoint.y - C2->midpoint.y) *
                         (C1->midpoint.y - C2->midpoint.y) +
                     (C1->midpoint.z - C2->midpoint.z) *
                         (C1->midpoint.z - C2->midpoint.z)) -
                C1->radius - C2->radius;
  int r;

  if (C1 == C2 || dist < 0 ||
      (C1->radius < C2->radius ? C2->radius : C1->radius) / dist >=
          hmatset->eta)
    r = 0;
  else
    r = 1;

  if (r == 0 && M - C1->level <= hmatset->min_bsize) r = 2;

  return r;
}

/**
 *  \brief         For two given clusters and this function generates the
 *                 corresponding product cluster subtree.
 *
 *  \param[in]     P              Point list.
 *  \param[in,out] pN1            Target node below diagonal.
 *  \param[in,out] pN2            Target node above diagonal.
 *  \param[in]     C1             Cluster 1.
 *  \param[in]     C2             Cluster 2.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 *  The function calls compare_cluster() to check the admissibility between the
 *  clusters C1 and C2.
 *
 *  If they are admissible, then the matrix block will be initialized as an
 *  rk-matrix using init_rkmat(). If C1 and C2 are inadmissible and \link
 *  min_bsize \endlink is reached, then the matrix block will be initialized
 *  as a full matrix using init_fmat().
 *
 *  If C1 and C2 are inadmissible and \link min_bsize \endlink is not reached,
 *  then 4*4 son clusters are assigned to pN and this function is called
 *  recursively.
 *
 */
int append_subtree(hmatrixfactory *hmatfac, vector3 *P, ct_node **pN1,
                   ct_node **pN2, et_node *pC1, et_node *pC2, int M) {
  int cc = 0;
  int i = 0;
  int j = 0;
  int ell;
  pdeproblem *pde = hmatfac->disc->pde;
  ct_node *(pn1[pde->nct]);
  ct_node *(pn2[pde->nct]);
  ct_rkmat *(prk1[pde->nct]);
  ct_rkmat *(prk2[pde->nct]);
  ct_fmat *(pf1[pde->nct]);
  ct_fmat *(pf2[pde->nct]);

  cc = compare_cluster(pC1, pC2, M, hmatfac->hmatset);
  if (cc) {
    /*
     * lowrank block applicable
     */
    if (1 == cc) {
      for (ell = 0; ell < pde->nct; ++ell) {
        pN1[ell]->rkmat = (ct_rkmat *)calloc(1, sizeof(ct_rkmat));
        prk1[ell] = pN1[ell]->rkmat;
        if (pde->symflags[ell] == 'N') {
          pN2[ell]->rkmat = (ct_rkmat *)calloc(1, sizeof(ct_rkmat));
          prk2[ell] = pN2[ell]->rkmat;
        }
      }
#pragma omp task
      init_rkmat(hmatfac, P, prk1, prk2, pC1, pC2, M);
    }
    /*
     * full matrix block applicable
     */
    else {
      /*
       * Check if we have to initialize a full matrix on the diagonal
       * or two matrices away from the diagonal
       */
      for (ell = 0; ell < pde->nct; ++ell) {
        pN1[ell]->fmat = (ct_fmat *)calloc(1, sizeof(ct_fmat));
        pf1[ell] = pN1[ell]->fmat;
      }
      if (pC1 != pC2) {
        for (ell = 0; ell < pde->nct; ++ell) {
          if (pde->symflags[ell] == 'N') {
            pN2[ell]->fmat = (ct_fmat *)calloc(1, sizeof(ct_fmat));
            pf2[ell] = pN2[ell]->fmat;
          }
        }
#pragma omp task
        init_fmat(hmatfac, P, pf1, pf2, pC1, pC2, M);
      } else {
#pragma omp task
        init_sym_fmat(hmatfac, P, pf1, pC1, pC2, M);
      }
    }
  }
  /*
   * if none of the above cases apply, there are children to handle
   */
  else {
    /*
     * pC1 == pC2 means, that the two clusters are on the diagonal
     */
    if (pC1 == pC2) {
      for (ell = 0; ell < pde->nct; ++ell) {
        pN1[ell]->nch = 16;
        pN1[ell]->pch = (ct_node **)calloc(1, 16 * sizeof(ct_node *));
      }
      for (i = 0; i < 4; ++i) {
        /*
         * initialization of the product cluster on the diagonal
         */
        for (ell = 0; ell < pde->nct; ++ell) {
          pN1[ell]->pch[i * 4 + i] = (ct_node *)calloc(1, sizeof(ct_node));
          pN1[ell]->pch[i * 4 + i]->cluster1 = i;
          pN1[ell]->pch[i * 4 + i]->cluster2 = i;
          pN1[ell]->pch[i * 4 + i]->level = pN1[ell]->level + 1;
          pN1[ell]->pch[i * 4 + i]->pdad = pN1[ell];
          pn1[ell] = pN1[ell]->pch[i * 4 + i];
        }
        append_subtree(hmatfac, P, pn1, pn1, pC1->son[i], pC2->son[i], M);
        for (j = 0; j < i; ++j) {
          /*
           * initialization of the product cluster away from
           * diagonal
           */
          for (ell = 0; ell < pde->nct; ++ell) {
            pN1[ell]->pch[i * 4 + j] = (ct_node *)calloc(1, sizeof(ct_node));
            pN1[ell]->pch[i * 4 + j]->cluster1 = i;
            pN1[ell]->pch[i * 4 + j]->cluster2 = j;
            pN1[ell]->pch[i * 4 + j]->level = pN1[ell]->level + 1;
            pN1[ell]->pch[i * 4 + j]->pdad = pN1[ell];
            pn1[ell] = pN1[ell]->pch[i * 4 + j];
            if (pde->symflags[ell] == 'N') {
              pN1[ell]->pch[j * 4 + i] = (ct_node *)calloc(1, sizeof(ct_node));
              pN1[ell]->pch[j * 4 + i]->cluster1 = j;
              pN1[ell]->pch[j * 4 + i]->cluster2 = i;
              pN1[ell]->pch[j * 4 + i]->level = pN1[ell]->level + 1;
              pN1[ell]->pch[j * 4 + i]->pdad = pN1[ell];
              pn2[ell] = pN1[ell]->pch[j * 4 + i];
            }
          }
          append_subtree(hmatfac, P, pn1, pn2, pC1->son[i], pC2->son[j], M);
        }
      }
    } else {
      for (ell = 0; ell < pde->nct; ++ell) {
        pN1[ell]->nch = 16;
        pN1[ell]->pch = (ct_node **)calloc(1, 16 * sizeof(ct_node *));
        if (pde->symflags[ell] == 'N') {
          pN2[ell]->nch = 16;
          pN2[ell]->pch = (ct_node **)calloc(1, 16 * sizeof(ct_node *));
        }
      }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j) {
          for (ell = 0; ell < pde->nct; ++ell) {
            pN1[ell]->pch[i * 4 + j] = (ct_node *)calloc(1, sizeof(ct_node));
            pN1[ell]->pch[i * 4 + j]->cluster1 = i;
            pN1[ell]->pch[i * 4 + j]->cluster2 = j;
            pN1[ell]->pch[i * 4 + j]->level = pN1[ell]->level + 1;
            pN1[ell]->pch[i * 4 + j]->pdad = pN1[ell];
            pn1[ell] = pN1[ell]->pch[i * 4 + j];
            if (pde->symflags[ell] == 'N') {
              pN2[ell]->pch[j * 4 + i] = (ct_node *)calloc(1, sizeof(ct_node));
              pN2[ell]->pch[j * 4 + i]->cluster1 = j;
              pN2[ell]->pch[j * 4 + i]->cluster2 = i;
              pN2[ell]->pch[j * 4 + i]->level = pN2[ell]->level + 1;
              pN2[ell]->pch[j * 4 + i]->pdad = pN2[ell];
              pn2[ell] = pN2[ell]->pch[j * 4 + i];
            }
          }
          append_subtree(hmatfac, P, pn1, pn2, pC1->son[i], pC2->son[j], M);
        }
    }
  }
  return 0;
}

/**
 *  \brief         Initializes the product cluster tree aka H-matrix structure.
 *
 *  \param[in]     P              Point list.
 *  \param[in,out] H              H[0] root of Single, H[1] root of Double
 *  \param[in]     E              Root of the corresponding element tree.
 *  \param[in]     p              Number of patches.
 *  \param[in]     M              Number of Levels of the cluster tree.
 *
 *  Generates p*p patches as children of H and applies append_subtree() to each
 *  of them.
 *
 */
void init_cluster_tree(hmatrixfactory *hmatfac, ct_root **H) {
  int i = 0;
  int j = 0;
  int ell = 0;
  void (**Tmom_left_phis)(double *, double, double);
  void (**Tmom_right_phis)(double *, double, double);
  discretization *disc = hmatfac->disc;
  pdeproblem *pde = disc->pde;
  meshdata *mesh = disc->mesh;
  const geometry &geom = *mesh->geom;
  const int p = geom.size();
  ct_node *(pNf[pde->nct]);
  ct_node *(pNs1[pde->nct]);
  ct_node *(pNs2[pde->nct]);

  init_Gauss_Square(&hmatfac->Q, g_max + 1);
  pde->init_randwerte(&hmatfac->RW, &hmatfac->Q[disc->g_far], geom,
                      1 << mesh->M);

  for (ell = 0; ell < pde->nct; ++ell) {
    H[0][ell].p = p;
    H[0][ell].M = mesh->M;
    H[0][ell].sym = pde->symflags[ell];
    H[0][ell].rank = hmatfac->hmatset->np_max * hmatfac->hmatset->np_max *
                     hmatfac->disc->pde->np_max_fac;
    H[0][ell].disc = hmatfac->disc;
    H[0][ell].hmatset = hmatfac->hmatset;
    init_transfermatrices(disc, hmatfac->hmatset, &H[0][ell].Ttr);
    Tmom_left_phis = pde->Tmom_left_phi(disc, ell);
    Tmom_right_phis = pde->Tmom_right_phi(disc, ell);
    init_momentmatrices(disc, hmatfac->hmatset, &H[0][ell].Tmom_left,
                        Tmom_left_phis);
    init_momentmatrices(disc, hmatfac->hmatset, &H[0][ell].Tmom_right,
                        Tmom_right_phis);
    free(Tmom_left_phis);
    free(Tmom_right_phis);
    pNf[ell] = H[0][ell].root = (ct_node *)calloc(1, sizeof(ct_node));
    pNf[ell]->nch = p * p;
    pNf[ell]->pch = (ct_node **)calloc(p * p, sizeof(ct_node *));
    assert(pde->symflags[ell] == 'N' || pde->symflags[ell] == 'S');
  }

#pragma omp parallel default(shared) firstprivate(i, j, ell, pNf, pNs1, pNs2, H)
  {
#pragma omp single
    {
      for (i = 0; i < p; ++i) {
        for (ell = 0; ell < pde->nct; ++ell) {
          pNs1[ell] = pNf[ell]->pch[i * p + i] =
              (ct_node *)calloc(1, sizeof(ct_node));
          pNs1[ell]->cluster1 = i;
          pNs1[ell]->cluster2 = i;
          pNs1[ell]->level = 0;
          pNs1[ell]->pdad = pNf[ell];
        }
        append_subtree(hmatfac, mesh->P, pNs1, pNs1, mesh->E.patch[i],
                       mesh->E.patch[i], mesh->M);
        for (j = 0; j < i; ++j) {
          for (ell = 0; ell < pde->nct; ++ell) {
            pNs1[ell] = pNf[ell]->pch[i * p + j] =
                (ct_node *)calloc(1, sizeof(ct_node));
            pNs1[ell]->cluster1 = i;
            pNs1[ell]->cluster2 = j;
            pNs1[ell]->level = 0;
            pNs1[ell]->pdad = pNf[ell];
            if (pde->symflags[ell] == 'N') {
              pNs2[ell] = pNf[ell]->pch[j * p + i] =
                  (ct_node *)calloc(1, sizeof(ct_node));
              pNs2[ell]->cluster1 = j;
              pNs2[ell]->cluster2 = i;
              pNs2[ell]->level = 0;
              pNs2[ell]->pdad = pNf[ell];
            }
          }
          append_subtree(hmatfac, mesh->P, pNs1, pNs2, mesh->E.patch[i],
                         mesh->E.patch[j], mesh->M);
        }
      }
    }
  }

  pde->free_randwerte(&hmatfac->RW, 1 << mesh->M, p);
  free_Gauss_Square(&hmatfac->Q, g_max + 1);

  pde->postproc(hmatfac, H);

  // printf("%26s Full matrices: %d\n", "", hmatfac->assemfmats);
  // printf("%26s Full symmetric matrices: %d\n", "", hmatfac->assemsfmats);
  // printf("%26s Low-rank matrices: %d\n", "", hmatfac->assemrkmats);

  return;
}

/**
 *  \brief         Writes two files that can be visualized with matlab.
 *
 *  \param[in]     x0             Enter zero.
 *  \param[in]     y0             Enter zero.
 *  \param[in]     pN             Root node of the Hmatrix (e.g. H.root).
 *  \param[in,out] f              File handle for the cluster-file.
 *  \param[in,out] g              File handle for the color-file.
 *
 *  The matlab routine is called print_Hmatrix.m and is in the matlab folder.
 *
 *  \note          You may want to call fast_print_cluster_tree() instead.
 *
 */
int print_cluster_tree(double x0, double y0, ct_node *pN, FILE *f, FILE *g) {
  double x = 0;
  double y = 0;
  int i = 0;

  for (i = 0; i < pN->nch; ++i) {
    if (pN->pch[i]) {
      x = x0 + 1. / (1 << pN->pch[i]->level * 2) * pN->pch[i]->cluster1;
      y = y0 + 1. / (1 << pN->pch[i]->level * 2) * pN->pch[i]->cluster2;
      fprintf(f, "%15.10f\t %15.10f\n", x, y);
      x = x0 + 1. / (1 << pN->pch[i]->level * 2) * (pN->pch[i]->cluster1 + 1);
      y = y0 + 1. / (1 << pN->pch[i]->level * 2) * pN->pch[i]->cluster2;
      fprintf(f, "%15.10f\t %15.10f\n", x, y);
      x = x0 + 1. / (1 << pN->pch[i]->level * 2) * (pN->pch[i]->cluster1 + 1);
      y = y0 + 1. / (1 << pN->pch[i]->level * 2) * (pN->pch[i]->cluster2 + 1);
      fprintf(f, "%15.10f\t %15.10f\n", x, y);
      x = x0 + 1. / (1 << pN->pch[i]->level * 2) * (pN->pch[i]->cluster1);
      y = y0 + 1. / (1 << pN->pch[i]->level * 2) * (pN->pch[i]->cluster2 + 1);
      fprintf(f, "%15.10f\t %15.10f\n", x, y);
      if (pN->pch[i]->fmat)
        fprintf(g, "%5.4f\t %5.4f\t %5.4f\n", 1., 0.,
                1. * pN->pch[i]->fmat->bs);
      else if (pN->pch[i]->rkmat) {
        if (pN->pch[i]->rkmat->ker)
          fprintf(g, "%5.4f\t %5.4f\t %5.4f\n",
                  1. * pN->pch[i]->rkmat->ker->rk / pN->pch[i]->rkmat->bs,
                  1. * pN->pch[i]->rkmat->ker->rk, 1. * pN->pch[i]->rkmat->bs);

        else
          fprintf(g, "%5.4f\t %5.4f\t %5.4f\n",
                  1. * pN->pch[i]->rkmat->k / pN->pch[i]->rkmat->bs,
                  1. * pN->pch[i]->rkmat->k, 1. * pN->pch[i]->rkmat->bs);
      } else
        fprintf(g, "%5.4f\t %5.4f\t %5.4f\n", 0., 0., 0.);
    }
  }

  for (i = 0; i < pN->nch; ++i)
    if (pN->pch[i])
      if (pN->pch[i]->pch) {
        x = x0 + 1. / (1 << pN->pch[i]->level * 2) * pN->pch[i]->cluster1;
        y = y0 + 1. / (1 << pN->pch[i]->level * 2) * pN->pch[i]->cluster2;
        print_cluster_tree(x, y, pN->pch[i], f, g);
      }

  return 0;
}

/**
 *  \brief         A fast way to use print_cluster_tree().
 *
 *  \param[in]     pN             E.g. H.root.
 *  \param[in]     namesuffix     The filenames will be this suffix with an
 *                                attached "_clusters" and "_colors"
 */
void fast_print_cluster_tree(ct_node *pN, char *namesuffix) {
  FILE *f;
  FILE *g;
  char fname[100];

  sprintf(fname, "%s_clusters", namesuffix);
  f = fopen(fname, "w");
  sprintf(fname, "%s_colors", namesuffix);
  g = fopen(fname, "w");

  print_cluster_tree(0, 0, pN, f, g);

  fclose(f);
  fclose(g);

  return;
}

/**
 *  \brief         Frees any tree structure similar to the one generated by
 *                 init_cluster_tree().
 *
 *  \param[in]     pN             E.g. H.root.
 *  \param[in]     p              Dummy argument.
 *
 */
int free_cluster_tree_rec(ct_node *pN, int p)
// ct_node *pN;
// int p;
{
  int i = 0;

  if (pN->fmat) {
    free_fmat(pN->fmat);
    free(pN->fmat);
    pN->fmat = NULL;
  } else if (pN->rkmat) {
    free_rkmat(pN->rkmat);
    free(pN->rkmat);
    pN->rkmat = NULL;
  } else if (pN->pch) {
    for (i = 0; i < pN->nch; ++i) free_cluster_tree_rec(pN->pch[i], p);
    free(pN->pch);
    pN->pch = NULL;
  } else {
    /*
     * nothing to do ?!??
     */
  }
  free(pN);
  pN = NULL;

  return 0;
}

/**
 *  \brief         Frees any tree structure similar to the one generated by
 *                 init_sym_cluster_tree().
 *
 *  \param[in]     pN             E.g. H.root.
 *  \param[in]     p              Number of patches.
 *
 */
int free_sym_cluster_tree_rec(ct_node *pN, int p) {
  int i = 0;
  int j = 0;

  if (pN->pch) {
    for (i = 0; i < p; ++i) {
      free_sym_cluster_tree_rec(pN->pch[i * p + i], 4);

      for (j = 0; j < i; ++j) free_cluster_tree_rec(pN->pch[i * p + j], 4);
    }
    free(pN->pch);
  }

  if (pN->fmat) {
    free_fmat(pN->fmat);
    free(pN->fmat);
  }

  free(pN);

  return 0;
}

int free_cluster_tree(ct_root *H) {
  free_momentmatrices(H->disc, H->hmatset, H->Tmom_left);
  free_momentmatrices(H->disc, H->hmatset, H->Tmom_right);
  free_transfermatrices(H->disc, H->hmatset, H->Ttr);
  if (H->sym == 'N' || H->sym == 'T')
    free_cluster_tree_rec(H->root, H->p);
  else
    free_sym_cluster_tree_rec(H->root, H->p);

  return 1;
}
}  // namespace Bembel
