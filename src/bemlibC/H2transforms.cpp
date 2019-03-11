// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "H2transforms.h"

namespace Bembel {
void new_H2tree(H2tree *t, int depth, int nsons, int rank) {
  int i;
  int j;
  H2tree *s;
  H2tree *sl;

  /*
   * assign stuff to root
   */
  t->nsons = nsons;
  t->depth = depth;
  t->data = (double *)calloc(rank, sizeof(double));

  /*
   * get storage for sons
   */
  s = t;
  for (i = depth; i > 0; --i) {
    s->sons =
        (H2tree *)calloc(nsons * (1 << (2 * (depth - i))), sizeof(H2tree));
    for (j = 0; j < nsons * (1 << (2 * (depth - i))); ++j)
      s->sons[j].data = (double *)calloc(rank, sizeof(double));
    s = &s->sons[0];
  }

  if (depth > 1) {
    s = &t->sons[0];
    sl = &s->sons[0];
    for (i = depth; i > 1; --i) {
      for (j = 0; j < nsons * (1 << (2 * (depth - i))); ++j) {
        s[j].nsons = 4;
        s[j].depth = i - 1;
      }
      for (j = 1; j < nsons * (1 << (2 * (depth - i))); ++j)
        s[j].sons = &sl[4 * j];
      s = sl;
      sl = &s->sons[0];
    }
  }

  return;
}

void free_H2tree(H2tree *t) {
  int i;
  int nsons;
  int depth;
  H2tree *s;

  free(t->data);
  nsons = t->nsons;
  depth = t->depth;
  t = &t->sons[0];

  for (; depth > 0; --depth) {
    for (i = 0; i < nsons; ++i) free(t[i].data);
    s = &t->sons[0];
    free(t);
    t = s;
    nsons *= 4;
  }

  return;
}

void add_H2tree(H2tree *t, H2tree *a, int rank) {
  int i;
  mydxpy(rank, a->data, t->data);
  for (i = 0; i < t->nsons; ++i) add_H2tree(&t->sons[i], &a->sons[i], rank);

  return;
}

/**
 *  \brief         Forward transform for H2-matrices
 *
 *  \param[in]     Tmom           Moment matrices on loweset level
 *  \param[in]     Ttr            Transfer matrices
 *  \param[in]     minlvl         Block size of smallest matrix block is
 *                                2^minlvl
 *  \param[in]     N             Length of nx
 *  \param[in]     M              Discretization level
 *  \param[in]     x              Vector to transform
 *  \param[out]    transx         H2 output
 *
 */
int ForwardTransform(ct_root *Hmat, double *x, H2tree *transx) {
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int aktn = 0;
  int minlvl = Hmat->hmatset->min_bsize;
  int a_bs = Hmat->disc->a_bs;
  int N = a_bs * Hmat->disc->mesh->nf;
  int rank = Hmat->rank;
  int np_max_fac = Hmat->disc->pde->np_max_fac;
  int np2 = rank / np_max_fac;
  int M = Hmat->disc->mesh->M;
  double **Tmom = Hmat->sym == 'T' ? Hmat->Tmom_left : Hmat->Tmom_right;
  double ***Ttr = Hmat->Ttr;
  H2tree *mytransx;
  H2tree *mytransx2;

  if (!Tmom) return 0;
  // assert(!(rank > (1 << 2 * minlvl) * a_bs));

  /*
   * find lowest level in H2tree
   */
  mytransx = transx;
  while (mytransx->depth > 0) mytransx = &mytransx->sons[0];

  /*
   * apply Moment Tmom to vector x and store everything in transx[0]
   */
  aktn = a_bs * (1 << (2 * (minlvl)));
  for (i = 0; i < N / aktn; ++i)
    for (j = 0; j < rank; ++j)
      mytransx[i].data[j] = myddot(aktn, Tmom[j], &x[i * aktn]);

  /*
   * successively determine Li to the vector x by use of the transforms Ttr
   */
  for (i = 1; i < M - minlvl + 1; ++i) {
    mytransx = transx;
    while (mytransx->depth > i) mytransx = &mytransx->sons[0];
    mytransx2 = &mytransx->sons[0];
    aktn = a_bs * (1 << (2 * (i + minlvl)));
    for (j = 0; j < N / aktn; ++j)
      for (k = 0; k < np_max_fac; ++k)
        for (l = 0; l < np2; ++l) {
          mytransx[j].data[k * np2 + l] +=
              myddot(np2, Ttr[0][l], mytransx2[4 * j].data + k * np2) +
              myddot(np2, Ttr[2][l], mytransx2[4 * j + 1].data + k * np2) +
              myddot(np2, Ttr[3][l], mytransx2[4 * j + 2].data + k * np2) +
              myddot(np2, Ttr[1][l], mytransx2[4 * j + 3].data + k * np2);
        }
  }

  return 0;
}

/**
 *  \brief         Backward transform for H2-matrices
 *
 *  \param[in]     Tmom           Moment matrices on loweset level
 *  \param[in]     Ttr            Transfer matrices
 *  \param[in]     minlvl         Block size of smallest matrix block is
 *                                2^minlvl
 *  \param[in]     N             Length of nx
 *  \param[in]     M              Discretization level
 *  \param[out]    x              Destination vector
 *  \param[in]    transx         H2 input
 *
 */
int BackwardTransform(ct_root *Hmat, double *x, H2tree *transx) {
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int aktn = 0;
  int minlvl = Hmat->hmatset->min_bsize;
  int a_bs = Hmat->disc->a_bs;
  int N = a_bs * Hmat->disc->mesh->nf;
  int M = Hmat->disc->mesh->M;
  int rank = Hmat->rank;
  int np_max_fac = Hmat->disc->pde->np_max_fac;
  int np2 = rank / np_max_fac;
  double **Tmom = Hmat->sym == 'T' ? Hmat->Tmom_right : Hmat->Tmom_left;
  double ***Ttr = Hmat->Ttr;
  H2tree *mytransx;

  if (!Tmom) return 0;
  // assert(!(rank > (1 << 2 * minlvl) * a_bs));

  mytransx = &transx->sons[0];
  for (i = M - minlvl; i > 0; --i) {
    aktn = a_bs * (1 << (2 * (i + minlvl)));
    for (j = 0; j < N / aktn; ++j)
      for (k = 0; k < np_max_fac; ++k)
        for (l = 0; l < np2; ++l) {
          mydaxpy(np2, mytransx[j].data[k * np2 + l], Ttr[0][l],
                  mytransx[j].sons[0].data + k * np2);
          mydaxpy(np2, mytransx[j].data[k * np2 + l], Ttr[2][l],
                  mytransx[j].sons[1].data + k * np2);
          mydaxpy(np2, mytransx[j].data[k * np2 + l], Ttr[3][l],
                  mytransx[j].sons[2].data + k * np2);
          mydaxpy(np2, mytransx[j].data[k * np2 + l], Ttr[1][l],
                  mytransx[j].sons[3].data + k * np2);
        }
    mytransx = &mytransx[0].sons[0];
  }

  aktn = a_bs * (1 << (2 * (minlvl)));
  for (i = 0; i < N / aktn; ++i)
    for (j = 0; j < rank; ++j)
      mydaxpy(aktn, mytransx[i].data[j], Tmom[j], &x[i * aktn]);

  return 0;
}
}  // namespace Bembel