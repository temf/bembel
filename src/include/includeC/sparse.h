// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_SPARSE__
#define __BEMBEL_C_SPARSE__

/**************
 *  sparsemats.h  *
 **************/

/*==================================================*
 *  Modul zur Rechnung mit duennbesetzten Matrizen  *
 *==================================================*/

#include <cstdlib>
#include <cstring>
namespace Bembel {

/*====================*
 *  Typendeklaration  *
 *====================*/

template <typename T>
struct sparsecore {
  int m, n;
  int *row_number, *max_row_number;
  int **index;
  T **value;
};

typedef sparsecore<double> sparse;
typedef sparsecore<int> sparsei;

/*==========================================*
 *  Suchalgorithmus gemaess binary-search:  *
 *  Liefert den Index des Elements mit      *
 *  Index >= dem gesuchten Index j.         *
 *==========================================*/

template <typename T>
int search_sparse(int *array, int rn, int j) noexcept {
  int low = 0;
  int high = rn;
  while (low < high) {
    int mid = (low + high) / 2;
    if (array[mid] < j)
      low = mid + 1;
    else if (array[mid] > j)
      high = mid;
    else
      return (mid);
  }
  return (low); /* low == high! */
}

/*==================================*
 *  memmove fuer die sparse-Matrix  *
 *==================================*/

template <typename T>
void memmove_sparse(sparsecore<T> *A, int i, int j) noexcept {
  int rn = A->row_number[i];

  if (rn == A->max_row_number[i]) { /* Anzahl Elemente in Zeile i muss erhoeht
                                       werden */
    int n_max = rn + 10;
    A->value[i] = (T *)realloc(A->value[i], n_max * sizeof(T));
    A->index[i] = (int *)realloc(A->index[i], (n_max + 1) * sizeof(int));
    A->max_row_number[i] = n_max;
  }

  memmove(&A->value[i][j + 1], &A->value[i][j], (rn - j) * sizeof(T));
  memmove(&A->index[i][j + 1], &A->index[i][j], (rn + 1 - j) * sizeof(int));
  A->row_number[i]++;
  return;
}

/*===================================*
 *  memmove fuer die Matrix-Pattern  *
 *===================================*/

template <typename T>
void memmove_pattern(sparsecore<T> *A, int i, int j) noexcept {
  int rn = A->row_number[i];

  if (rn == A->max_row_number[i]) { /* Anzahl Elemente in Zeile i muss erhoeht
                                       werden */
    int n_max = rn + 10;
    A->index[i] = (int *)realloc(A->index[i], (n_max + 1) * sizeof(int));
    A->max_row_number[i] = n_max;
  }

  memmove(&A->index[i][j + 1], &A->index[i][j], (rn + 1 - j) * sizeof(int));
  A->row_number[i]++;
  return;
}

/*=====================================*
 *  Initialisierung der sparse-Matrix  *
 *=====================================*/

template <typename T>
void init_sparse(sparsecore<T> *A, int m, int n, int n_max) noexcept {
  A->m = m;
  A->n = n;

  A->value = (T **)malloc(m * sizeof(T *));
  A->index = (int **)malloc(m * sizeof(int *));
  A->row_number = (int *)calloc(m, sizeof(int));
  A->max_row_number = (int *)malloc(m * sizeof(int));

  for (int i = 0; i < m; i++) { /* 1 Dummy-Element: Waechter- + Pufferelement */
    A->max_row_number[i] = n_max;
    A->value[i] = (T *)malloc(n_max * sizeof(T));
    A->index[i] = (int *)malloc((n_max + 1) * sizeof(int));
    A->index[i][0] = n; /* Waechterelement */
  }

  return;
}

/*======================================*
 *  Initialisierung der Matrix-Pattern  *
 *======================================*/

template <typename T>
void init_pattern(sparsecore<T> *A, int m, int n, int n_max) noexcept {
  A->m = m;
  A->n = n;

  A->value = (T **)malloc(m * sizeof(T *));
  A->index = (int **)malloc(m * sizeof(int *));
  A->row_number = (int *)calloc(m, sizeof(int));
  A->max_row_number = (int *)malloc(m * sizeof(int));

  for (int i = 0; i < m; i++) { /* 1 Dummy-Element: Waechter- + Pufferelement */
    A->max_row_number[i] = n_max;
    A->index[i] = (int *)malloc((n_max + 1) * sizeof(int));
    A->index[i][0] = n; /* Waechterelement */
  }

  return;
}

/*===============================*
 *  Freigeben der sparse-Matrix  *
 *===============================*/

template <typename T>
inline void free_sparse(sparsecore<T> *A) noexcept {
  for (int i = 0; i < A->m; i++) {
    free(A->value[i]);
    free(A->index[i]);
  }

  free(A->value);
  free(A->index);
  free(A->row_number);
  free(A->max_row_number);

  return;
}

/*==================================*
 *  Fuege A(i,j) den Pattern hinzu  *
 *==================================*/

template <typename T>
void set_pattern(sparsecore<T> *A, int i, int j) noexcept {
  int k = search_sparse(A->index[i], A->row_number[i], j);

  if (A->index[i][k] != j) /* Eintrag noch nicht vorhanden */
  {
    memmove_pattern(A, i, k);
    A->index[i][k] = j;
  }

  return;
}

/*===============*
 *  A(i,j) := z  *
 *===============*/

template <typename T>
void set_sparse(sparsecore<T> *A, int i, int j, T z) noexcept {
  int k = search_sparse<T>(A->index[i], A->row_number[i], j);

  if (A->index[i][k] != j) /* Eintrag noch nicht vorhanden */
  {
    memmove_sparse(A, i, k);
    A->index[i][k] = j;
  }

  A->value[i][k] = z;
  return;
}

/*===============*
 *  A(i,j) += z  *
 *===============*/

template <typename T>
void add_sparse(sparsecore<T> *A, int i, int j, T z) noexcept {
  if (z == 0) return;

  int k = search_sparse<T>(A->index[i], A->row_number[i], j);
  if (A->index[i][k] == j)
    A->value[i][k] += z; /* Eintrag schon vorhanden      */
  else                   /* Eintrag noch nicht vorhanden */
  {
    memmove_sparse(A, i, k);
    A->index[i][k] = j;
    A->value[i][k] = z;
  }
  return;
}

/*===============*
 *  z := A(i,j)  *
 *===============*/

template <typename T>
T get_sparse(sparsecore<T> *A, int i, int j) noexcept {
  int k = search_sparse<T>(A->index[i], A->row_number[i], j);
  if (A->index[i][k] == j)
    return (A->value[i][k]); /* Eintrag vorhanden       */
  else
    return (0); /* Eintrag nicht vorhanden */
}

/*====================*
 *  nz = non_zero(A)  *
 *====================*/

template <typename T>
int non_zero(sparsecore<T> *A) noexcept {
  int nz = 0;
  for (int i = 0; i < A->m; i++) nz += A->row_number[i];

  return (nz);
}

/*=====================*
 *  Garbage-Collector  *
 *=====================*/

template <typename T>
void garbage_collect(sparsecore<T> *A) noexcept {
  for (int i = 0; i < A->m; i++) {
    int rn = A->row_number[i];
    A->value[i] = (T *)realloc(A->value[i], rn * sizeof(T));
    A->index[i] = (int *)realloc(A->index[i], (rn + 1) * sizeof(int));
    A->max_row_number[i] = rn;
  }

  return;
}

/*==================*
 *  Finish-Pattern  *
 *==================*/

template <typename T>
void finish_pattern(sparsecore<T> *A) noexcept {
  for (int i = 0; i < A->m; i++) {
    int rn = A->row_number[i];
    A->value[i] = (T *)calloc(rn, sizeof(T));
    A->index[i] = (int *)realloc(A->index[i], (rn + 1) * sizeof(int));
    A->max_row_number[i] = rn;
  }

  return;
}
}  // namespace Bembel

#endif
