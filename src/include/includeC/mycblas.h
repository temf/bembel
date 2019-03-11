// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_MYCBLAS__
#define __BEMBEL_C_MYCBLAS__

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
namespace Bembel {
/**
 *  \brief         Computes the euclidean norm of x.
 *
 *  \param[in]     n              Length of x.
 *  \param[in]     x              Vector of which the norm will be computed.
 *
 *  \return        Euclidean norm of x.
 *
 */
inline double mydnrm2(int n, double *x) {
  double norm2 = 0;

  for (int i = 0; i < n; ++i) norm2 += x[i] * x[i];

  return sqrt(norm2);
}

/**
 *  \brief         Scale x by alpha.
 *
 *  \param[in]     n              Length of x.
 *  \param[in]     alpha          Scale factor.
 *  \param[in,out] x              Vector which should be scaled.
 *
 */
inline int mydscal(int n, double alpha, double *x) {
  for (int i = 0; i < n; ++i) x[i] *= alpha;

  return 0;
}

/**
 *  \brief         Compute euclidean scalar product of x and y.
 *
 *  \param[in]     n              Length of x and y.
 *  \param[in]     x              First vector for scalar product.
 *  \param[in]     y              Second vector for scalar product.
 *
 *  \return        Euclidean scalar product of x and y.
 */
inline double myddot(int n, double *x, double *y) {
  double ddot = 0;

  for (int i = 0; i < n; ++i) ddot += x[i] * y[i];

  return ddot;
}

/**
 *  \brief         Compute vector y plus x.
 *
 *  \param[in]     n              Length of x and y.
 *  \param[in]     x              Vector to add.
 *  \param[in,out] y              Vector to which x will be added.
 */
inline int mydxpy(int n, double *x, double *y) {
  for (int i = 0; i < n; ++i) y[i] += x[i];

  return 0;
}

/**
 *  \brief         Compute vector y plus alpha*x.
 *
 *  \param[in]     n              Length of x and y.
 *  \param[in]     alpha          Scale factor of x.
 *  \param[in]     x              Vector to add.
 *  \param[in,out] y              Vector to which alpha*x will be added.
 */
inline int mydaxpy(int n, double alpha, double *x, double *y) {
  for (int i = 0; i < n; ++i) y[i] += alpha * x[i];

  return 0;
}

/**
 *  \brief         Compute vector y plus matrix*vector A*x.
 *
 *  \param[in]     trans          Take A transposed or not, N or T
 *  \param[in]     n              Length of x and y.
 *  \param[in]     A              Matrix in row major format.
 *  \param[in]     x              Vector x.
 *  \param[in,out] y              Vector to which A*x will be added.
 */
inline int myqdgemv(char trans, int n, double *A, double *x, double *y) {
  if (trans == 'N') {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) y[i] += A[i * n + j] * x[j];
  } else if (trans == 'T') {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) y[j] += A[i * n + j] * x[i];
  }

  return 0;
}

/**
 *  \brief Header for LAPACK function, Documentation under
 *         http://www.netlib.no/netlib/lapack/double/dgetrf.f
 */
extern void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);
/**
 *  \brief Header for LAPACK function, Documentation under
 *         http://www.netlib.no/netlib/lapack/double/dgetri.f
 */
extern void dgetri_(int *n, double *A, int *lda, int *ipiv, double *work,
                    int *lwork, int *info);

inline int mydinv(int n, double *A) {
  // int NB = 3200;
  // int *ipiv;
  // double *work;
  // int lwork;
  // int info;
  //
  // /* prepare lapack parameters */
  // ipiv = (int*)calloc (1, n * sizeof(int));
  // work = (double*)calloc (1, n * NB * sizeof(double));
  // lwork = n * NB;
  //
  // /* invert matrix with LAPACK */
  // dgetrf_ (&n, &n, A, &n, ipiv, &info);
  //
  // assert(!info || "dgetrf of full matrix failed");
  //
  // dgetri_ (&n, A, &n, ipiv, work, &lwork, &info);
  //
  // assert(!info || "dgetri of full matrix failed\n");
  //
  // free(ipiv);
  // free(work);

  Eigen::Map<Eigen::MatrixXd> mat(A, n, n);
  mat = mat.inverse();

  return 0;
}
}  // namespace Bembel
#endif
