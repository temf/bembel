// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_EvalSolution_
#define __BEMBEL_EvalSolution_

#include "Conversions.hpp"
#include "Spline.hpp"
#include "pot.h"

/**
 * @brief Eigen-wrapper for the evaluation of the solution.
 */
namespace Bembel {
namespace Sol {
// LaplaceSingleLayerCase

inline Eigen::Matrix<double, Eigen::Dynamic, 1> evalSolution(
    const Eigen::MatrixXd &gridpoints, Eigen::VectorXd &rho,
    Discretization<LaplaceSingle> &disc) {
  assert(gridpoints.cols() == 3 &&
         "Must be a Matrix with a point in each row!");
  const int gps = gridpoints.rows();
  vector3 *pts;
  pts = (vector3 *)malloc(sizeof(vector3) * gps);

  for (int i = 0; i < gps; i++) {
    pts[i] = vector3_make(gridpoints.row(i)(0), gridpoints.row(i)(1),
                          gridpoints.row(i)(2));
  }

  double *tmppot;

  pdeproblem *pde = disc.get_disc().pde;

  pot(rho.data(), &tmppot, pts, gps, &(disc.get_disc()), pde);

  Eigen::Matrix<double, Eigen::Dynamic, 1> potential(gps);

  for (int i = 0; i < gps; i++) {
    potential(i) = (tmppot[i]);
  }

  free(tmppot);
  free(pts);
  return potential;
}

// HelmholtzSingleLayerCase

inline Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> evalSolution(
    const Eigen::MatrixXd &gridpoints, Eigen::VectorXcd &rho,
    Discretization<HelmholtzSingle> &disc) {
  assert(gridpoints.cols() == 3 &&
         "Must be a Matrix with a point in each row!");
  vector3 *Q;
  double *Pot;
  double *ptrrho;

  const int nq = gridpoints.rows();
  const int rhosz = rho.rows();

  ptrrho = (double *)malloc(sizeof(double) * rhosz * 2);
  Q = (vector3 *)malloc(sizeof(vector3) * nq);

  for (int i = 0; i < rhosz; i++) {
    ptrrho[i] = rho(i).real();
    ptrrho[i + rhosz] = rho(i).imag();
  }

  for (int i = 0; i < nq; i++) {
    Q[i] = vector3_make(gridpoints.row(i)(0), gridpoints.row(i)(1),
                        gridpoints.row(i)(2));
  }

  potcomplex(ptrrho, &Pot, Q, nq, &(disc.get_disc()), &(disc.get_pde()._pde));

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> out(nq);

  for (int i = 0; i < nq; i++) {
    out(i) = (std::complex<double>(Pot[i], Pot[i + nq]));
  }

  free(Pot);
  free(Q);
  free(ptrrho);
  return out;
}

// MaxwellSingleLayerCase

inline Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 3> evalSolution(
    const Eigen::MatrixXd &gridpoints, Eigen::VectorXcd &rho,
    Discretization<MaxwellSingle> &disc) {
  assert(gridpoints.cols() == 3 &&
         "Must be a Matrix with a point in each row!");
  double *tmppot;

  vector3 *tmpgrid;
  const int gps = gridpoints.rows();
  tmpgrid = (vector3 *)malloc(sizeof(vector3) * gps);

  for (int i = 0; i < gps; i++) {
    tmpgrid[i] = vector3_make(gridpoints.row(i)(0), gridpoints.row(i)(1),
                              gridpoints.row(i)(2));
  }

  double *rhoptr = eigen2cmplxptr(rho);

  potMaxwell(rhoptr, &tmppot, tmpgrid, gps, &(disc.get_disc()),
             &(disc.get_pde()._pde));

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 3> potential(gps, 3);

  for (int i = 0; i < gps; i++) {
    potential.row(i) = Eigen::RowVector3cd(
        std::complex<double>(tmppot[i], tmppot[i + gps]),
        std::complex<double>(tmppot[i + 2 * gps], tmppot[i + 3 * gps]),
        std::complex<double>(tmppot[i + 4 * gps], tmppot[i + 5 * gps]));
  }

  free(rhoptr);
  free(tmpgrid);
  free(tmppot);

  return potential;
}
}  // namespace Sol
}  // namespace Bembel
#endif
