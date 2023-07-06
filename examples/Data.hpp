// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef EXAMPLES_DATA_HPP_
#define EXAMPLES_DATA_HPP_

#include <Eigen/Dense>
#include <cmath>

namespace Bembel {
namespace Data {

/*	@brief This function implements a harmonic function, which, in the
 * interior dormain, satisfies the Laplace equation.
 *
 */
inline double HarmonicFunction(Eigen::Vector3d in) {
  return (4 * in(0) * in(0) - 3 * in(1) * in(1) - in(2) * in(2));
}

/*	@brief This function implements the Helmholtz fundamental solution,
 * which, if center is placed in the interior domain, satisfies the Helmholtz
 * equation and radiation condition in the exterior domain.
 *
 * Note that the sign in the fundamental solution is linked to the sign of kappa
 * in this function.
 */
inline std::complex<double> HelmholtzFundamentalSolution(
    Eigen::Vector3d pt, std::complex<double> kappa,
    Eigen::Vector3d center = Eigen::Vector3d(0, 0, 0)) {
  return std::exp(-std::complex<double>(0, 1) * kappa * (pt - center).norm()) /
         (pt - center).norm();
}

/* @brief This function implement a Hertz Dipole as in page 411 of J.D.Jacksons
 * "Classical Electrodynamics", 3rd ed., which, if the dipole axis given by
 * position and length remains in the interior domain, satisfies the curl-curl
 * equation and the radiation condition in the exterior domain.
 *
 * Note that the sign in the fundamental solution is linked to the sign of kappa
 * in this function.
 */
inline Eigen::Vector3cd Dipole(Eigen::Vector3d pt, std::complex<double> kappa,
                               Eigen::Vector3d position,
                               Eigen::Vector3d length) {
  constexpr std::complex<double> i(0, 1);
  const Eigen::Vector3d c = pt - position;
  double r = c.norm();
  const Eigen::Vector3d n = c / r;
  const std::complex<double> expc = std::exp(-i * kappa * r);
  const Eigen::Vector3cd E =
      (kappa * kappa * expc / r) * (n.cross(length).cross(n));
  const Eigen::Vector3d h = 3 * n.dot(length) * n - length;
  const std::complex<double> ec =
      ((1 / (r * r * r)) - (i * (-kappa) / (r * r))) * expc;
  return E + ec * h;
}

}  // namespace Data
}  // namespace Bembel
#endif  // EXAMPLES_DATA_HPP_
