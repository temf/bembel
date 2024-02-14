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

#include <Bembel/HomogenisedLaplace>
#include <Eigen/Dense>

#include <iostream>
#include <functional>

#include "tests/Test.hpp"

inline double k_per(Eigen::Vector3d in, Eigen::VectorXd coeffs,
    unsigned int deg);

inline Eigen::Vector3d Dk_per(Eigen::Vector3d in, Eigen::VectorXd coeffs,
    unsigned int deg);

Eigen::Matrix<double, 2, Eigen::Dynamic> tensorise(Eigen::ArrayXd xs);

int main() {
  using Eigen::VectorXd;
  using Eigen::Vector3d;
  using Eigen::ArrayXd;
  using Eigen::MatrixXd;

  double precision = 1e-7;
  unsigned int deg = Bembel::getDegree(precision);

  VectorXd cs = Bembel::getCoefficients(precision);

  unsigned int Npoints = 101;

  ArrayXd xs = ArrayXd::LinSpaced(Npoints, -0.5, 0.5);
  MatrixXd ys = tensorise(xs);

  double err = 0.0;

  Vector3d ex(1.0, 0.0, 0.0);
  Vector3d ey(0.0, 1.0, 0.0);
  Vector3d ez(0.0, 0.0, 1.0);

  std::function<double(Vector3d)> u = [cs, deg](Vector3d in) {
    return k_per(in, cs, deg);
  };
  std::function<double(Vector3d, unsigned int)> Du = [cs, deg](Vector3d in,
      unsigned int d) {
    return Dk_per(in, cs, deg)(d);
  };

  unsigned int k;
  Vector3d v;
  for (k = 0; k < Npoints * Npoints; k++) {
    v = Vector3d(-0.5, ys(0, k), ys(1, k));
    err += abs(k_per(v, cs, deg) - k_per(v + ex, cs, deg));
    err += abs(Du(v, 0) - Du(v + ex, 0));

    v = Vector3d(ys(0, k), -0.5, ys(1, k));
    err += abs(u(v) - u(v + ey));
    err += abs(Du(v, 1) - Du(v + ey, 1));

    v = Vector3d(ys(0, k), ys(1, k), -0.5);
    err += abs(u(v) - u(v + ez));
    err += abs(Du(v, 2) - Du(v + ez, 2));
  }

  err /= (6 * Npoints * Npoints);

  /* test the average pointwise error */
  assert(Bembel::Constants::isAlmostZero(err));

  /* test the gradient with the differential quotient */

  double h = 1e-8;
  Vector3d p(0.1, 0.2, -0.3);
  Vector3d dir = Vector3d::Random(3).normalized();
  Vector3d t = p + h * dir;
  double dq = (k_per(t, cs, deg) - k_per(p, cs, deg)) / h;
  double gr = Dk_per(p, cs, deg).dot(dir);

  /* test the error between the function and the difference quotient */
  assert(Bembel::Constants::isAlmostZero(dq - gr));

  return 0;
}

inline double k_per(Eigen::Vector3d in, Eigen::VectorXd coeffs,
    unsigned int deg) {

  return Bembel::k_mod(in)
      + Bembel::evaluate_solid_sphericals(in, coeffs, deg, false);
}

inline Eigen::Vector3d Dk_per(Eigen::Vector3d in, Eigen::VectorXd coeffs,
    unsigned int deg) {

  if (deg == 0) {
    return Bembel::Dk_mod(in);
  }

  Eigen::Vector3d res = Bembel::Dk_mod(in);
  res += (in / in.norm())
      * Bembel::evaluate_solid_sphericals(in, coeffs, deg, true);
  res += Bembel::functionalMatrix(in)
      * Bembel::evaluate_dsolid_sphericals(in, coeffs, deg);

  return res;
}

Eigen::Matrix<double, 2, Eigen::Dynamic> tensorise(Eigen::ArrayXd xs) {
  unsigned int len = xs.rows();
  Eigen::MatrixXd ys(2, len * len);

  unsigned int k;
  for (k = 0; k < len; k++) {
    ys.block(0, k * len, 1, len) = xs(k) * Eigen::MatrixXd::Ones(1, len);
    ys.block(1, k * len, 1, len) = xs.transpose();
  }

  return ys;
}
