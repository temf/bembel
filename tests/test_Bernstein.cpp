// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/Spline>

#include "tests/Test.hpp"

int main() {
  using namespace Bembel;

  // We test the Bernstein with given control points against the deBoor code
  constexpr int P = Bembel::Constants::MaxP;
  double* coefs = new double[P + 1];

  for (int i = 0; i < P + 1; ++i) {
    coefs[i] = (i + 1) / static_cast<double>(P + 1);
  }

  Eigen::Map<Eigen::MatrixXd> coefs_vector(coefs, 1, P + 1);

  for (int p = 0; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      double result1 = Basis::ShapeFunctionHandler::evalCoef(p, coefs, x);

      std::vector<double> v = {x};
      double result2 =
          Spl::DeBoor(Eigen::MatrixXd(coefs_vector.leftCols(p + 1)),
                      Spl::MakeBezierKnotVector(p + 1), v)(0);

      BEMBEL_TEST_IF(std::abs(result1 - result2) <
                     Test::Constants::coefficient_accuracy);
    }
  }

  // Now, we do the same for the derivatives
  for (int p = 1; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      double result1 = Basis::ShapeFunctionHandler::evalDerCoef(p, coefs, x);

      std::vector<double> v = {x};
      double result2 =
          Spl::DeBoorDer(Eigen::MatrixXd(coefs_vector.leftCols(p + 1)),
                         Spl::MakeBezierKnotVector(p + 1), v)(0);

      BEMBEL_TEST_IF(std::abs(result1 - result2) <
                     Test::Constants::coefficient_accuracy);
    }
  }

  delete[] coefs;
  return 0;
}
