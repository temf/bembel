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
  using namespace Eigen;

  // We test the Bernstein with given control points against the deBoor code
  constexpr int P = Bembel::Constants::MaxP;
  VectorXd buffer(P + 1);
  VectorXd coefs =
      VectorXd::LinSpaced(P + 1, 1, P + 1) / static_cast<double>(P + 1);

  for (int p = 0; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      buffer.setZero();
      Basis::ShapeFunctionHandler::evalBasis(p, buffer.data(), x);
      double result1 = buffer.dot(coefs);

      std::vector<double> v = {x};
      double result2 = Spl::DeBoor(MatrixXd(coefs.head(p + 1).transpose()),
                                   Spl::MakeBezierKnotVector(p + 1), v)(0);

      BEMBEL_TEST_IF(std::abs(result1 - result2) <
                     Test::Constants::coefficient_accuracy);
    }
  }

  // Now, we do the same for the derivatives
  for (int p = 1; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      buffer.setZero();
      Basis::ShapeFunctionHandler::evalDerBasis(p, buffer.data(), x);
      double result1 = buffer.dot(coefs);

      std::vector<double> v = {x};
      double result2 = Spl::DeBoorDer(MatrixXd(coefs.head(p + 1).transpose()),
                                      Spl::MakeBezierKnotVector(p + 1), v)(0);

      BEMBEL_TEST_IF(std::abs(result1 - result2) <
                     Test::Constants::coefficient_accuracy);
    }
  }

  return 0;
}
