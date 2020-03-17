// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include "Test.hpp"
#include <Bembel/Spline>


int main() {
  using namespace Bembel;

  // We test the Bernstein Basis of the BasisHandler against the deBoor code
  for (int p = 0; p < 18; ++p) {
    for (auto x : Test::Constants::eq_points) {
      Eigen::VectorXd result1 = Eigen::VectorXd::Zero(p + 1);
      Basis::BasisHandler<double>::phi(p, &result1, 1, x);

      Eigen::VectorXd result2 = Eigen::VectorXd::Zero(p + 1);
      Eigen::MatrixXd coef = Eigen::VectorXd::Zero(p + 1).transpose();
      for (int i = 0; i < p + 1; ++i) {
        coef(i) = 1;
        std::vector<double> v = {x};
        result2(i) = Spl::DeBoor(coef, Spl::MakeBezierKnotVector(p + 1), v)(0);
        coef(i) = 0;
      }
      BEMBEL_TEST_IF((result1 - result2).norm() < Test::Constants::coefficient_accuracy);
    }
  }

  // Now, we do the same for the derivatives
  for (int p = 1; p < 18; ++p) {
    for (auto x : Test::Constants::eq_points) {
      Eigen::VectorXd result1 = Eigen::VectorXd::Zero(p + 1);
      Basis::BasisHandler<double>::phiDx(p, &result1, 1, x);

      Eigen::VectorXd result2 = Eigen::VectorXd::Zero(p + 1);
      Eigen::MatrixXd coef = Eigen::VectorXd::Zero(p + 1).transpose();
      for (int i = 0; i < p + 1; ++i) {
        coef(i) = 1;
        std::vector<double> v = {x};
        result2(i) =
            Spl::DeBoorDer(coef, Spl::MakeBezierKnotVector(p + 1), v)(0);
        coef(i) = 0;
      }
      BEMBEL_TEST_IF((result1 - result2).norm() < Test::Constants::coefficient_accuracy);
    }
  }

  return 0;
}
