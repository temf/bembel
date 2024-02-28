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

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/LinearForm>

#include "tests/TestGeometries.hpp"

class TestOperatorDivC;
template <>
struct Bembel::LinearOperatorTraits<TestOperatorDivC> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum { OperatorOrder = 0, Form = DifferentialForm::DivConforming };
};

int main() {
  Test::TestGeometryWriter::writeScreen();
  Bembel::Geometry geometry("test_Screen.dat");

  const int refinement_level = 0;
  const int polynomial_degree = 1;
  Bembel::AnsatzSpace<TestOperatorDivC> ansatz_space(geometry, refinement_level,
                                                     polynomial_degree);

  std::function<Eigen::Vector3d(Eigen::Vector3d)> fun = [](Eigen::Vector3d in) {
    return Eigen::Vector3d(in(0), in(1), in(2));
  };

  Bembel::DiscreteLinearForm<Bembel::TangentialTrace<double>, TestOperatorDivC>
      disc_lf(ansatz_space);
  disc_lf.get_linear_form().set_function(fun);
  disc_lf.compute();

  Eigen::VectorXd ref_sol(4);
  ref_sol << 0.25, 0.25, -0.25, -0.25;

  BEMBEL_TEST_IF((disc_lf.get_discrete_linear_form() - ref_sol).norm() <
                 Test::Constants::coefficient_accuracy);

  return 0;
}
