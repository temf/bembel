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

#include <Bembel/IO>
#include <Bembel/LaplaceBeltrami>
#include <Bembel/LinearForm>
#include <Eigen/Dense>
#include <iostream>

#include "Bembel/src/util/surfaceL2error.hpp"
#include "examples/Error.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;
  IO::Stopwatch sw;

  int polynomial_degree_max = 4;
  int refinement_level_max = 5;

  std::function<double(const Vector3d &)> fun = [](const Vector3d &in) {
    // return 0.5 / sqrt(M_PI);
    return sqrt(3 / M_PI) * in(2);
  };
  std::function<double(const Vector3d &)> refsol = [](const Vector3d &in) {
    // return 0.5 / sqrt(M_PI);
    return 0.5 * sqrt(3 / M_PI) * in(2);
  };

  std::function<double(const Vector3d &)> const_fun = [](const Vector3d &in) {
    return 1.0;
  };
  Geometry geometry("sphere.dat");
  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  for (int polynomial_degree = 1; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level << "\t\t";
      // Build ansatz space
      AnsatzSpace<LaplaceBeltramiOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up and compute discrete operator
      DiscreteLocalOperator<LaplaceBeltramiOperator> disc_op(ansatz_space);
      disc_op.compute();

      DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator>
          disc_cf(ansatz_space);
      disc_cf.get_linear_form().set_function(const_fun);
      disc_cf.compute();

      DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();
      sw.tic();
      MatrixXd B = MatrixXd(disc_op.get_discrete_operator()) +
                          disc_cf.get_discrete_linear_form() *
                              disc_cf.get_discrete_linear_form().transpose();
      LLT<MatrixXd> solver(B);
      VectorXd x = solver.solve(disc_lf.get_discrete_linear_form());
      error(refinement_level) = surfaceL2error(ansatz_space, x, refsol);
      std::cout << error(refinement_level) << std::endl;
    }

    // estimate rate of convergence and check whether it is at least 90% of the
    // expected value
    assert(
        checkRateOfConvergence(error.tail(2), polynomial_degree + 1, 0.9));

    std::cout << std::endl;
    }
  std::cout << "\n" << std::string(60, '=') << "\n";
  return 0;
}
