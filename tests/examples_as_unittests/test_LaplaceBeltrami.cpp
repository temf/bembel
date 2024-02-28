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

  int polynomial_degree = 1;
  int refinement_level = 0;

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
  std::cout << "Degree " << polynomial_degree << " Level " << refinement_level
            << "\t\t";
  // Build ansatz space
  AnsatzSpace<LaplaceBeltramiOperator> ansatz_space(geometry, refinement_level,
                                                    polynomial_degree);

  // Set up and compute discrete operator
  DiscreteLocalOperator<LaplaceBeltramiOperator> disc_op(ansatz_space);
  disc_op.compute();

  DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator> disc_cf(
      ansatz_space);
  disc_cf.get_linear_form().set_function(const_fun);
  disc_cf.compute();

  DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator> disc_lf(
      ansatz_space);
  disc_lf.get_linear_form().set_function(fun);
  disc_lf.compute();
  sw.tic();
  MatrixXd B = MatrixXd(disc_op.get_discrete_operator()) +
               disc_cf.get_discrete_linear_form() *
                   disc_cf.get_discrete_linear_form().transpose();
  LLT<MatrixXd> solver(B);
  VectorXd x = solver.solve(disc_lf.get_discrete_linear_form());
  double error = surfaceL2error(ansatz_space, x, refsol);
  std::cout << error << std::endl;

  std::cout << "\n" << std::string(60, '=') << "\n";
  return 0;
}
