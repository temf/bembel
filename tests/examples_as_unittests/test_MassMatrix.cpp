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

#include <Bembel/Geometry>
#include <Bembel/IO>
#include <Bembel/Identity>
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
    return in(0);
  };
  Geometry geometry("sphere.dat");
  std::cout << "\n" << std::string(60, '=') << "\n";
  std::cout << "Degree " << polynomial_degree << " Level " << refinement_level
            << "\t\t";
  // Build ansatz space
  AnsatzSpace<MassMatrixScalarCont> ansatz_space(geometry, refinement_level,
                                                 polynomial_degree);

  // Set up and compute discrete operator
  sw.tic();
  DiscreteLocalOperator<MassMatrixScalarCont> disc_op(ansatz_space);
  disc_op.compute();
  auto difft = sw.toc();
  DiscreteLinearForm<DirichletTrace<double>, MassMatrixScalarCont> disc_lf(
      ansatz_space);
  disc_lf.get_linear_form().set_function(fun);
  disc_lf.compute();
  sw.tic();
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  // Compute the ordering permutation vector from the structural pattern of
  solver.analyzePattern(disc_op.get_discrete_operator());
  // Compute the numerical factorization
  solver.factorize(disc_op.get_discrete_operator());
  // Use the factors to solve the linear system
  auto x = solver.solve(disc_lf.get_discrete_linear_form());
  difft = sw.toc();
  double error = surfaceL2error(ansatz_space, x, fun);
  std::cout << error << std::endl;

  // The VTKwriter sets up initial geomety information.
  std::cout << "\n" << std::string(60, '=') << "\n";

  return 0;
}
