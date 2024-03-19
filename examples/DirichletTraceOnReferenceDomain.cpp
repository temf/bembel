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

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  std::function<double(const Vector3d &)> fun = [](const Vector3d &in) {
    return in(0);
  };
  Geometry geometry("sphere.dat");
  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  for (int polynomial_degree = 0; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level << std::endl;
      // Build ansatz space
      AnsatzSpace<MassMatrixScalarDisc> ansatz_space(geometry, refinement_level,
                                                     polynomial_degree);

      // Set up and compute discrete operator
      DiscreteLocalOperator<MassMatrixScalarDisc> disc_op(ansatz_space);
      disc_op.compute();
      DiscreteLinearForm<DirichletTrace<double>, MassMatrixScalarDisc> disc_lf(
          ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up solver
      SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
      solver.analyzePattern(disc_op.get_discrete_operator());
      solver.factorize(disc_op.get_discrete_operator());
      VectorXd x = solver.solve(disc_lf.get_discrete_linear_form());

      // Set up function on reference domain
      FunctionEvaluator<MassMatrixScalarDisc> evaluator(ansatz_space);
      evaluator.set_function(x);
      std::function<double(int, const Vector2d &)> fun_disc =
          [&](int p, const Vector2d &x) {
            return evaluator.evaluateOnPatch(p, x)(0);
          };

      // Set up linear form on reference domain
      DiscreteLinearForm<DirichletTraceOnReferenceDomain<double>,
                         MassMatrixScalarDisc>
          disc_lf2(ansatz_space);
      disc_lf2.get_linear_form().set_function(fun_disc);
      disc_lf2.compute();

      // check correctness
      assert(std::abs((disc_lf.get_discrete_linear_form() -
                       disc_lf2.get_discrete_linear_form())
                          .norm() /
                      disc_lf.get_discrete_linear_form().norm()) < 1e-10);
    }

    std::cout << std::endl;
  }
  // The VTKwriter sets up initial geomety information.
  std::cout << "\n" << std::string(60, '=') << "\n";

  return 0;
}
