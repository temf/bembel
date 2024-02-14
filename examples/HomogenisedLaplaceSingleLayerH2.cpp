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

#include <Eigen/Dense>

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/HomogenisedLaplace>
#include <Bembel/LinearForm>

#include <iostream>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

int main() {
  using Bembel::HomogenisedLaplaceSingleLayerOperator;
  using Bembel::HomogenisedLaplaceSingleLayerPotential;

  using Eigen::Vector3d;
  using Eigen::VectorXd;
  using Eigen::H2Matrix;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  // Load geometry from file "cube_small.dat", which must be placed in the same
  // directory as the executable
  Bembel::Geometry geometry("cube_small.dat");

  // Define evaluation points for potential field, a tensor product grid of
  // 10*10*10 points in [0, 0.25]^3
  Eigen::MatrixXd gridpoints = Bembel::Util::makeTensorProductGrid(
      VectorXd::LinSpaced(10, 0.11, 0.14), VectorXd::LinSpaced(10, 0.11, 0.14),
      VectorXd::LinSpaced(10, 0.11, 0.14));

  // Define analytical solution using lambda function, in this case a harmonic
  // function, see Data.hpp
  std::function<double(Vector3d)> fun = [](Vector3d in) {
    return Bembel::Data::HarmonicFunction(in);
  };

  // Define the precision of the fundamental solution of the
  // potential and the operator
  double precision = 1e-7;
  HomogenisedLaplaceSingleLayerOperator::setPrecision(precision);

  std::cout << "\n============================================================="
      "==========\n";
  // Iterate over polynomial degree.
  for (int polynomial_degree = 0; polynomial_degree < polynomial_degree_max + 1;
      polynomial_degree++) {
    VectorXd error(refinement_level_max + 1);
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
        refinement_level++) {
      std::cout << "Degree " << polynomial_degree << " Level "
          << refinement_level << "\t\t";
      // Build ansatz space
      Bembel::AnsatzSpace<HomogenisedLaplaceSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      Bembel::DiscreteLinearForm<Bembel::DirichletTrace<double>,
      HomogenisedLaplaceSingleLayerOperator> disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up and compute discrete operator
      Bembel::DiscreteOperator<H2Matrix<double>,
      HomogenisedLaplaceSingleLayerOperator> disc_op(ansatz_space);
      disc_op.compute();

      // solve system
      Eigen::ConjugateGradient<H2Matrix<double>, Eigen::Lower | Eigen::Upper,
      Eigen::IdentityPreconditioner> cg;
      cg.compute(disc_op.get_discrete_operator());
      VectorXd rho = cg.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      Bembel::DiscretePotential<
      HomogenisedLaplaceSingleLayerPotential<
      HomogenisedLaplaceSingleLayerOperator>,
      HomogenisedLaplaceSingleLayerOperator> disc_pot(ansatz_space);
      disc_pot.set_cauchy_data(rho);
      VectorXd pot = disc_pot.evaluate(gridpoints);

      // print error
      error(refinement_level) = Bembel::maxPointwiseError<double>(pot,
          gridpoints, fun);
      std::cout << error(refinement_level) << std::endl;
    }

    // estimate the rate of convergence and check whether it is at least
    // 90% of the expected value
    assert(Bembel::checkRateOfConvergence(error.tail(2),
        2*polynomial_degree + 3, 0.9));

    std::cout << std::endl;
  }
  std::cout << "============================================================="
      "==========" << std::endl;

  return 0;
}
