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

#include "./Data.hpp"
#include "./Error.hpp"
#include "./Grids.hpp"

int main() {
  using Bembel::HomogenisedLaplaceSingleLayerOperator;
  using Bembel::HomogenisedLaplaceSingleLayerPotential;

  using Eigen::Vector3d;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;

  // Load geometry from file "cube_small.dat", which must be placed in the same
  // directory as the executable
  Bembel::Geometry geometry("cube_small.dat");

  // Define evaluation points for potential field, a tensor product grid of
  // 10*10*10 points in [0, 0.25]^3
  MatrixXd gridpoints = Bembel::Util::makeTensorProductGrid(
      VectorXd::LinSpaced(10, 0.05, 0.20), VectorXd::LinSpaced(10, 0.05, 0.20),
      VectorXd::LinSpaced(10, 0.05, 0.20));

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
  for (unsigned int polynomial_degree : { 0, 1, 2, 3 }) {
    // Iterate over refinement levels
    for (unsigned int refinement_level : { 0, 1, 2, 3 }) {
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
      Bembel::DiscreteOperator<MatrixXd, HomogenisedLaplaceSingleLayerOperator>
        disc_op(ansatz_space);
      disc_op.compute();

      // solve system
      Eigen::LLT<MatrixXd> llt;
      llt.compute(disc_op.get_discrete_operator());
      VectorXd rho = llt.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      Bembel::DiscretePotential<
          HomogenisedLaplaceSingleLayerPotential<
              HomogenisedLaplaceSingleLayerOperator>,
          HomogenisedLaplaceSingleLayerOperator> disc_pot(ansatz_space);
      disc_pot.set_cauchy_data(rho);
      VectorXd pot = disc_pot.evaluate(gridpoints);

      // print error
      std::cout << Bembel::maxPointwiseError<double>(pot, gridpoints, fun)
          << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "============================================================="
      "==========" << std::endl;

  return 0;
}

