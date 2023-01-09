/*
 * example_HomogenisedLaplaceSingleLayerFull.cpp
 *
 *  Created on: 9 Jan 2023
 *      Author: ricrem00
 */

#include <iostream>

#include <Eigen/Dense>

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/HomogenisedLaplace>
#include <Bembel/LinearForm>

#include "Data.hpp"
#include "Error.hpp"
#include "Grids.hpp"

int main() {
	using namespace Bembel;
	using namespace Eigen;

	// Load geometry from file "sphere.dat", which must be placed in the same
	// directory as the executable
	Geometry geometry("../geo/cube_small.dat");

	// Define evaluation points for potential field, a tensor product grid of
	// 7*7*7 points in [0, 0.25]^3
	MatrixXd gridpoints = Util::makeTensorProductGrid(
			VectorXd::LinSpaced(10, 0.05, 0.20), VectorXd::LinSpaced(10, 0.05, 0.20),
			VectorXd::LinSpaced(10, 0.05, 0.20));

	// Define analytical solution using lambda function, in this case a harmonic
	// function, see Data.hpp
	std::function<double(Vector3d)> fun = [](Vector3d in) {
		return Data::HarmonicFunction(in);
	};

	// Define the precision of the fundamental solution of the potential and the operator
	double precision = 1e-7;
	HomogenisedLaplaceSingleLayerOperator::setPrecision(precision);

	std::cout << "\n============================================================="
			"==========\n";
	// Iterate over polynomial degree.
	for (auto polynomial_degree : {0, 1, 2, 3}) {
		// Iterate over refinement levels
		for (auto refinement_level : {0, 1, 2, 3}) {
			std::cout << "Degree " << polynomial_degree << " Level "
					<< refinement_level << "\t\t";
			// Build ansatz space
			AnsatzSpace<HomogenisedLaplaceSingleLayerOperator> ansatz_space(
					geometry, refinement_level, polynomial_degree);

			// Set up load vector
			DiscreteLinearForm<DirichletTrace<double>, HomogenisedLaplaceSingleLayerOperator>
			disc_lf(ansatz_space);
			disc_lf.get_linear_form().set_function(fun);
			disc_lf.compute();

			// Set up and compute discrete operator
			DiscreteOperator<MatrixXd, HomogenisedLaplaceSingleLayerOperator> disc_op(
					ansatz_space);
			disc_op.compute();

			// solve system
			LLT<MatrixXd> llt;
			llt.compute(disc_op.get_discrete_operator());
			auto rho = llt.solve(disc_lf.get_discrete_linear_form());

			// evaluate potential
			DiscretePotential<HomogenisedLaplaceSingleLayerPotential<HomogenisedLaplaceSingleLayerOperator>,
			HomogenisedLaplaceSingleLayerOperator>
			disc_pot(ansatz_space);
			disc_pot.set_cauchy_data(rho);
			auto pot = disc_pot.evaluate(gridpoints);

			// print error
			std::cout << maxPointwiseError<double>(pot, gridpoints, fun) << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << "============================================================="
			"=========="
			<< std::endl;

	return 0;
}

