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
#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/Laplace>
#include <Bembel/LinearForm>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("torus.dat");

  // Define 100 evaluation points for potential field, on a circle in xy-plane
  // in the middle of the torus
  int n_gridpoints = 100;
  double h = 2. * BEMBEL_PI / n_gridpoints;
  Eigen::Matrix<double, Eigen::Dynamic, 3> gridpoints = MatrixXd::Zero(100, 3);
  for (int i = 0; i < n_gridpoints; ++i) {
    gridpoints(i, 0) = 2. * cos(h * i);
    gridpoints(i, 1) = 2. * sin(h * i);
  }

  // Define analytical solution using lambda function, in this case a
  // harmonic function, see Data.hpp
  std::function<double(Vector3d)> fun = [](Vector3d in) {
    return Data::HarmonicFunction(in);
  };

  // Iterate over polynomial degree.
  for (int polynomial_degree = 0; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);
    IO::Logger<12> logger("log_LaplaceSingle_" +
                          std::to_string(polynomial_degree) + ".log");
    logger.both("P", "M", "error");
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level << "\t\t";
      // Build ansatz space
      AnsatzSpace<LaplaceSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      DiscreteLinearForm<DirichletTrace<double>, LaplaceSingleLayerOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up and compute discrete operator
      DiscreteOperator<H2Matrix<double>, LaplaceSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.compute();

      // solve system
      ConjugateGradient<H2Matrix<double>, Lower | Upper, IdentityPreconditioner>
          cg;
      cg.compute(disc_op.get_discrete_operator());
      auto rho = cg.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<LaplaceSingleLayerPotential<LaplaceSingleLayerOperator>,
                        LaplaceSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // compute reference and compute error
      VectorXd pot_ref(gridpoints.rows());
      for (int i = 0; i < gridpoints.rows(); ++i)
        pot_ref(i) = fun(gridpoints.row(i));
      error(refinement_level) = (pot - pot_ref).cwiseAbs().maxCoeff();

      logger.both(polynomial_degree, refinement_level, error(refinement_level));

      // we only need one visualization per polynomial degree
      if (refinement_level == 3) {
        // export geometry with density
        VTKSurfaceExport writer_geo(geometry, 5);
        FunctionEvaluator<LaplaceSingleLayerOperator> evaluator(ansatz_space);
        evaluator.set_function(rho);
        std::function<double(int, const Eigen::Vector2d &)> density =
            [&](int patch_number,
                const Eigen::Vector2d &reference_domain_point) {
              return evaluator.evaluateOnPatch(patch_number,
                                               reference_domain_point)(0);
            };
        writer_geo.addDataSet("Density", density);
        writer_geo.writeToFile("LaplaceSingleDensity.vtp");

        // export point evaluations, can be visualized using glyphs in paraview
        VTKPointExport writer_points(gridpoints);
        writer_points.addDataSet("Potential", pot);
        writer_points.writeToFile("LaplaceSinglePoints.vtp");
      }
    }

    std::cout << std::endl;
  }
  std::cout << "============================================================="
               "=========="
            << std::endl;

  return 0;
}
