
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/Laplace>
#include <Bembel/LinearForm>

#include "Data.hpp"
#include "Error.hpp"
#include "Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("torus.dat");

  // Define evaluation points for potential field, a tensor product grid of
  // 7*7*7 points in [-.1,.1]^3
  MatrixXd gridpoints = Util::makeTensorProductGrid(
      VectorXd::LinSpaced(10, -.1, .1), VectorXd::LinSpaced(10, -2.1, -1.9),
      VectorXd::LinSpaced(10, -.1, .1));

  // Define analytical solution using lambda function, in this case a harmonic
  // function, see Data.hpp
  std::function<double(Vector3d)> fun = [](Vector3d in) {
    return Data::HarmonicFunction(in);
  };

  // Iterate over polynomial degree.
  for (auto polynomial_degree : {0, 1, 2}) {
    // Iterate over refinement levels
    IO::Logger<12> logger("log_LaplaceSingle_" +
                          std::to_string(polynomial_degree) + ".log");
    logger.both("P", "M", "error");

    for (auto refinement_level : {0, 1, 2, 3, 4}) {
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

      auto error = maxPointwiseError<double>(pot, gridpoints, fun);

      logger.both(polynomial_degree, refinement_level, error);

      // we only need one visualization
      if (refinement_level == 3 && polynomial_degree == 2) {
        VTKSurfaceExport writer(geometry, 5);

        FunctionEvaluator<LaplaceSingleLayerOperator> evaluator(ansatz_space);
        evaluator.set_function(rho);

        std::function<double(int, const Eigen::Vector2d &)> density =
            [&](int patch_number,
                const Eigen::Vector2d &reference_domain_point) {
              return evaluator.evaluateOnPatch(patch_number,
                                               reference_domain_point)(0);
            };
        writer.addDataSet("Density", density);
        writer.writeToFile("LaplaceSingle.vtp");
      }
    }
  }

  return 0;
}
