// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "Discretization.hpp"
#include "Geometry.hpp"
#include "HierarchicalMatrix.hpp"
#include "PDEproblem.hpp"

#include "Data.hpp"
#include "Error.hpp"
#include "EvalSolution.hpp"
#include "Grids.hpp"
#include "Logger.hpp"
#include "Rhs.hpp"
#include "Stopwatch.hpp"

/**
 *  @brief         This is a demonstration how to use the high-level classes &
 * templates of bembel to run a Laplace-BEM based on the single layer operator.
 *
 */
int main() {
  using namespace Bembel;

  // Set some parameters, see documentation
  const int multipoleDegree = 16;
  const int knotRepetition = 1;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry myGeom("sphere.dat");

  // Define single layer operator
  LaplaceSingle myLap;

  // Define evaluation points for potential field, a tensor product grid of
  // 7*7*7 points in [-.1,.1]^3
  Eigen::MatrixXd gridpoints =
      Util::makeTensorProductGrid(Eigen::VectorXd::LinSpaced(7, -.1, .1),
                                  Eigen::VectorXd::LinSpaced(7, -.1, .1),
                                  Eigen::VectorXd::LinSpaced(7, -.1, .1));

  // Define analytical solution using lambda function, in this case a harmonic
  // function, see Data.hpp
  std::function<double(Eigen::Vector3d)> fun = [](Eigen::Vector3d in) {
    return Data::HarmonicFunction(in);
  };

  std::cout << "\n============================================================="
               "==========\n";
  // Iterate over polynomial degree.
  for (auto P : {0, 1, 2}) {
    std::cout << std::endl;
    // Initialize timer and log-file
    Util::Stopwatch sw;
    Util::Logger<7> log("test/LaplaceSingle_" + std::to_string(P) + ".log");
    log.both("P", "M", "t_disc", "t_rhs", "t_mat", "t_solve", "t_tot", "error");
    // Iterate over refinement levels
    for (auto M : {0, 1, 2, 3}) {
      sw.start();

      // Build discretization
      Discretization<LaplaceSingle> myDisc(myGeom, myLap, P, knotRepetition, M);
      const double timeDisc = sw.lap();

      // Build right-hand side
      Eigen::VectorXd rhs = Rhs::computeRhs(myDisc, fun);
      const double timeRhs = sw.lap();

      // Assemble system matrix
      Eigen::HierarchicalMatrix<LaplaceSingle> myH(myDisc, multipoleDegree);
      const double timeMat = sw.lap();

      // Solve myH*rho=rhs using CG from Eigen
      Eigen::ConjugateGradient<Eigen::HierarchicalMatrix<LaplaceSingle>,
                               Eigen::Lower | Eigen::Upper,
                               Eigen::IdentityPreconditioner>
          cg;
      cg.compute(myH);
      Eigen::VectorXd rho = cg.solve(rhs);
      const double timeSolve = sw.lap();

      // Evaluate potential
      Eigen::MatrixXd pot = Sol::evalSolution(gridpoints, rho, myDisc);
      const double timeTot = sw.stop();

      // Write data to log-file
      log.both(P, M, timeDisc, timeRhs, timeMat, timeSolve, timeTot,
               maxPointwiseError(pot, gridpoints, fun));
    }
  }
  std::cout << "============================================================="
               "=========="
            << std::endl;
  return (0);
}
