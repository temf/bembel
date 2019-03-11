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
 * templates of bembel to run a Maxwell-BEM based on the EFIE.
 *
 */
int main() {
  using namespace Bembel;

  // Set some parameters, see documentation
  const int multipoleDegree = 16;
  const int knotRepetition = 1;
  const std::complex<double> wavenumber(1, 0);

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry myGeom("sphere.dat");

  // Define single layer operator
  MaxwellSingle myMax(wavenumber);

  // Define evaluation points for scattered field, sphere of radius 2, 10*10
  // points.
  Eigen::MatrixXd gridpoints = Util::makeSphereGrid(2, 10);

  // Define analytical solution using lambda function, in this case a dipole,
  // see Data.hpp
  const std::function<Eigen::Vector3cd(Eigen::Vector3d, std::complex<double>)>
      fun = [](Eigen::Vector3d pt, std::complex<double> kappa) {
        const Eigen::Vector3d position(0.2, 0.2, 0.2);
        const Eigen::Vector3d length(0, 0.1, 0.1);
        return Data::Dipole(pt, kappa, position, length);
      };

  std::cout << "\n============================================================="
               "==========\n";
  // Iterate over polynomial degree.
  for (auto P : {1, 2}) {
    std::cout << std::endl;
    // Initialize timer and log-file
    Util::Stopwatch sw;
    Util::Logger<7> log("test/MaxwellSingle_" + std::to_string(P) + ".log");
    log.both("P", "M", "t_disc", "t_rhs", "t_mat", "t_solve", "t_tot", "error");
    // Iterate over refinement levels
    for (auto M : {0, 1, 2}) {
      sw.start();

      // Build discretization
      Discretization<MaxwellSingle> myDisc(myGeom, myMax, P, knotRepetition, M);
      const double timeDisc = sw.lap();

      // Build right-hand side
      Eigen::VectorXcd rhs = Rhs::computeRhs(myDisc, fun);
      const double timeRhs = sw.lap();

      // Assemble system matrix
      Eigen::HierarchicalMatrix<MaxwellSingle> myH(myDisc, multipoleDegree);
      const double timeMat = sw.lap();

      // Solve myH*rho=rhs using GMRES from Eigen
      Eigen::GMRES<Eigen::HierarchicalMatrix<MaxwellSingle>,
                   Eigen::IdentityPreconditioner>
          gmres;
      gmres.setTolerance(1e-20);
      gmres.set_restart(1000);
      gmres.compute(myH);
      Eigen::VectorXcd rho = gmres.solve(rhs);
      const double timeSolve = sw.lap();

      // Evaluate scattered field
      Eigen::MatrixXcd pot = Sol::evalSolution(gridpoints, rho, myDisc);
      const double timeTot = sw.stop();

      // Write data to log-file
      log.both(P, M, timeDisc, timeRhs, timeMat, timeSolve, timeTot,
               maxPointwiseError(pot, gridpoints, fun, myMax.get_wavenumber()));
    }
  }
  std::cout << "============================================================="
               "=========="
            << std::endl;
  return (0);
}
