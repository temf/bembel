// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>
#include "Discretization.hpp"
#include "Geometry.hpp"
#include "HierarchicalMatrix.hpp"
#include "PDEproblem.hpp"

#include "Data.hpp"
#include "Grids.hpp"
#include "Logger.hpp"
#include "Rhs.hpp"
#include "Stopwatch.hpp"
#include "Visualize.hpp"

/**
 *  @brief         This is a demonstration on how to use the visualization
 * routine. Apart from the commented lines for visualization, this example
 * follows the lines of example_HelmholtzSingle.cpp
 *
 */
int main() {
  using namespace Bembel;

  const int multipoleDegree = 16;
  const int knotRepetition = 1;
  const std::complex<double> wavenumber = std::complex<double>(1, 0);

  Geometry myGeom("sphere.dat");

  HelmholtzSingle myHelm(wavenumber);

  const std::function<std::complex<double>(Eigen::Vector3d,
                                           std::complex<double>)>
      fun = [](Eigen::Vector3d pt, std::complex<double> kappa) {
        return Data::HelmholtzFundamentalSolution(pt, kappa);
      };

  for (auto P : {2}) {
    for (auto M : {2}) {
      Discretization<HelmholtzSingle> myDisc(myGeom, myHelm, P, knotRepetition,
                                             M);
      Eigen::VectorXcd rhs = Rhs::computeRhs(myDisc, fun);
      Eigen::HierarchicalMatrix<HelmholtzSingle> myH(myDisc, multipoleDegree);
      Eigen::GMRES<Eigen::HierarchicalMatrix<HelmholtzSingle>,
                   Eigen::IdentityPreconditioner>
          gmres;
      gmres.setTolerance(1e-20);
      gmres.set_restart(1000);
      gmres.compute(myH);
      Eigen::VectorXcd rho = gmres.solve(rhs);

      // This just print the geometry, with the default name "geometry.vtk"
      Bembel::Vis::plotDiscretizationToVTK(myDisc);

      // This visualizes rho on the surface as well. It is stored under the
      // default name "<pde>.vtk", in this case "helmholtz.vtk"
      Bembel::Vis::plotDiscretizationToVTK(myDisc, rho);

      // Here, we explicitly invoke the two optional parameters. The first
      // allows us to name the output file, the second says which level of
      // uniform refinement should be applied to the visualization. The default
      // is 5. This is independent of M of the computation, which can be seen
      // nicely if locally constant functions are used for the discretization.
      // In the case of Maxwell, rho is visualized as a vector field on the
      // surface.
      Bembel::Vis::plotDiscretizationToVTK(myDisc, rho, "fine-grid.vtk", 7);
    }
  }
  return (0);
}
