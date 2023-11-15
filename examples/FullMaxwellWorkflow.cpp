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
#include <Bembel/LinearForm>
#include <Bembel/Maxwell>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
//
#include <iostream>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;
  std::complex<double> wavenumber(2., 0.);

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2, 100
  // equispaced points.
  MatrixXd gridpoints = Util::makeSphereGrid(2., 100);

  // Define analytical solution using lambda function, in this case a dipole
  // centered on 0, see Data.hpp
  const std::function<VectorXcd(Vector3d)> fun = [wavenumber](Vector3d pt) {
    return Data::Dipole(pt, wavenumber, Vector3d(0.2, 0.2, 0.2),
                        Vector3d(0., 0.1, 0.1));
  };

  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  for (int polynomial_degree = 1; polynomial_degree < polynomial_degree_max + 1;
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
      AnsatzSpace<MaxwellSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      DiscreteLinearForm<RotatedTangentialTrace<std::complex<double>>,
                         MaxwellSingleLayerOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up and compute discrete operator
      DiscreteOperator<H2Matrix<std::complex<double>>,
                       MaxwellSingleLayerOperator>
          disc_op(ansatz_space);
      disc_op.get_linear_operator().set_wavenumber(wavenumber);
      disc_op.compute();

      // solve system
      GMRES<H2Matrix<std::complex<double>>, IdentityPreconditioner> gmres;
      gmres.compute(disc_op.get_discrete_operator());
      gmres.set_restart(2000);
      VectorXcd rho = gmres.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<MaxwellSingleLayerPotential<MaxwellSingleLayerOperator>,
                        MaxwellSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.get_potential().set_wavenumber(wavenumber);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // compute reference, print time, and compute error
      MatrixXcd pot_ref(gridpoints.rows(), 3);
      for (int i = 0; i < gridpoints.rows(); ++i)
        pot_ref.row(i) = fun(gridpoints.row(i));
      error(refinement_level) = (pot - pot_ref).rowwise().norm().maxCoeff();

      logger.both(polynomial_degree, refinement_level, error(refinement_level));

      // we only need one visualization
      if (refinement_level == 3) {
        // export geometry with surface current
        VTKSurfaceExport writer(geometry, 5);
        writer.addDataSet("SurfaceCurrent", ansatz_space, rho);
        writer.writeToFile("MaxwellSingle.vtp");

        // export point evaluations, can be visualized using glyphs in paraview
        VTKPointExport writer_points(gridpoints);
        writer_points.addComplexDataSet("Potential", pot);
        writer_points.writeToFile("MaxwellSinglePoints.vtp");
      }
    }

    std::cout << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;

  return 0;
}
