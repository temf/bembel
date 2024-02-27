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

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/Helmholtz>
#include <Bembel/IO>
#include <Bembel/Identity>
#include <Bembel/LinearForm>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  Bembel::IO::Stopwatch sw;

  int polynomial_degree_max = 2;
  int refinement_level_max = 3;
  std::complex<double> wavenumber(2., 0.);

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2, 10*10
  // points.
  MatrixXd gridpoints = Util::makeSphereGrid(2., 10);

  // Define analytical solution using lambda function, in this case the
  // Helmholtz fundamental solution centered on 0, see Data.hpp
  const std::function<std::complex<double>(Vector3d)> fun =
      [wavenumber](Vector3d pt) {
        return Data::HelmholtzFundamentalSolution(pt, wavenumber,
                                                  Vector3d(0., 0., 0.));
      };
  const std::function<Vector3cd(Vector3d)> funGrad = [wavenumber](Vector3d pt) {
    return Data::HelmholtzFundamentalSolutionGrad(pt, wavenumber,
                                                  Vector3d(0., 0., 0.));
  };

  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree
  for (int polynomial_degree = 1; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      sw.tic();
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level;
      // Build ansatz space
      AnsatzSpace<HelmholtzHypersingularOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      DiscreteLinearForm<NeumannTrace<std::complex<double>>,
                         HelmholtzHypersingularOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(funGrad);
      disc_lf.compute();

      // Set up and compute discrete operator
      DiscreteOperator<H2Matrix<std::complex<double>>,
                       HelmholtzHypersingularOperator>
          disc_op(ansatz_space);
      disc_op.get_linear_operator().set_wavenumber(wavenumber);
      disc_op.compute();

      // solve system
      GMRES<H2Matrix<std::complex<double>>, IdentityPreconditioner> gmres;
      gmres.compute(disc_op.get_discrete_operator());
      auto rho = gmres.solve(-disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<
          HelmholtzDoubleLayerPotential<HelmholtzHypersingularOperator>,
          HelmholtzHypersingularOperator>
          disc_pot(ansatz_space);
      disc_pot.get_potential().set_wavenumber(wavenumber);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // compute reference, print time, and compute error
      VectorXcd pot_ref(gridpoints.rows());
      for (int i = 0; i < gridpoints.rows(); ++i)
        pot_ref(i) = fun(gridpoints.row(i));
      error(refinement_level) = (pot - pot_ref).cwiseAbs().maxCoeff();
      std::cout << " time " << std::setprecision(4) << sw.toc() << "s\t\t";
      std::cout << error(refinement_level) << std::endl;
    }

    // estimate rate of convergence and check whether it is at least 90% of the
    // expected value
    assert(
        checkRateOfConvergence(error.tail(3), 2 * polynomial_degree + 1, 0.9));

    std::cout << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;

  return 0;
}
