// This file is part of Bembel, the higher order C++ boundary element library.
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
#include <Bembel/LinearForm>
#include <Bembel/IO>
#include <Eigen/Dense>
#include <iostream>

#include "Data.hpp"
#include "Error.hpp"
#include "Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;
  Bembel::IO::Stopwatch sw;
  std::complex<double> wavenumber(1., 0.);

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2, 10*10
  // points.
  Eigen::MatrixXd gridpoints = Util::makeSphereGrid(2, 10);

  // Define analytical solution using lambda function, in this case the
  // Helmholtz fundamental solution centered on 0, see Data.hpp
  const std::function<std::complex<double>(Eigen::Vector3d)> fun =
      [wavenumber](Eigen::Vector3d pt) {
        return Data::HelmholtzFundamentalSolution(pt, wavenumber,
                                                  Eigen::Vector3d(0., 0., 0.));
      };

  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  for (auto polynomial_degree : {0, 1, 2, 3}) {
    // Iterate over refinement levels
    for (auto refinement_level : {0, 1, 2, 3}) {
      sw.tic();
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level;
      // Build ansatz space
      AnsatzSpace<HelmholtzSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      DiscreteLinearForm<DirichletTrace<std::complex<double>>,
                         HelmholtzSingleLayerOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up and compute discrete operator
      DiscreteOperator<MatrixXcd, HelmholtzSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.get_linear_operator().set_wavenumber(wavenumber);
      disc_op.compute();

      // solve system
      PartialPivLU<MatrixXcd> lu;
      lu.compute(disc_op.get_discrete_operator());
      auto rho = lu.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<
          HelmholtzSingleLayerPotential<HelmholtzSingleLayerOperator>,
          HelmholtzSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.get_potential().set_wavenumber(wavenumber);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // print error
      std::cout << " time " << std::setprecision(4) << sw.toc() << "s\t\t";
      std::cout << maxPointwiseError<std::complex<double>>(pot, gridpoints, fun)
                << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;

  return 0;
}
