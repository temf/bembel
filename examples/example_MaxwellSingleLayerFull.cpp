
#include <iostream>

#include <Eigen/Dense>

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/LinearForm>
#include <Bembel/Maxwell>

#include "Data.hpp"
#include "Error.hpp"
#include "Grids.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  std::complex<double> wavenumber(1., 0.);

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2, 10*10
  // points.
  MatrixXd gridpoints = Util::makeSphereGrid(2, 10);

  // Define analytical solution using lambda function, in this case a dipole
  // centered on 0, see Data.hpp
  const std::function<VectorXcd(Vector3d)> fun = [wavenumber](Vector3d pt) {
    return Data::Dipole(pt, wavenumber, Vector3d(0.2, 0.2, 0.2),
                        Vector3d(0., 0.1, 0.1));
  };

  std::cout << "\n============================================================="
               "==========\n";
  // Iterate over polynomial degree.
  for (auto polynomial_degree : {1, 2, 3}) {
    // Iterate over refinement levels
    for (auto refinement_level : {0, 1, 2, 3}) {
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
      DiscreteOperator<MatrixXcd, MaxwellSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.get_linear_operator().set_wavenumber(wavenumber);
      disc_op.compute();

      // solve system
      PartialPivLU<MatrixXcd> lu;
      lu.compute(disc_op.get_discrete_operator());
      auto rho = lu.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<MaxwellSingleLayerPotential<MaxwellSingleLayerOperator>,
                        MaxwellSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.get_potential().set_wavenumber(wavenumber);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // print error
      std::cout << maxPointwiseError(pot, gridpoints, fun) << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "============================================================="
               "=========="
            << std::endl;

  return 0;
}
