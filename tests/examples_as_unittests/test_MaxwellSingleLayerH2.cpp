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
#include <Bembel/IO>
#include <Bembel/LinearForm>
#include <Bembel/Maxwell>
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

  int polynomial_degree = 1;
  int refinement_level = 2;
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
  sw.tic();
  std::cout << "Degree " << polynomial_degree << " Level " << refinement_level;
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
  DiscreteOperator<H2Matrix<std::complex<double>>, MaxwellSingleLayerOperator>
      disc_op(ansatz_space);
  disc_op.get_linear_operator().set_wavenumber(wavenumber);
  disc_op.compute();

  // solve system
  GMRES<H2Matrix<std::complex<double>>, IdentityPreconditioner> gmres;
  gmres.compute(disc_op.get_discrete_operator());
  gmres.set_restart(2000);
  auto rho = gmres.solve(disc_lf.get_discrete_linear_form());

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
  double error = (pot - pot_ref).rowwise().norm().maxCoeff();
  std::cout << " time " << std::setprecision(4) << sw.toc() << "s\t\t";
  std::cout << error << std::endl;
  std::cout << std::string(60, '=') << std::endl;

  return 0;
}
