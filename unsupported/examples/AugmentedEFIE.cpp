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
#include <Bembel/IO>
#include <unsupported/Bembel/AugmentedEFIE>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

/**
 * \ingroup Examples
 * \defgroup AugmentedEFIEDipole Augmented Electric Field Integral Equation
 *
 * \brief This example uses the A-EFIE method to solve a scattering problem.
 *
 * This example sets up an academic example to apply the A-EFIE method to a
 * problem where the analytical solution is known. A dipole is placed inside
 * the unit sphere.
 * The scattered field in the exterior is computed with the A-EFIE formulation.
 * This is a mixed problem which is discretized with higher-order B-splines
 * whereas the geometry is represented with NURBS functions.
 */

int main() {
  using namespace Bembel;
  using namespace Eigen;

  const int refinement_level_max = 4;
  const int polynomial_degree_max = 3;

  // set up excitation
  double wavenumber_re =
      2 * BEMBEL_PI * 30e6 * std::sqrt(Constants::mu0 * Constants::eps0);
  std::complex<double> wavenumber = std::complex<double>(wavenumber_re, 0.0);

  Vector3cd k(0, 0, wavenumber);
  Vector3cd E0(1, 0, 0);
  // Define analytical solution using lambda function, in this case a dipole
  // centered on 0, see Data.hpp
  const std::function<VectorXcd(Vector3d)> fun = [wavenumber](Vector3d pt) {
    return Data::Dipole(pt, wavenumber, Vector3d(0.2, 0.2, 0.2),
                        Vector3d(0., 0.1, 0.1));
  };

  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2
  MatrixXd gridpoints = Util::makeSphereGrid(2., 5);

  std::cout << "\n" << std::string(60, '=') << "\n";
  for (int polynomial_degree = 1; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);

    IO::Logger<10> logger(
        "Dipole_P" + std::to_string(polynomial_degree) + ".txt", ",");
    logger.file("refinement_level", "dofs", "error");
    for (int refinement_level = 0; refinement_level <= refinement_level_max;
         ++refinement_level) {
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level << "\t\t";

      // set up vector valued ansatz space
      AnsatzSpace<InductanceMatrix> ansatz_space_vector(
          geometry, refinement_level, polynomial_degree);

      // set up matrix
      AugmentedEFIE<MatrixXcd, InductanceMatrix> augmented_efie(
          ansatz_space_vector, geometry);
      augmented_efie.set_wavenumber(wavenumber);
      augmented_efie.compute();

      // assemble the right hand side
      AugmentedEFIEExcitation<RotatedTangentialTrace<std::complex<double>>,
                              InductanceMatrix>
          excitation(ansatz_space_vector, augmented_efie.get_dofs_scalar());
      excitation.get_linear_form().set_function(fun);
      excitation.compute();

      // solve system
      PartialPivLU<MatrixXcd> lu;
      lu.compute(augmented_efie.get_system_matrix());
      VectorXcd u = lu.solve(excitation.get_excitation());
      VectorXcd j = u.head(augmented_efie.get_dofs_vector());

      // evaluate potential
      DiscretePotential<EFIE<InductanceMatrix>, InductanceMatrix> disc_pot(
          ansatz_space_vector);
      disc_pot.get_potential().set_wavenumber(wavenumber);
      disc_pot.set_cauchy_data(j);
      auto pot = disc_pot.evaluate(gridpoints);

      // compute reference, print time, and compute error
      MatrixXcd pot_ref(gridpoints.rows(), 3);
      for (int i = 0; i < gridpoints.rows(); ++i)
        pot_ref.row(i) = fun(gridpoints.row(i));
      error(refinement_level) = (pot - pot_ref).rowwise().norm().maxCoeff();
      const int dofs_total =
          augmented_efie.get_dofs_scalar() + augmented_efie.get_dofs_vector();
      logger.file(refinement_level, dofs_total, error(refinement_level));

      std::cout << error(refinement_level) << std::endl;
    }
    assert(checkRateOfConvergence(error.tail(2), 2 * polynomial_degree, 0.9));
    std::cout << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;
  return 0;
}
