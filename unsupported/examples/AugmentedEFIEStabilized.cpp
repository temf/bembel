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
 *  \ingroup Examples
 *  \defgroup AugmentedEFIE Augmented Electric Field Integral Equation Method
 *
 *  \brief This example uses the A-EFIE method to solve a scattering problem.
 *
 *  This example sets up an academic example to apply the A-EFIE method to a
 *  problem where the analytical solution is known. A dipole is placed inside
 *  the unit sphere.
 *  Low frequency stability is demonstrated by applying the deflation method.
 */

template <typename t>
class logspace {
 private:
  t curvalue, base;

 public:
  logspace(t first, t base) : curvalue(first), base(base) {}

  t operator()() {
    t retval = curvalue;
    curvalue *= base;
    return retval;
  }
};

int main() {
  using namespace Bembel;
  using namespace Eigen;

  const int refinement_level = 3;
  const int polynomial_degree_max = 3;

  const double f_min = 1e-9;
  const double base = 3;
  const int N = 40;
  std::vector<double> frequencies;
  std::generate_n(std::back_inserter(frequencies), N + 1,
                  logspace<double>(f_min, base));

  std::cout << "f_min: " << f_min << std::endl;
  std::cout << "f_max: " << frequencies[frequencies.size() - 1] << std::endl;
  std::cout << "size:  " << frequencies.size() << std::endl;

  Geometry geometry("sphere.dat");

  // Define evaluation points for scattered field, sphere of radius 2
  MatrixXd gridpoints = Util::makeSphereGrid(2., 5);

  std::cout << "\n" << std::string(60, '=') << "\n";
  for (int polynomial_degree = 1; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    std::string filename = "Dipole_P" + std::to_string(polynomial_degree) +
                           "M" + std::to_string(refinement_level);

    IO::Logger<10> logger(filename + ".txt", ",");
    logger.file("f", "original_condition", "original_pw_error",
                "stabilized_condition", "stabilized_pw_error");
    std::cout << "Degree " << polynomial_degree << " Level " << refinement_level
              << std::endl;
    for (auto f = frequencies.begin(); f != frequencies.end(); ++f) {
      std::cout << std::setprecision(1) << std::scientific << *f << "Hz\t";

      // set up excitation
      double omega = 2 * BEMBEL_PI * *f;
      double wavenumber_re =
          omega * std::sqrt(Constants::mu0 * Constants::eps0);
      std::complex<double> wavenumber =
          std::complex<double>(wavenumber_re, 0.0);

      Vector3cd k(0, 0, wavenumber);
      Vector3cd E0(1, 0, 0);
      // Define analytical solution using lambda function, in this case a dipole
      // centered on 0, see Data.hpp
      const std::function<VectorXcd(Vector3d)> fun = [wavenumber](Vector3d pt) {
        return Data::Dipole(pt, wavenumber, Vector3d(0.2, 0.2, 0.2),
                            Vector3d(0., 0.1, 0.1));
      };

      // set up vector valued ansatz space
      AnsatzSpace<InductanceMatrix> ansatz_space_vector(
          geometry, refinement_level, polynomial_degree);

      // set up matrix
      AugmentedEFIE<MatrixXcd, InductanceMatrix> augmented_efie(
          ansatz_space_vector, geometry);
      augmented_efie.set_wavenumber(wavenumber);
      augmented_efie.set_omega(omega);
      augmented_efie.compute();

      // assemble the right hand side
      AugmentedEFIEExcitation<RotatedTangentialTrace<std::complex<double>>,
                              InductanceMatrix>
          excitation(ansatz_space_vector, augmented_efie.get_dofs_scalar());
      excitation.get_linear_form().set_function(fun);
      excitation.compute();

      // prepare post processing
      DiscreteLocalOperator<MassMatrixScalarDisc> mass_matrix(
          augmented_efie.get_ansatz_space_mass());
      mass_matrix.compute();
      MatrixXcd M = MatrixXcd(mass_matrix.get_discrete_operator());
      DiscreteOperator<MatrixXcd, HelmholtzSingleLayerOperator> cap_matrix(
          augmented_efie.get_ansatz_space_scalar());
      cap_matrix.get_linear_operator().set_wavenumber(wavenumber);
      cap_matrix.compute();
      MatrixXcd P = cap_matrix.get_discrete_operator() / Constants::eps0;
      PartialPivLU<MatrixXcd> lu_P;
      lu_P.compute(P);

      // solve the system, compute error and condition number
      double condition_original, condition_stabilized;
      double error_original, error_stabilized;
      for (auto simulation = 0; simulation < 2; ++simulation) {
        bool stabilized = false;
        if (simulation == 1) {
          stabilized = true;
        }

        // apply stabilization
        if (stabilized) {
          augmented_efie.stabilize();
        }

        // solve system
        PartialPivLU<MatrixXcd> lu;
        lu.compute(augmented_efie.get_system_matrix());
        VectorXcd u = lu.solve(excitation.get_excitation());
        VectorXcd j = u.head(augmented_efie.get_dofs_vector());
        VectorXcd phi = u.tail(augmented_efie.get_dofs_scalar());

        // correct stabilization in right hand side
        if (stabilized) {
          j *= 1. / (omega * std::complex<double>(0., 1.) * Constants::mu0);
        }

        // post processing
        VectorXcd q_post = lu_P.solve(-M * phi);

        MixedEFIE<InductanceMatrix, HelmholtzSingleLayerOperator> efie(
            ansatz_space_vector, augmented_efie.get_ansatz_space_scalar());

        efie.set_omega(omega);
        efie.set_wavenumber(wavenumber);
        efie.set_current(j);
        efie.set_charges(q_post);

        auto pot = efie.evaluate(gridpoints);

        // compute reference solution
        MatrixXcd pot_ref(gridpoints.rows(), 3);
        for (int i = 0; i < gridpoints.rows(); ++i)
          pot_ref.row(i) = fun(gridpoints.row(i));
        double error = (pot - pot_ref).rowwise().norm().maxCoeff();

        // compute condition number
        JacobiSVD<MatrixXcd> svd(augmented_efie.get_system_matrix());
        double cond =
            std::abs(svd.singularValues()(0) /
                     svd.singularValues()(svd.singularValues().size() - 1));

        if (stabilized) {
          std::cout << "\t\t(stabilized) Error:\t";
          error_stabilized = error;
          condition_stabilized = cond;
        } else {
          std::cout << "(original)   Error:\t";
          error_original = error;
          condition_original = cond;
        }

        std::cout << error << "\tCondition number:\t" << cond << std::endl;
      }
      logger.file(*f, condition_original, error_original, condition_stabilized,
                  error_stabilized);
    }
  }
  std::cout << std::string(60, '=') << std::endl;
  return 0;
}
