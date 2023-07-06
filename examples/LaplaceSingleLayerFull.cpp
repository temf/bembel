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
#include <Bembel/Laplace>
#include <Bembel/LinearForm>
#include <Eigen/Dense>
#include <iostream>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

/**
 *  \ingroup Examples
 *  \defgroup LaplaceSingleLayerFull Laplace Single Layer Full
 *  \brief This example computes the solution of a Laplace problem.
 *
 *  Consider the following Dirichlet problem.
 *  Let \f$\Omega\f$ be a compact domain with boundary \f$\Gamma\f$.
 *  Given a function \f$g\f$ on \f$\Gamma\f$, find a function \f$u\f$ on
 *  \f$\Omega\f$ such that,
 *
 *  \f{aligned}{
 *    \Delta u &= 0\quad\mathrm{in}\, \Omega \\
 *    u &= f\quad\mathrm{on}\, \Gamma.
 *  \f}
 *
 *  The solution can be computed with the representation formulae
 *
 *  \f{aligned}{
 *    u(x) = \tilde V(\omega)(x)\quad\mathrm{in}\,\Omega,
 *  \f}
 *
 *  where
 *
 *  \f{aligned}{
 *    \tilde V(\mu)(x) = \int_\Gamma \frac{\mu(y)}{4\pi\|x - y\|}\,\mathrm{d}y
 *  \f}
 *
 *  is the Laplace single layer potential with some density \f$\mu\f$.
 *  Utilizing the trace operator \f$\gamma\f$, the single layer operator is
 *  defined by \f$V = \gamma\circ\tilde V\f$. When discretized in a conformal
 *  finite-dimensional function space \f$\mathbb{S}\f$, the variational
 *  formulation is as follows: Find \f$\omega\in\mathbb{S}\f$, such that
 *
 *  \f{aligned}{
 *    \int_\Gamma \mu \int_\Gamma V(\omega)(x)\,\mathrm{d}x = \int_\Gamma \mu
 *  f\,\mathrm{d}x\quad\forall\mu\in \mathbb{S} \\ \f}
 *
 *holds, which is realized in this example.
 */

int main() {
  using namespace Bembel;
  using namespace Eigen;
  Bembel::IO::Stopwatch sw;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Define evaluation points for potential field, a tensor product grid of
  // 7*7*7 points in [-.1,.1]^3
  MatrixXd gridpoints = Util::makeTensorProductGrid(
      VectorXd::LinSpaced(10, -.25, .25), VectorXd::LinSpaced(10, -.25, .25),
      VectorXd::LinSpaced(10, -.25, .25));

  // Define analytical solution using lambda function, in this case a harmonic
  // function, see Data.hpp
  std::function<double(Vector3d)> fun = [](Vector3d in) {
    return Data::HarmonicFunction(in);
  };

  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  for (int polynomial_degree = 0; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    VectorXd error(refinement_level_max + 1);
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      sw.tic();

      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level << "\t\t";
      // Build ansatz space
      AnsatzSpace<LaplaceSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up load vector
      DiscreteLinearForm<DirichletTrace<double>, LaplaceSingleLayerOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      // Set up and compute discrete operator
      DiscreteOperator<MatrixXd, LaplaceSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.compute();

      // solve system
      LLT<MatrixXd> llt;
      llt.compute(disc_op.get_discrete_operator());
      auto rho = llt.solve(disc_lf.get_discrete_linear_form());

      // evaluate potential
      DiscretePotential<LaplaceSingleLayerPotential<LaplaceSingleLayerOperator>,
                        LaplaceSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      // compute reference, print time, and compute error
      VectorXd pot_ref(gridpoints.rows());
      for (int i = 0; i < gridpoints.rows(); ++i)
        pot_ref(i) = fun(gridpoints.row(i));
      error(refinement_level) = (pot - pot_ref).cwiseAbs().maxCoeff();
      std::cout << " time " << std::setprecision(4) << sw.toc() << "s\t\t";
      std::cout << error(refinement_level) << std::endl;
    }

    // estimate rate of convergence and check whether it is at least 90% of the
    // expected value
    assert(
        checkRateOfConvergence(error.tail(2), 2 * polynomial_degree + 3, 0.9));

    std::cout << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;

  return 0;
}
