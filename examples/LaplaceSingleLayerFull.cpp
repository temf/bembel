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
 * \ingroup Examples
 * \defgroup LaplaceSingleLayerOperatorFull Laplace Single Layer Operator Full
 * \brief This examples solves the Laplace problem with dense matrices.
 *
 * Let \f$\Omega\in\mathbb{R}^3\f$ be a bounded domain with Lipschitz boundary
 * \f$\Gamma = \partial\Omega\f$. We solve the Laplace problem
 *
 * \f{eqnarray*}{
 * -\Delta u &=& 0,\quad \textrm{in}\,\Omega, \\
 * u &=& g,\quad\textrm{on}\,\Gamma.
 * \f}
 *
 * The solution \f$u\f$ can be computed by a single layer potential ansatz
 *
 * \f{eqnarray*}{
 * u(\mathbf{x}) &=& \tilde{\mathcal{V}}(\rho) = \int_\Gamma
 * \frac{\rho(\mathbf{y})}{4\pi\|\mathbf{x} -
 * \mathbf{y}\|}\,\mathrm{d}\sigma(\mathbf{y}),\quad\mathbf{x}\in\Omega.
 * \f}
 *
 * Applying the Dirichlet trace operator yields the integral equation
 *
 * \f{eqnarray*}{
 * \mathcal{V}(\rho) = \gamma_0\tilde{\mathcal{V}}(\rho) &=& g
 * \f}
 *
 * on \f$\Gamma\f$.
 *
 * Let \f$\langle\cdot,\cdot\rangle\f$ denote the \f$L^2\f$-scalar product.
 * Discretizing \f$\rho^h = \in\mathbb{S}^0_{p,m}(\Gamma)\f$ with
 * discontinuous tensorized B-splines with polynomial degree \f$p\f$ and refine
 * \f$m\f$ times uniformly  yields the variational formulation: Find \f$\rho^h =
 * \in\mathbb{S}^0_{p,m}(\Gamma)\f$ such that
 *
 * \f{eqnarray*}{
 * \langle\mathcal{V}(\rho^h), \varphi\rangle &=& \langle g, \varphi\rangle,
 * \quad\forall \varphi\in\mathbb{S}^0_{p,m}(\Gamma).
 * \f}
 *
 * This system is solved by the Eigen framework. The computed density
 * \f$\rho^h\f$ is then inserted in the single layer potential ansatz to compute
 * the numerical solution of the Laplace problem.
 *
 * See \cite Dolz_2018aa for details.
 */

/**
 * \ingroup LaplaceSingleLayerOperatorFull
 * The procedure is as follows:
 */
int main() {
  using namespace Bembel;
  using namespace Eigen;
  Bembel::IO::Stopwatch sw;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  /// Load geometry from file "sphere.dat"
  Geometry geometry("sphere.dat");

  /// Define evaluation points \f$\mathbf{x}_i\f$ for potential field
  MatrixXd gridpoints = Util::makeTensorProductGrid(
      VectorXd::LinSpaced(10, -.25, .25), VectorXd::LinSpaced(10, -.25, .25),
      VectorXd::LinSpaced(10, -.25, .25));

  /// Define analytical solution \f$g\f$ using a lambda function
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

      /// Build ansatz space \f$\mathbb{S}^0_{p,m}(\Gamma)\f$ with discontinuous
      /// tensorized B-splines.
      AnsatzSpace<LaplaceSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      /// Set up linear form \f$\langle g, \varphi\rangle\f$
      DiscreteLinearForm<DirichletTrace<double>, LaplaceSingleLayerOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      /// Assemble system matrix \f$\langle\mathcal{V}(\rho^h),
      /// \varphi\rangle\f$
      DiscreteOperator<MatrixXd, LaplaceSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.compute();

      /// Solve system with Eigen
      LLT<MatrixXd> llt;
      llt.compute(disc_op.get_discrete_operator());
      auto rho = llt.solve(disc_lf.get_discrete_linear_form());

      /// Evaluate single layer potential \f$u(\mathbf{x}_i) =
      /// \tilde{\mathcal{V}}(\rho^h)\f$
      DiscretePotential<LaplaceSingleLayerPotential<LaplaceSingleLayerOperator>,
                        LaplaceSingleLayerOperator>
          disc_pot(ansatz_space);
      disc_pot.set_cauchy_data(rho);
      auto pot = disc_pot.evaluate(gridpoints);

      /// Compute error to analytical solution
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
