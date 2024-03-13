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

// [main]
int main() {
  using namespace Bembel;
  using namespace Eigen;

  // import geometry
  Geometry geometry("external_geometry_file.dat");

  // define ansatz space
  AnsatzSpace<Operator> ansatz_space(geometry, refinement_level,
                                     polynomial_degree);

  // define linear form
  DiscreteLinearForm<DirichletTrace<double>, Operator> disc_lf(ansatz_space);
  disc_lf.get_linear_form().set_function(fun);
  disc_lf.compute();

  // define linear operator
  DiscreteOperator<H2Matrix<double>, Operator> disc_op(ansatz_space);
  disc_op.get_linear_operator().set_wavenumber(wavenumber);
  disc_op.compute();

  // iterative solution using Eigen -> rho

  // evaluate potential
  DiscretePotential<Potential<Operator>, Operator> disc_pot(ansatz_space);
  disc_pot.set_cauchy_data(rho);
  auto pot = disc_pot.evaluate(gridpoints);

  return 0;
}
// [main]
