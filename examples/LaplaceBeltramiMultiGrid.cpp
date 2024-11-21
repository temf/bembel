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

#include <Bembel/IO>
#include <Bembel/Identity>
#include <Bembel/LaplaceBeltrami>
#include <Bembel/LinearForm>
#include <Bembel/MultiGrid>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

#include <Bembel/src/util/Macros.hpp>
#include <Bembel/src/util/surfaceL2error.hpp>

// MG solver for the Laplace Beltrami with Neumann condition.
int mmg(const Eigen::VectorXd &constant_one,
        const Eigen::SparseMatrix<double> &S,
        const Eigen::Matrix<double, Eigen::Dynamic, 1> &L,
        const Eigen::Matrix<double, Eigen::Dynamic, 1> &f,
        Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
        const std::vector<Eigen::SparseMatrix<double>> &Ps, int lvl,
        double tol = 1e-6) {
  // stiffness matrices
  std::vector<Eigen::SparseMatrix<double>> Ss(lvl + 1);
  // the stiffness matrices on each level
  Ss[lvl] = S;
  for (int i = lvl - 1; i >= 0; --i)
    Ss[i] = Ps[i].transpose() * Ss[i + 1] * Ps[i];

  // Multi Grid Part
  int iter = 0;
  // initialize the x with zeros
  x.resize(f.size());
  x.setZero();
  Eigen::VectorXd res;
  do {
    // repetitively call multiplicativeMultiGrid()
    Bembel::MG::multiplicativeMultiGrid(Ss, f, &x, Ps, lvl);
    // minus regularization terms such that the integral of function represented
    // by x is 0. L.dot(x) is integral of x and L.dot(constant_one) is integral
    // of constant one function.
    x = x - L.dot(x) / L.dot(constant_one) * constant_one;
    res = f - S * x;
    ++iter;
  } while (res.norm() > tol);
  return iter;
}

int main() {
  using namespace Bembel;
  using namespace Eigen;
  IO::Stopwatch sw;
  std::function<double(const Vector3d &)> fun = [](const Vector3d &in) {
    // Spherical harmonics
    return sqrt(3 / BEMBEL_PI) * in(2);
  };
  std::function<double(const Vector3d &)> refsol = [](const Vector3d &in) {
    return 0.5 * sqrt(3 / BEMBEL_PI) * in(2);
  };

  std::function<double(const Vector3d &)> const_fun = [](const Vector3d &in) {
    return 1.0;
  };
  Geometry geometry("sphere.dat");
  std::cout << "\n============================================================="
               "==========\n";
  // Iterate over polynomial degree.
  for (auto polynomial_degree : {1, 2, 3, 4}) {
    std::cout << "Degree " << polynomial_degree << std::endl;
    int MAX_LVL = 7;

    std::vector<Eigen::SparseMatrix<double>> Ps;

    // Iterate over refinement levels
    std::cout << std::left << std::setw(8) << "Level" << std::left
              << std::setw(8) << "Iters" << std::left << std::setw(15)
              << "L2 norm" << std::left << std::setw(15) << "Time/sec"
              << std::endl;

    for (int refinement_level = 0; refinement_level <= MAX_LVL;
         ++refinement_level) {
      // Build ansatz space
      AnsatzSpace<LaplaceBeltramiOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);
      if (refinement_level != 0)
        Ps.push_back(MG::prolongationMatrix(ansatz_space, refinement_level));
      // Set up and compute discrete operator
      DiscreteLocalOperator<LaplaceBeltramiOperator> disc_op(ansatz_space);
      disc_op.compute();

      DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator>
          disc_lf_const(ansatz_space);
      disc_lf_const.get_linear_form().set_function(const_fun);
      disc_lf_const.compute();

      DiscreteLinearForm<DirichletTrace<double>, LaplaceBeltramiOperator>
          disc_lf(ansatz_space);
      disc_lf.get_linear_form().set_function(fun);
      disc_lf.compute();

      AnsatzSpace<MassMatrixScalarCont> temp_ansatz_space(
          geometry, refinement_level, polynomial_degree);
      DiscreteLocalOperator<MassMatrixScalarCont> disc_identity_op(
          temp_ansatz_space);
      disc_identity_op.compute();
      disc_identity_op.get_discrete_operator();

      // prepare solver for the mass matrix
      Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
      // Compute the numerical factorization
      solver.compute(disc_identity_op.get_discrete_operator());
      Eigen::VectorXd constant_one =
          solver.solve(disc_lf_const.get_discrete_linear_form());

      // solve linear system of equation with MG
      Eigen::Matrix<double, Eigen::Dynamic, 1> x;
      sw.tic();
      int iters = mmg(constant_one, disc_op.get_discrete_operator(),
                      disc_lf_const.get_discrete_linear_form(),
                      disc_lf.get_discrete_linear_form(), x, Ps,
                      refinement_level, 1e-10);
      auto difft = sw.toc();
      // print INFO
      auto err = surfaceL2error(ansatz_space, x, refsol);
      std::cout << std::left << std::setw(8) << refinement_level << std::left
                << std::setw(8) << iters << std::left << std::setw(15) << err
                << std::left << std::setw(15) << difft << std::endl;
    }
  }
  std::cout << "============================================================="
               "=========="
            << std::endl;

  return 0;
}
