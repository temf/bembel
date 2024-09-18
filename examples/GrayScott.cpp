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
#include <Eigen/Dense>
#include <iostream>

#include <Bembel/src/util/surfaceL2error.hpp>

int main() {
  using namespace Bembel;
  using namespace Eigen;
  IO::Stopwatch sw;
  std::function<double(const Vector3d &)> init_u = [](const Vector3d &in) {
    // small perturbation
    return 0.6581 + 0.01 * (double)rand() / RAND_MAX - 0.005;
    // return 0.6581 + 0.005;
    // return 1.0;
  };

  std::function<double(const Vector3d &)> init_v = [](const Vector3d &in) {
    // small perturbation
    return 0.2279 + 0.01 * (double)rand() / RAND_MAX - 0.005;
    // return 0.2279 - 0.005;
  };

  std::function<double(const Vector3d &)> init_cf = [](const Vector3d &in) {
    return 1.0;
  };

  std::function<double(const double &u_val, const double &v_val)> reaction =
      [](const double &u_val, const double &v_val) {
        return u_val * v_val * v_val;
      };

  Eigen::VectorXd u_init;
  Eigen::VectorXd v_init;

  Geometry geometry("sphere.dat");
  std::cout << "\n" << std::string(60, '=') << "\n";
  // Iterate over polynomial degree.
  unsigned int polynomial_degree = 2;
  // Iterate over refinement levels
  unsigned int refinement_level = 4;
  std::cout << "Degree " << polynomial_degree << " Level " << refinement_level
            << "\n";
  sw.tic();
  // Build ansatz space
  AnsatzSpace<LaplaceBeltramiOperator> ansatz_space_lb(
      geometry, refinement_level, polynomial_degree);

  AnsatzSpace<MassMatrixScalarCont> ansatz_space_mass(
      geometry, refinement_level, polynomial_degree);
  // Gray Scott Model
  // u_t = r_u * \laplacian u -uv^2 + f(1-u)
  // v_y = r_v * \laplacian v +uv^2 - (f+k)v
  double delta_t = 0.1;
  double diffusion_rate_u = 0.01;
  double diffusion_rate_v =
      0.001;  // for sake of spots on the surface, please set the
              // diffusion_rate_v to 0.0005 or less
  double f = 0.1;
  double k = 0.05;
  // Build discrete operator for Laplace Beltrami and Identity
  DiscreteLocalOperator<LaplaceBeltramiOperator> disc_op_lb(ansatz_space_lb);
  disc_op_lb.compute();
  DiscreteLocalOperator<MassMatrixScalarCont> disc_op_mass(ansatz_space_mass);
  disc_op_mass.compute();
  // Initial solutions of u and v
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  solver.compute(disc_op_mass.get_discrete_operator());
  DiscreteLinearForm<DirichletTrace<double>, MassMatrixScalarCont> disc_lf(
      ansatz_space_mass);
  disc_lf.get_linear_form().set_function(init_u);
  disc_lf.compute();
  u_init = solver.solve(disc_lf.get_discrete_linear_form());
  disc_lf.get_linear_form().set_function(init_v);
  disc_lf.compute();
  v_init = solver.solve(disc_lf.get_discrete_linear_form());
  disc_lf.get_linear_form().set_function(init_cf);
  disc_lf.compute();
  // Left hand side
  const Eigen::SparseMatrix<double> &stiff_mat =
      disc_op_lb.get_discrete_operator();
  const Eigen::SparseMatrix<double> &mass_mat =
      disc_op_mass.get_discrete_operator();

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver_u;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver_v;
  solver_u.compute((1 + delta_t * f) * mass_mat +
                   delta_t * diffusion_rate_u * stiff_mat);
  solver_v.compute((1 + (f + k) * delta_t) * mass_mat +
                   delta_t * diffusion_rate_v * stiff_mat);
  sw.tic();
  Eigen::VectorXd u = u_init;
  Eigen::VectorXd v = v_init;
  // Solve the solutions of u and v at the end of the time
  for (int i = 0; i < 20000; ++i) {
    // Right hand side
    Eigen::VectorXd reaction_term;
    reactionLinearFrom(ansatz_space_lb, u, v, reaction, &reaction_term);
    Eigen::Matrix<double, Eigen::Dynamic, 1> rhs_u =
        (delta_t * f) * disc_lf.get_discrete_linear_form() +
        disc_op_mass.get_discrete_operator() * u - delta_t * reaction_term;

    Eigen::Matrix<double, Eigen::Dynamic, 1> rhs_v =
        disc_op_mass.get_discrete_operator() * v + delta_t * reaction_term;
    u = solver_u.solve(rhs_u);
    v = solver_v.solve(rhs_v);
  }
  auto difft = sw.toc();
  std::cout << "RD is solved in " << difft << " sec.\n";

  // Visualization
  FunctionEvaluator<MassMatrixScalarCont> evaluator(ansatz_space_mass);
  evaluator.set_function(u);
  VTKSurfaceExport writer(geometry, 6);
  std::function<double(int, const Eigen::Vector2d &)> fun1 =
      [&](int patch_number, const Eigen::Vector2d &reference_domain_point) {
        auto retval =
            evaluator.evaluateOnPatch(patch_number, reference_domain_point);
        return double(retval(0, 0));
      };
  writer.addDataSet("u", fun1);
  evaluator.set_function(v);
  writer.addDataSet("v", fun1);
  writer.writeToFile("sol.vtp");
  std::cout << "\n" << std::string(60, '=') << "\n";

  return 0;
}
